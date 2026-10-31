"""\
=====================================================================
LAMMPS data files --- :mod:`polyphys.simulate.data`
=====================================================================

The :mod:`~polyphys.simulate.data` module reads and writes LAMMPS data
files, the input of the ``read_data`` command, for molecular systems of
atoms, bonds, and angles. A data file built here serves two ends: a start
configuration that LAMMPS reads and runs, and a topology that analysis tools
such as MDAnalysis read with a dump.

.. versionadded:: 0.6.0

Contents
========

:class:`DataFile` holds the contents of one data file as numpy arrays, its
box as a :class:`~polyphys.simulate.core.Box`, and checks them.
:func:`read_data` reads a file into it, raising :class:`BrokenDataError`, a
:class:`~polyphys.simulate.core.LammpsFormatError`, for a file that does not
follow the format; :func:`write_data` writes one out. :func:`with_frame`
takes the positions and the box of a dump frame
(:class:`~polyphys.simulate.dump.DumpFrame`) and the rest from a template, to
give a dump its matching topology.

Format
======

A data file is a header and a body. The first line is a title and is
skipped. The header lines give counts and the box, in any order::

    8 atoms
    1 atom types
    8 bonds
    1 bond types
    0.0 6.0 xlo xhi
    0.0 6.0 ylo yhi
    0.0 6.0 zlo zhi
    1.0 0.5 -0.5 xy xz yz

The body is a sequence of sections, each a keyword line followed by one line
per entry: ``Masses`` (``type mass``), ``Atoms`` (columns by atom style, see
:class:`DataFile`), ``Velocities`` (``atom-ID vx vy vz``), ``Bonds``
(``ID type atom1 atom2``), and ``Angles`` (``ID type atom1 atom2 atom3``).
``#`` starts a comment and blank lines are skipped. ``Atoms # bond`` names the
atom style in a comment. A file whose name ends in ``.gz`` is read through
gzip. Atom lines may end with the three image flags, on
all lines or on none; LAMMPS sets missing image flags and velocities to zero.

Supported atom styles: ``atomic``, ``bond``, ``angle``, ``molecular``, and
``full``. Sections with names ending in ``Coeffs`` (force-field coefficients)
are skipped when read and never written; the input script sets them, as
``write_data nocoeff`` assumes. Dihedrals, impropers, type labels, and
general triclinic boxes are not supported.

The module converts no units; values are in the ``units`` of the run.

Complexity
==========

:func:`read_data` and :func:`write_data` take :math:`O(n_a + n_b + n_t)`
time and memory for :math:`n_a` atoms, :math:`n_b` bonds, and :math:`n_t`
angles; :func:`with_frame` takes :math:`O(n_a \\log n_a)` time to match atom
IDs.

References
==========

The format follows the LAMMPS documentation of `read_data
<https://docs.lammps.org/read_data.html>`_ and `write_data
<https://docs.lammps.org/write_data.html>`_.

Dependencies
============
- Python standard library: `dataclasses`, `logging`, `pathlib`
- Third-party: `numpy`
- PolyPhys: :mod:`polyphys.simulate.core`, :mod:`polyphys.simulate.dump`
"""

import dataclasses
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import numpy.typing as npt

from polyphys.simulate.core import (
    Box,
    FloatArray,
    IntArray,
    LammpsFormatError,
    _frozen,
    _open_text,
)
from polyphys.simulate.dump import DumpFrame

__all__ = [
    "BrokenDataError",
    "DataFile",
    "read_data",
    "with_frame",
    "write_data",
]

logger = logging.getLogger(__name__)


class _Style(NamedTuple):
    """The per-atom columns and the topology an atom style allows."""

    molecule: bool
    charge: bool
    bonds: bool
    angles: bool


# Atoms-section columns per atom style (LAMMPS read_data, "Atoms section").
_STYLES = {
    "atomic": _Style(molecule=False, charge=False, bonds=False, angles=False),
    "bond": _Style(molecule=True, charge=False, bonds=True, angles=False),
    "angle": _Style(molecule=True, charge=False, bonds=True, angles=True),
    "molecular": _Style(molecule=True, charge=False, bonds=True, angles=True),
    "full": _Style(molecule=True, charge=True, bonds=True, angles=True),
}
_COUNTS = {
    "atoms": "atoms", "bonds": "bonds", "angles": "angles",
    "atom types": "atom_types", "bond types": "bond_types",
    "angle types": "angle_types",
}
_UNSUPPORTED_COUNTS = (
    "dihedrals", "impropers", "dihedral types", "improper types",
    "ellipsoids", "lines", "triangles", "bodies",
)
_AXES = (("xlo", "xhi"), ("ylo", "yhi"), ("zlo", "zhi"))
_DEFAULT_TITLE = "LAMMPS data file written by polyphys.simulate.data"


class BrokenDataError(LammpsFormatError):
    """A LAMMPS data file does not follow the data-file format.

    .. versionadded:: 0.6.0

    Examples
    --------
    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "empty.data"
    ...     _ = path.write_text("title\\n\\n2 atoms\\n")
    ...     try:
    ...         read_data(path)
    ...     except BrokenDataError as err:
    ...         print(str(err).split(": ", 1)[1])
    the file has a header but no body
    """


def _empty(width: int) -> IntArray:
    """Return an empty ``(0, width)`` int64 array."""
    return np.empty((0, width), dtype=np.int64)


@dataclass(frozen=True, eq=False)
class DataFile:
    """The contents of one LAMMPS data file.

    The arrays are copied, made read-only, and checked on construction; use
    :func:`dataclasses.replace` to change a field.

    .. versionadded:: 0.6.0

    Attributes
    ----------
    atom_style : {'atomic', 'bond', 'angle', 'molecular', 'full'}
        The LAMMPS atom style. It sets the Atoms columns: ``atom-ID
        atom-type x y z`` for ``atomic``; ``atom-ID molecule-ID atom-type x y
        z`` for ``bond``, ``angle``, and ``molecular``; and ``atom-ID
        molecule-ID atom-type q x y z`` for ``full``. ``atomic`` has no bonds,
        ``bond`` no angles.
    box : Box
        The simulation box. A triclinic box writes its ``xy xz yz`` line;
        the boundary flags are not part of a data file.
    ids : numpy.ndarray
        Shape ``(n_atoms,)``, int64: unique atom IDs, from 1.
    types : numpy.ndarray
        Shape ``(n_atoms,)``, int64: atom types, from 1 to `n_atom_types`.
    positions : numpy.ndarray
        Shape ``(n_atoms, 3)``, float64: coordinates in the length unit of
        the run. The writer does not check that they lie in the box: LAMMPS
        remaps an atom outside it across a periodic boundary and loses one
        outside a non-periodic boundary.
    n_atom_types : int
        Number of atom types, 1 or more; a type may be unused.
    masses : numpy.ndarray or None
        Shape ``(n_atom_types,)``: the mass of type :math:`k` at index
        :math:`k - 1`, or None to write no ``Masses`` section and leave the
        masses to the input script.
    mol_ids : numpy.ndarray or None
        Shape ``(n_atoms,)``, int64: molecule IDs, required by every style
        but ``atomic`` and None for it.
    charges : numpy.ndarray or None
        Shape ``(n_atoms,)``: charges, required by ``full`` and None for the
        others.
    image_flags : numpy.ndarray or None
        Shape ``(n_atoms, 3)``, int64: image flags, or None to write none
        (LAMMPS then sets them to zero).
    velocities : numpy.ndarray or None
        Shape ``(n_atoms, 3)``: velocities, or None to write no
        ``Velocities`` section (LAMMPS then sets them to zero).
    bonds : numpy.ndarray
        Shape ``(n_bonds, 3)``, int64: ``type, atom1, atom2`` per bond, the
        atoms by ID. Bond IDs are not kept: the writer numbers bonds from 1
        in row order. A ring is a chain with one more bond.
    angles : numpy.ndarray
        Shape ``(n_angles, 4)``, int64: ``type, atom1, atom2, atom3`` per
        angle, with ``atom2`` the vertex.
    n_bond_types, n_angle_types : int
        Numbers of bond and angle types; 1 or more when there are bonds or
        angles.

    Raises
    ------
    TypeError
        If `box` is not a :class:`~polyphys.simulate.core.Box`.
    ValueError
        If the atom style is unknown, an array has the wrong shape or
        non-integer values where integers are needed, an atom ID repeats or
        is below 1, a type is out of range, a bond or angle names an atom ID
        that does not exist, or a field is given that the style does not
        allow or omitted where it needs one.

    Examples
    --------
    A ring of three beads with a closing bond:

    >>> ring = DataFile(
    ...     atom_style="bond", box=Box(lo=np.zeros(3), hi=np.full(3, 10.0)),
    ...     ids=np.arange(1, 4), types=np.ones(3, dtype=int),
    ...     positions=np.array([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0],
    ...                         [5.5, 5.9, 5.0]]),
    ...     n_atom_types=1, mol_ids=np.ones(3, dtype=int),
    ...     bonds=np.array([[1, 1, 2], [1, 2, 3], [1, 3, 1]]), n_bond_types=1)
    >>> ring.n_atoms, ring.n_bonds, ring.n_angles
    (3, 3, 0)
    """

    atom_style: str
    box: Box
    ids: IntArray
    types: IntArray
    positions: FloatArray
    n_atom_types: int
    masses: FloatArray | None = None
    mol_ids: IntArray | None = None
    charges: FloatArray | None = None
    image_flags: IntArray | None = None
    velocities: FloatArray | None = None
    bonds: IntArray = field(default_factory=lambda: _empty(3))
    angles: IntArray = field(default_factory=lambda: _empty(4))
    n_bond_types: int = 0
    n_angle_types: int = 0

    def __post_init__(self) -> None:
        if self.atom_style not in _STYLES:
            raise ValueError(
                f"'atom_style' must be one of {tuple(_STYLES)}, "
                f"got {self.atom_style!r}"
            )
        if not isinstance(self.box, Box):
            raise TypeError(f"'box' must be a Box, got {type(self.box)}")
        style = _STYLES[self.atom_style]
        ids = self._set_int("ids", (-1,))
        n_atoms = ids.shape[0]
        self._set_float("positions", (n_atoms, 3))
        types = self._set_int("types", (n_atoms,))
        self._set_optional_float("velocities", (n_atoms, 3))
        self._set_optional_int("image_flags", (n_atoms, 3))
        self._set_needed("mol_ids", style.molecule, (n_atoms,), integer=True)
        self._set_needed("charges", style.charge, (n_atoms,), integer=False)
        bonds = self._set_int("bonds", (-1, 3))
        angles = self._set_int("angles", (-1, 4))
        if np.unique(ids).size != n_atoms or np.any(ids < 1):
            raise ValueError("atom IDs must be unique and 1 or more")
        _check_types("atom", types, self.n_atom_types, minimum=1)
        self._set_optional_float("masses", (self.n_atom_types,))
        if bonds.shape[0] and not style.bonds:
            raise ValueError(f"atom style {self.atom_style!r} has no bonds")
        if angles.shape[0] and not style.angles:
            raise ValueError(f"atom style {self.atom_style!r} has no angles")
        _check_types("bond", bonds[:, 0], self.n_bond_types,
                     minimum=1 if bonds.shape[0] else 0)
        _check_types("angle", angles[:, 0], self.n_angle_types,
                     minimum=1 if angles.shape[0] else 0)
        for name, members in (("bond", bonds), ("angle", angles)):
            missing = np.setdiff1d(members[:, 1:], ids)
            if missing.size:
                raise ValueError(
                    f"a {name} names atom IDs that do not exist: "
                    f"{missing.tolist()}"
                )

    @property
    def n_atoms(self) -> int:
        """The number of atoms."""
        return int(self.ids.shape[0])

    @property
    def n_bonds(self) -> int:
        """The number of bonds."""
        return int(self.bonds.shape[0])

    @property
    def n_angles(self) -> int:
        """The number of angles."""
        return int(self.angles.shape[0])

    def unwrapped_positions(self) -> FloatArray:
        """Return the positions unwrapped by the image flags.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_atoms, 3)``, float64: the positions moved by the image
            flags along the box edge vectors, by
            :meth:`~polyphys.simulate.core.Box.unwrap`; a copy of the
            positions when there are no image flags, which LAMMPS reads as
            zero. Time and memory are :math:`O(n_a)`.

        Examples
        --------
        >>> pair = DataFile(
        ...     atom_style="atomic",
        ...     box=Box(lo=np.zeros(3), hi=np.full(3, 10.0)),
        ...     ids=np.array([1, 2]), types=np.array([1, 1]),
        ...     positions=np.array([[9.5, 5.0, 5.0], [0.5, 5.0, 5.0]]),
        ...     n_atom_types=1, image_flags=np.array([[-1, 0, 0], [0, 0, 0]]))
        >>> pair.unwrapped_positions()
        array([[-0.5,  5. ,  5. ],
               [ 0.5,  5. ,  5. ]])
        """
        if self.image_flags is None:
            return self.positions.copy()
        return self.box.unwrap(self.positions, self.image_flags)

    def _set_float(self, name: str, shape: tuple[int, ...]) -> FloatArray:
        value = _frozen(getattr(self, name), np.float64)
        _check_shape(name, value, shape)
        object.__setattr__(self, name, value)
        return value

    def _set_int(self, name: str, shape: tuple[int, ...]) -> IntArray:
        value = _integers(name, getattr(self, name))
        if value.size == 0 and len(shape) == 2 and shape[1] > 0:
            value = value.reshape(0, shape[1])  # no bonds or no angles
        _check_shape(name, value, shape)
        value = _frozen(value, np.int64)
        object.__setattr__(self, name, value)
        return value

    def _set_optional_float(self, name: str, shape: tuple[int, ...]) -> None:
        if getattr(self, name) is not None:
            self._set_float(name, shape)

    def _set_optional_int(self, name: str, shape: tuple[int, ...]) -> None:
        if getattr(self, name) is not None:
            self._set_int(name, shape)

    def _set_needed(
        self, name: str, needed: bool, shape: tuple[int, ...], integer: bool
    ) -> None:
        given = getattr(self, name) is not None
        if needed and not given:
            raise ValueError(f"atom style {self.atom_style!r} needs '{name}'")
        if given and not needed:
            raise ValueError(
                f"atom style {self.atom_style!r} has no '{name}'"
            )
        if given and integer:
            self._set_int(name, shape)
        elif given:
            self._set_float(name, shape)


def _integers(name: str, value: Any) -> IntArray:
    """Return `value` as an int64 array, refusing non-integer values."""
    array = np.array(value)
    if array.size and not np.issubdtype(array.dtype, np.integer):
        if not np.issubdtype(array.dtype, np.number) or not np.all(
            array == np.round(array)
        ):
            raise ValueError(f"'{name}' must hold integers")
    return array.astype(np.int64)


def _check_shape(
    name: str, value: npt.NDArray[Any], shape: tuple[int, ...]
) -> None:
    """Check `value` against `shape`, where -1 matches any length."""
    if value.ndim != len(shape) or any(
        want not in (-1, got) for want, got in zip(shape, value.shape)
    ):
        raise ValueError(
            f"'{name}' must have shape {shape}, got {value.shape}"
        )


def _check_types(
    kind: str, types: IntArray, n_types: int, minimum: int
) -> None:
    """Check that every type lies from 1 to `n_types`."""
    if n_types < minimum:
        raise ValueError(f"'n_{kind}_types' must be {minimum} or more, "
                         f"got {n_types}")
    if types.size and (types.min() < 1 or types.max() > n_types):
        raise ValueError(f"{kind} types must be from 1 to {n_types}")


def write_data(
    data: DataFile, path: str | Path, *, title: str | None = None
) -> Path:
    """Write a LAMMPS data file.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    data : DataFile
        The contents to write.
    path : str or pathlib.Path
        The file to write; an existing file is overwritten.
    title : str, optional
        The first line, which LAMMPS skips; one line, defaulting to a line
        that names this module.

    Returns
    -------
    pathlib.Path
        `path`.

    Raises
    ------
    ValueError
        If `title` has a line break.

    Notes
    -----
    The header gives ``atoms`` and ``atom types``; ``bonds`` and ``bond
    types`` for every style with bonds; ``angles`` and ``angle types`` for
    every style with angles; the box; and ``xy xz yz`` when its box is
    triclinic. Then come ``Masses`` if given, ``Atoms # <style>`` with
    image flags if given, ``Velocities`` if given, and ``Bonds`` and
    ``Angles`` if there are any, numbered from 1. Floats are written as the
    shortest decimal that reads back to the same float64, so
    ``read_data(write_data(data, path))`` returns `data` unchanged. No
    coefficient section is written. Time and memory are :math:`O(n_a + n_b +
    n_t)`.

    Examples
    --------
    >>> import tempfile
    >>> pair = DataFile(
    ...     atom_style="atomic", box=Box(lo=np.zeros(3), hi=np.full(3, 5.0)),
    ...     ids=np.array([1, 2]), types=np.array([1, 2]),
    ...     positions=np.array([[2.0, 2.5, 2.5], [3.1, 2.5, 2.5]]),
    ...     n_atom_types=2)
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = write_data(pair, Path(tmp) / "pair.data", title="Two beads")
    ...     print(path.read_text())
    Two beads
    <BLANKLINE>
    2 atoms
    2 atom types
    <BLANKLINE>
    0.0 5.0 xlo xhi
    0.0 5.0 ylo yhi
    0.0 5.0 zlo zhi
    <BLANKLINE>
    Atoms # atomic
    <BLANKLINE>
    1 1 2.0 2.5 2.5
    2 2 3.1 2.5 2.5
    <BLANKLINE>
    """
    title = _DEFAULT_TITLE if title is None else title
    if "\n" in title or "\r" in title:
        raise ValueError("'title' must be one line")
    style = _STYLES[data.atom_style]
    lines = [title, "", f"{data.n_atoms} atoms",
             f"{data.n_atom_types} atom types"]
    if style.bonds:
        lines += [f"{data.n_bonds} bonds", f"{data.n_bond_types} bond types"]
    if style.angles:
        lines += [f"{data.n_angles} angles",
                  f"{data.n_angle_types} angle types"]
    lines.append("")
    for (lo_name, hi_name), lo, hi in zip(_AXES, data.box.lo.tolist(),
                                          data.box.hi.tolist()):
        lines.append(f"{lo!r} {hi!r} {lo_name} {hi_name}")
    if data.box.triclinic:
        lines.append(" ".join(map(repr, data.box.tilt.tolist()))
                     + " xy xz yz")
    if data.masses is not None:
        lines += ["", "Masses", ""]
        lines += [f"{k} {mass!r}"
                  for k, mass in enumerate(data.masses.tolist(), start=1)]
    lines += ["", f"Atoms # {data.atom_style}", ""]
    lines += _atom_lines(data)
    if data.velocities is not None:
        lines += ["", "Velocities", ""]
        lines += [f"{atom} " + " ".join(map(repr, v))
                  for atom, v in zip(data.ids.tolist(),
                                     data.velocities.tolist())]
    for name, rows in (("Bonds", data.bonds), ("Angles", data.angles)):
        if rows.shape[0]:
            lines += ["", name, ""]
            lines += [f"{k} " + " ".join(map(str, row))
                      for k, row in enumerate(rows.tolist(), start=1)]
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _atom_lines(data: DataFile) -> list[str]:
    """Return the Atoms-section lines in the column order of the style."""
    columns: list[list[str]] = [[str(atom)] for atom in data.ids.tolist()]
    if data.mol_ids is not None:
        for row, mol in zip(columns, data.mol_ids.tolist()):
            row.append(str(mol))
    for row, atom_type in zip(columns, data.types.tolist()):
        row.append(str(atom_type))
    if data.charges is not None:
        for row, charge in zip(columns, data.charges.tolist()):
            row.append(repr(charge))
    for row, xyz in zip(columns, data.positions.tolist()):
        row.extend(map(repr, xyz))
    if data.image_flags is not None:
        for row, image in zip(columns, data.image_flags.tolist()):
            row.extend(map(str, image))
    return [" ".join(row) for row in columns]


class _Section(NamedTuple):
    """One body section: its name, its comment, and its numbered rows."""

    name: str
    comment: str
    line: int
    rows: list[tuple[int, list[str]]]


def read_data(path: str | Path, atom_style: str | None = None) -> DataFile:
    """Read a LAMMPS data file.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    path : str or pathlib.Path
        The data file, read as UTF-8; a name ending in ``.gz`` is read
        through gzip.
    atom_style : str, optional
        The atom style, one of those of :class:`DataFile`. Needed when the
        ``Atoms`` line has no style comment such as ``Atoms # bond``; when
        both are given they must agree.

    Returns
    -------
    DataFile
        Atoms sorted by ID; bonds and angles sorted by their IDs, which are
        then dropped.

    Raises
    ------
    BrokenDataError
        If the file does not follow the format, uses a feature this module
        does not support (see the module documentation), or gives no atom
        style, or two that disagree. The message names the file and, where
        it can, the line, from 1.

    Notes
    -----
    Sections whose names end in ``Coeffs`` are skipped and logged at level
    DEBUG. Time and memory are :math:`O(L)` in the length :math:`L` of the
    file.

    Examples
    --------
    >>> import tempfile
    >>> text = '''Two beads and one bond
    ...
    ... 2 atoms
    ... 1 atom types
    ... 1 bonds
    ... 1 bond types
    ... 0 5 xlo xhi
    ... 0 5 ylo yhi
    ... 0 5 zlo zhi
    ...
    ... Bond Coeffs # fene
    ...
    ... 1 30.0 1.5 1.0 1.0
    ...
    ... Atoms # bond
    ...
    ... 2 1 1 3.0 2.5 2.5 0 0 0
    ... 1 1 1 2.0 2.5 2.5 0 0 0
    ...
    ... Bonds
    ...
    ... 1 1 1 2
    ... '''
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "pair.data"
    ...     _ = path.write_text(text)
    ...     data = read_data(path)
    >>> data.atom_style, data.ids, data.bonds
    ('bond', array([1, 2]), array([[1, 1, 2]]))
    >>> data.positions[:, 0]
    array([2., 3.])
    """
    path = Path(path)
    with _open_text(path) as stream:
        lines = stream.read().splitlines()
    header, body = _read_header(lines, path)
    sections = _read_sections(lines, body, path)
    atoms = next((s for s in sections if s.name == "Atoms"), None)
    if atoms is None:
        raise _broken(path, None, "no Atoms section")
    style = _atom_style(atoms, atom_style, path)
    try:
        fields = _read_body(sections, style, header, path)
        return DataFile(
            box=header["box"], n_atom_types=header["atom_types"],
            n_bond_types=header["bond_types"],
            n_angle_types=header["angle_types"], **fields,
        )
    except BrokenDataError:
        raise
    except ValueError as err:
        raise _broken(path, None, str(err)) from None


def _broken(path: Path, line: int | None, message: str) -> BrokenDataError:
    """Return an error naming the file and, when known, the line."""
    where = f"{path}" if line is None else f"{path}, line {line}"
    return BrokenDataError(f"{where}: {message}")


def _read_body(
    sections: list["_Section"], style: str, header: dict[str, Any],
    path: Path,
) -> dict[str, Any]:
    """Read the sections into the fields of :class:`DataFile`."""
    fields: dict[str, Any] = {"atom_style": style}
    for section in sections:
        if section.name == "Atoms":
            fields.update(_read_atoms(section, style, header, path))
        elif section.name == "Velocities":
            fields["velocities"] = _read_velocities(section, header, path)
        elif section.name == "Masses":
            fields["masses"] = _read_masses(section, header, path)
        elif section.name in ("Bonds", "Angles"):
            width = 4 if section.name == "Bonds" else 5
            count = header["bonds" if width == 4 else "angles"]
            rows = _table(section, width, count, path)
            order = np.argsort(rows[:, 0], kind="stable")
            fields[section.name.lower()] = _integers(
                section.name, rows[order, 1:]
            )
        elif section.name.endswith("Coeffs"):
            logger.debug("%s: skipped section %r (line %d)",
                         path, section.name, section.line)
        else:
            raise _broken(path, section.line,
                          f"section {section.name!r} is not supported")
    if "velocities" in fields:
        fields["velocities"] = _by_atom(fields.pop("velocities"),
                                        fields["ids"], "Velocities", path)
    if header["bonds"] and "bonds" not in fields:
        raise _broken(path, None, f"{header['bonds']} bonds but no Bonds")
    if header["angles"] and "angles" not in fields:
        raise _broken(path, None, f"{header['angles']} angles but no Angles")
    return fields


def _content(line: str) -> tuple[str, str]:
    """Split a line into its content and its comment, both stripped."""
    content, _, comment = line.partition("#")
    return content.strip(), comment.strip()


def _read_header(lines: list[str], path: Path) -> tuple[dict[str, Any], int]:
    """Read the header; return its values and the index of the body."""
    header: dict[str, Any] = {key: 0 for key in _COUNTS.values()}
    lo, hi = np.zeros(3), np.zeros(3)
    tilt: FloatArray | None = None
    seen: set[str] = set()
    for index in range(1, len(lines)):
        content, _ = _content(lines[index])
        if not content:
            continue
        words = content.split()
        keyword = " ".join(words[1:])
        try:
            if words[-3:] == ["xy", "xz", "yz"] and len(words) == 6:
                tilt = np.array(words[:3], dtype=np.float64)
            elif len(words) == 4 and tuple(words[2:]) in _AXES:
                axis = _AXES.index((words[2], words[3]))
                lo[axis] = float(words[0])
                hi[axis] = float(words[1])
                seen.add(words[2])
            elif keyword in _COUNTS:
                header[_COUNTS[keyword]] = int(words[0])
            elif keyword in _UNSUPPORTED_COUNTS:
                if int(words[0]):
                    raise ValueError(f"{keyword} are not supported")
            elif len(words) > 2 and words[1] == "extra":
                logger.debug("%s: skipped header line %r", path, content)
            elif words[-2:] == ["abc", "origin"] or words[-1] in (
                "avec", "bvec", "cvec"
            ):
                raise ValueError("general triclinic boxes are not supported")
            else:
                break
        except ValueError as err:
            raise _broken(path, index + 1, str(err)) from None
    else:
        raise _broken(path, None, "the file has a header but no body")
    missing = [name for name, _ in _AXES if name not in seen]
    if missing:
        raise _broken(path, None, f"the header gives no {missing} bounds")
    try:
        header["box"] = Box(lo=lo, hi=hi,
                            tilt=np.zeros(3) if tilt is None else tilt,
                            triclinic=tilt is not None)
    except ValueError as err:
        raise _broken(path, None, str(err)) from None
    return header, index


def _is_number(word: str) -> bool:
    """Return True if `word` reads as a number."""
    try:
        float(word)
    except ValueError:
        return False
    return True


def _read_sections(
    lines: list[str], start: int, path: Path
) -> list[_Section]:
    """Split the body into sections; a line that starts with a number is an
    entry of the section above it."""
    sections: list[_Section] = []
    names: set[str] = set()
    for index in range(start, len(lines)):
        content, comment = _content(lines[index])
        if not content:
            continue
        words = content.split()
        if _is_number(words[0]):
            if not sections:
                raise _broken(path, index + 1,
                              "an entry before any section keyword")
            sections[-1].rows.append((index + 1, words))
            continue
        if content in names:
            raise _broken(path, index + 1, f"section {content!r} repeats")
        names.add(content)
        sections.append(_Section(content, comment, index + 1, []))
    return sections


def _atom_style(atoms: _Section, given: str | None, path: Path) -> str:
    """Return the atom style from the argument and the Atoms comment."""
    named = atoms.comment.split()[0] if atoms.comment else None
    if given is not None:
        if named is not None and given != named:
            raise _broken(path, atoms.line, f"'atom_style' is {given!r} but "
                          f"the file says {named!r}")
        style = given
    elif named is not None:
        style = named
    else:
        raise _broken(path, atoms.line,
                      "the Atoms line names no atom style; pass 'atom_style'")
    if style not in _STYLES:
        raise _broken(path, atoms.line, f"atom style {style!r} is not "
                      f"supported; use one of {tuple(_STYLES)}")
    return style


def _table(
    section: _Section, width: int | tuple[int, ...], count: int, path: Path
) -> FloatArray:
    """Return a section's rows as a float64 array, checking their number and
    width."""
    widths = (width,) if isinstance(width, int) else width
    if len(section.rows) != count:
        raise _broken(path, section.line, f"{section.name} has "
                      f"{len(section.rows)} entries, the header says {count}")
    if not section.rows:
        return np.empty((0, widths[0]))
    found = len(section.rows[0][1])
    for line, words in section.rows:
        if len(words) != found or found not in widths:
            raise _broken(path, line, f"a {section.name} entry needs "
                          f"{' or '.join(map(str, widths))} fields, the same "
                          f"on every line")
    try:
        return np.array([words for _, words in section.rows],
                        dtype=np.float64)
    except ValueError:
        raise _broken(path, section.line, f"a {section.name} entry holds a "
                      f"value that is not a number") from None


def _read_atoms(
    section: _Section, style: str, header: dict[str, Any], path: Path
) -> dict[str, Any]:
    """Read the Atoms section into the fields of :class:`DataFile`."""
    spec = _STYLES[style]
    names = ["ids"] + (["mol_ids"] if spec.molecule else []) + ["types"]
    names += (["charges"] if spec.charge else []) + ["x", "y", "z"]
    rows = _table(section, (len(names), len(names) + 3), header["atoms"],
                  path)
    rows = rows[np.argsort(rows[:, 0], kind="stable")]
    fields: dict[str, Any] = {}
    for column, name in enumerate(names[:-3]):
        values = rows[:, column]
        fields[name] = values if name == "charges" else _integers(
            name, values)
    fields["positions"] = rows[:, len(names) - 3:len(names)]
    if rows.shape[1] > len(names):
        fields["image_flags"] = _integers("image flags",
                                          rows[:, len(names):])
    return fields


def _read_velocities(
    section: _Section, header: dict[str, Any], path: Path
) -> FloatArray:
    """Read the Velocities section as ``atom-ID vx vy vz`` rows."""
    return _table(section, 4, header["atoms"], path)


def _by_atom(
    rows: FloatArray, ids: IntArray, name: str, path: Path
) -> FloatArray:
    """Order ``atom-ID value...`` rows like `ids`, which are sorted."""
    keys = _integers(name, rows[:, 0])
    order = np.argsort(keys, kind="stable")
    if not np.array_equal(keys[order], ids):
        raise _broken(path, None,
                      f"the {name} atom IDs differ from the Atoms atom IDs")
    ordered: FloatArray = rows[order, 1:]
    return ordered


def _read_masses(
    section: _Section, header: dict[str, Any], path: Path
) -> FloatArray:
    """Read the Masses section, one ``type mass`` line per atom type."""
    rows = _table(section, 2, header["atom_types"], path)
    types = _integers("Masses", rows[:, 0])
    order = np.argsort(types, kind="stable")
    expected = np.arange(1, header["atom_types"] + 1)
    if not np.array_equal(types[order], expected):
        raise _broken(path, section.line, f"Masses must give each atom type "
                      f"from 1 to {header['atom_types']} once")
    masses: FloatArray = rows[order, 1]
    return masses


def with_frame(template: DataFile, frame: DumpFrame) -> DataFile:
    """Return `template` with the positions and the box of a dump frame.

    Use it to write the topology that goes with a dump, for analysis: the
    atoms, types, masses, bonds, and angles come from `template`, and the
    positions and box from `frame`.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    template : DataFile
        A data file of the same system, such as the one the run started
        from.
    frame : DumpFrame
        A frame with an ``id`` column and coordinate columns that
        :meth:`~polyphys.simulate.dump.DumpFrame.positions` can use,
        scaled or not.

    Returns
    -------
    DataFile
        `template` with, matched by atom ID: the frame's wrapped positions;
        its image flags if it has ``ix iy iz``, else none; its velocities if
        it has ``vx vy vz``, else none; and its box.

    Raises
    ------
    ValueError
        If the frame has no ``id`` column, its atom IDs differ from those of
        `template`, or its ``type`` column, when present, disagrees with the
        template's types.
    KeyError
        If the frame has no coordinate columns.

    Notes
    -----
    The template's image flags and velocities belong to its own positions,
    so they are never carried over. Time is :math:`O(n_a \\log n_a)`, memory
    :math:`O(n_a)`.

    Examples
    --------
    >>> template = DataFile(
    ...     atom_style="bond", box=Box(lo=np.zeros(3), hi=np.full(3, 10.0)),
    ...     ids=np.array([1, 2]), types=np.array([1, 1]),
    ...     positions=np.array([[1.0, 5.0, 5.0], [2.0, 5.0, 5.0]]),
    ...     n_atom_types=1, mol_ids=np.array([1, 1]),
    ...     bonds=np.array([[1, 1, 2]]), n_bond_types=1)
    >>> frame = DumpFrame(timestep=1000,
    ...                   box=Box(lo=np.zeros(3), hi=np.full(3, 8.0)),
    ...                   columns=("id", "type", "xs", "ys", "zs"),
    ...                   data=np.array([[2, 1, 0.5, 0.5, 0.5],
    ...                                  [1, 1, 0.375, 0.5, 0.5]]))
    >>> moved = with_frame(template, frame)
    >>> moved.positions
    array([[3., 4., 4.],
           [4., 4., 4.]])
    >>> moved.box.hi, moved.bonds
    (array([8., 8., 8.]), array([[1, 1, 2]]))
    """
    if "id" not in frame.columns:
        raise ValueError("the frame has no 'id' column")
    frame_ids = frame.column("id")
    frame_order = np.argsort(frame_ids, kind="stable")
    template_order = np.argsort(template.ids, kind="stable")
    if not np.array_equal(frame_ids[frame_order],
                          template.ids[template_order]):
        raise ValueError("the frame's atom IDs differ from the template's")
    # Row k of the template takes the frame row with the same atom ID.
    rows = np.empty_like(frame_order)
    rows[template_order] = frame_order
    if "type" in frame.columns and not np.array_equal(
        frame.column("type")[rows], template.types
    ):
        raise ValueError("the frame's atom types differ from the template's")
    columns = set(frame.columns)
    image: IntArray | None = None
    if {"ix", "iy", "iz"} <= columns:
        image = np.column_stack(
            [frame.column(n) for n in ("ix", "iy", "iz")]
        ).astype(np.int64)[rows]
    velocities: FloatArray | None = None
    if {"vx", "vy", "vz"} <= columns:
        velocities = np.column_stack(
            [frame.column(n) for n in ("vx", "vy", "vz")]
        ).astype(np.float64)[rows]
    return dataclasses.replace(
        template,
        positions=frame.positions()[rows], box=frame.box,
        image_flags=image, velocities=velocities,
    )
