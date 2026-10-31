"""\
=====================================================================
LAMMPS dump files --- :mod:`polyphys.simulate.dump`
=====================================================================

The :mod:`~polyphys.simulate.dump` module reads the text dump files
written by the LAMMPS ``dump atom`` and ``dump custom`` commands into numpy
arrays, one frame at a time. It needs numpy and the standard library only, so
an install without MDAnalysis can read dumps.

.. versionadded:: 0.6.0

Contents
========

:func:`iter_dump` yields the frames of a file one at a time, and
:func:`read_dump_frame` returns one of them. Each frame is a
:class:`DumpFrame`: its timestep, its box as a
:class:`~polyphys.simulate.core.Box`, the names of the per-atom columns, and
the per-atom values. :meth:`DumpFrame.positions` returns Cartesian coordinates
from whichever coordinate columns the frame has. :class:`BrokenDumpError`, a
:class:`~polyphys.simulate.core.LammpsFormatError`, reports a file that is cut
off or does not follow the format.

Format
======

A dump is a sequence of frames, each a header and one line per atom::

    ITEM: TIMESTEP
    1000
    ITEM: NUMBER OF ATOMS
    2
    ITEM: BOX BOUNDS pp pp pp
    0.0 10.0
    0.0 10.0
    0.0 10.0
    ITEM: ATOMS id type x y z
    1 1 1.0 5.0 5.0
    2 1 2.0 5.0 5.0

- ``ITEM: UNITS`` and ``ITEM: TIME``, each followed by its value, may come
  before ``ITEM: TIMESTEP``; ``dump_modify units yes`` writes the first in
  the first frame only, and ``dump_modify time yes`` the second in every
  frame.
- ``ITEM: BOX BOUNDS`` ends with three boundary flags, such as ``pp ff pp``,
  and is followed by one ``lo hi`` line per dimension. For a restricted
  triclinic box it reads ``ITEM: BOX BOUNDS xy xz yz pp pp pp`` and its lines
  are ``xlo_bound xhi_bound xy``, ``ylo_bound yhi_bound xz``, and
  ``zlo_bound zhi_bound yz``: the bounds of the orthogonal box that encloses
  the triclinic box, then a tilt factor. The general triclinic form,
  ``ITEM: BOX BOUNDS abc origin``, is not supported.
- ``ITEM: ATOMS`` names the per-atom columns. Unless ``dump_modify sort`` is
  set, the atom lines are in no fixed order.

Coordinates
===========

``x y z`` are Cartesian coordinates wrapped into the box, and ``xu yu zu``
are unwrapped by the image flags ``ix iy iz``. ``xs ys zs`` and ``xsu ysu
zsu`` are the same coordinates scaled: fractions of the box edge vectors,
converted as :mod:`polyphys.simulate.core` describes. ``dump atom`` writes
``id type xs ys zs`` unless ``dump_modify scale no`` is set. The module
converts no units; values are in the ``units`` of the run.

Complexity
==========

Reading a file takes :math:`O(n_f n_a n_c)` time for :math:`n_f` frames of
:math:`n_a` atoms with :math:`n_c` columns. :func:`iter_dump` holds one frame
at a time, :math:`O(n_a n_c)` memory, plus the frames the caller keeps.

References
==========

The format follows the LAMMPS documentation of `dump
<https://docs.lammps.org/dump.html>`_ and `Howto_triclinic
<https://docs.lammps.org/Howto_triclinic.html>`_, and these lines of the
LAMMPS source, release ``stable_30Sep2026_update1``; the coordinate formulas
are cited in :mod:`polyphys.simulate.core`:

- Frame headers: `dump_custom.cpp, lines 660 to 725
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/dump_custom.cpp#L660-L725>`_,
  and the same in `dump_atom.cpp, lines 388 to 450
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/dump_atom.cpp#L388-L450>`_;
  the bounds of a triclinic box: `dump.cpp, lines 327 to 344
  <https://github.com/lammps/lammps/blob/098ffdefa525cabc82272925de02e9884f58c657/src/dump.cpp#L327-L344>`_.

Dependencies
============
- Python standard library: `contextlib`, `dataclasses`, `itertools`,
  `logging`, `pathlib`
- Third-party: `numpy`
- PolyPhys: :mod:`polyphys.simulate.core`
"""

import logging
from contextlib import closing
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import IO, Generator

import numpy as np
import numpy.typing as npt

from polyphys.simulate.core import (
    Box,
    FloatArray,
    LammpsFormatError,
    _frozen,
    _open_text,
)

__all__ = [
    "BrokenDumpError",
    "DumpFrame",
    "iter_dump",
    "read_dump_frame",
]

logger = logging.getLogger(__name__)

# Per-atom attributes of dump custom that hold integers (LAMMPS dump command,
# "Per-atom attributes used as arguments to the custom and cfg styles").
_INTEGER_COLUMNS = frozenset(
    {"id", "mol", "proc", "procp1", "type", "ix", "iy", "iz"}
)
_WRAPPED = (("x", "y", "z"), ("xs", "ys", "zs"))
_UNWRAPPED = (("xu", "yu", "zu"), ("xsu", "ysu", "zsu"))
_SCALED = frozenset({("xs", "ys", "zs"), ("xsu", "ysu", "zsu")})
_IMAGE = ("ix", "iy", "iz")
_ITEM = "ITEM:"


class BrokenDumpError(LammpsFormatError):
    """A LAMMPS dump is cut off or does not follow the dump format.

    .. versionadded:: 0.6.0

    Examples
    --------
    A run killed while writing its second frame leaves that frame short of
    atom lines:

    >>> import tempfile
    >>> frame = '''ITEM: TIMESTEP
    ... {step}
    ... ITEM: NUMBER OF ATOMS
    ... 2
    ... ITEM: BOX BOUNDS pp pp pp
    ... 0 10
    ... 0 10
    ... 0 10
    ... ITEM: ATOMS id type x y z
    ... 1 1 1.0 5.0 5.0
    ... '''
    >>> second_atom = "2 1 2.0 5.0 5.0\\n"
    >>> text = frame.format(step=0) + second_atom + frame.format(step=100)
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "killed.lammpstrj"
    ...     _ = path.write_text(text)
    ...     steps = []
    ...     try:
    ...         for dump_frame in iter_dump(path):
    ...             steps.append(dump_frame.timestep)
    ...     except BrokenDumpError as err:
    ...         print(str(err).split(": ", 1)[1])
    frame 1 (timestep 100), line 21: 1 of 2 atom lines; the dump may be cut off
    >>> steps
    [0]
    """


@dataclass(frozen=True, eq=False)
class DumpFrame:
    """One frame of a LAMMPS dump.

    .. versionadded:: 0.6.0

    Attributes
    ----------
    timestep : int
        The ``ITEM: TIMESTEP`` value.
    box : Box
        The simulation box, with the boundary flags of the frame; for a
        triclinic box, its own bounds, not those the dump writes.
    columns : tuple of str
        The per-atom column names of ``ITEM: ATOMS``, in file order, such as
        ``('id', 'type', 'xs', 'ys', 'zs')``.
    data : numpy.ndarray
        Shape ``(n_atoms, n_columns)``, float64: one row per atom, one column
        per name in `columns`. Integer columns such as ``id`` are exact up to
        :math:`2^{53}`. Read-only.
    time : float or None
        The ``ITEM: TIME`` value (simulated time), or None if the frame has
        none.
    units : str or None
        The ``ITEM: UNITS`` value, such as ``'lj'``, or None. LAMMPS writes it
        in the first frame only; :func:`iter_dump` gives it to every frame.

    Raises
    ------
    ValueError
        If `data` is not two-dimensional with one column per name in
        `columns`, or a column name repeats.

    Examples
    --------
    >>> box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
    >>> frame = DumpFrame(timestep=0, box=box,
    ...                   columns=("id", "type", "xs", "ys", "zs"),
    ...                   data=np.array([[1, 1, 0.1, 0.5, 0.5],
    ...                                  [2, 1, 0.2, 0.5, 0.5]]))
    >>> frame.n_atoms
    2
    >>> frame.column("id")
    array([1, 2])
    >>> frame.positions()
    array([[1., 5., 5.],
           [2., 5., 5.]])
    """

    timestep: int
    box: Box
    columns: tuple[str, ...]
    data: FloatArray
    time: float | None = None
    units: str | None = None

    def __post_init__(self) -> None:
        data = np.array(self.data, dtype=np.float64)
        if data.ndim != 2 or data.shape[1] != len(self.columns):
            raise ValueError(
                f"'data' must have shape (n_atoms, {len(self.columns)}), "
                f"got {data.shape}"
            )
        if len(set(self.columns)) != len(self.columns):
            raise ValueError(f"column names repeat: {self.columns}")
        object.__setattr__(self, "columns", tuple(self.columns))
        object.__setattr__(self, "data", _frozen(data, np.float64))

    @property
    def n_atoms(self) -> int:
        """The number of atoms in the frame."""
        return int(self.data.shape[0])

    def column(self, name: str) -> npt.NDArray[np.float64 | np.int64]:
        """Return the values of one per-atom column.

        Parameters
        ----------
        name : str
            A name in :attr:`columns`.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_atoms,)``: int64 for the integer attributes ``id``,
            ``mol``, ``proc``, ``procp1``, ``type``, ``ix``, ``iy``, and
            ``iz``, else a read-only float64 view of :attr:`data`.

        Raises
        ------
        KeyError
            If the frame has no column `name`.
        """
        try:
            index = self.columns.index(name)
        except ValueError:
            raise KeyError(
                f"no column {name!r}; the frame has {self.columns}"
            ) from None
        values = self.data[:, index]
        if name in _INTEGER_COLUMNS:
            return values.astype(np.int64)
        return values

    def positions(self, unwrapped: bool = False) -> FloatArray:
        r"""Return the Cartesian coordinates of the atoms.

        Parameters
        ----------
        unwrapped : bool
            False for coordinates wrapped into the box; True for coordinates
            unwrapped by the image flags.

        Returns
        -------
        numpy.ndarray
            Shape ``(n_atoms, 3)``, float64, in the length unit of the run.

        Raises
        ------
        KeyError
            If the frame lacks the columns for the requested coordinates.

        Notes
        -----
        Wrapped coordinates come from ``x y z``, else from ``xs ys zs``.
        Unwrapped ones come from ``xu yu zu``, else from ``xsu ysu zsu``,
        else from the wrapped coordinates and the image flags ``ix iy iz``.
        Scaled coordinates are converted by
        :meth:`~polyphys.simulate.core.Box.to_cartesian` and image flags
        applied by :meth:`~polyphys.simulate.core.Box.unwrap`, as LAMMPS's
        ``Domain::lamda2x`` and ``Domain::unmap`` compute them. Time and
        memory are :math:`O(n_a)` for :math:`n_a` atoms.

        Examples
        --------
        >>> box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
        >>> frame = DumpFrame(timestep=0, box=box,
        ...                   columns=("id", "x", "y", "z", "ix", "iy", "iz"),
        ...                   data=np.array([[1, 9.5, 5.0, 5.0, -1, 0, 0]]))
        >>> frame.positions()
        array([[9.5, 5. , 5. ]])
        >>> frame.positions(unwrapped=True)
        array([[-0.5,  5. ,  5. ]])
        """
        options = _UNWRAPPED if unwrapped else _WRAPPED
        for names in options:
            if set(names) <= set(self.columns):
                return self._coordinates(names)
        if unwrapped and set(_IMAGE) <= set(self.columns):
            for names in _WRAPPED:
                if set(names) <= set(self.columns):
                    image = np.column_stack(
                        [self.column(name) for name in _IMAGE])
                    return self.box.unwrap(self._coordinates(names), image)
        wanted = [" ".join(names) for names in options]
        if unwrapped:
            wanted.append("x y z or xs ys zs, with ix iy iz")
        raise KeyError(
            f"no coordinate columns for unwrapped={unwrapped}: need one of "
            f"{wanted}; the frame has {self.columns}"
        )

    def _coordinates(self, names: tuple[str, str, str]) -> FloatArray:
        """Return three coordinate columns as Cartesian coordinates."""
        values = np.column_stack([self.data[:, self.columns.index(name)]
                                  for name in names])
        if names in _SCALED:
            return self.box.to_cartesian(values)
        return values


def iter_dump(
    path: str | Path, *, sort_by_id: bool = True
) -> Generator[DumpFrame, None, None]:
    """Yield the frames of a LAMMPS dump file one at a time.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    path : str or pathlib.Path
        A text dump of ``dump atom`` or ``dump custom``; a name ending in
        ``.gz`` is read through gzip.
    sort_by_id : bool
        True to order each frame's rows by the ``id`` column; LAMMPS writes
        atoms in no fixed order unless ``dump_modify sort`` is set. A frame
        without an ``id`` column keeps its file order.

    Yields
    ------
    DumpFrame
        The frames in file order, each read only when it is requested.

    Raises
    ------
    BrokenDumpError
        If a frame is cut off or does not follow the format; the frames
        before it have been yielded. The message names the file, the frame
        (from 0), its timestep when known, and the line (from 1).
    FileNotFoundError
        If `path` does not exist, raised when the first frame is requested.

    Notes
    -----
    A frame whose timestep is not larger than the one before, as when a
    restarted run writes its restart step again, is yielded and logged as a
    warning.

    Time is :math:`O(n_a n_c)` per frame for :math:`n_a` atoms and
    :math:`n_c` columns; memory is that of one frame, plus the frames the
    caller keeps. Each frame's atom lines are parsed by one split and one
    conversion to float64, with no loop over atoms in Python.

    Examples
    --------
    >>> import tempfile
    >>> text = '''ITEM: TIMESTEP
    ... 0
    ... ITEM: NUMBER OF ATOMS
    ... 2
    ... ITEM: BOX BOUNDS pp pp pp
    ... 0.0 10.0
    ... 0.0 10.0
    ... 0.0 10.0
    ... ITEM: ATOMS id type xs ys zs
    ... 2 1 0.2 0.5 0.5
    ... 1 1 0.1 0.5 0.5
    ... '''
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "chain.lammpstrj"
    ...     _ = path.write_text(text)
    ...     frames = list(iter_dump(path))
    >>> [(frame.timestep, frame.n_atoms) for frame in frames]
    [(0, 2)]
    >>> frames[0].column("id")
    array([1, 2])
    >>> frames[0].positions()
    array([[1., 5., 5.],
           [2., 5., 5.]])
    """
    path = Path(path)
    with _open_text(path) as stream:
        lines = _Lines(stream, path)
        units: str | None = None
        previous: int | None = None
        index = 0
        while True:
            frame = _read_frame(lines, index, units)
            if frame is None:
                return
            units = frame.units
            if previous is not None and frame.timestep <= previous:
                logger.warning(
                    "%s: frame %d has timestep %d, not after the previous "
                    "frame's %d", path, index, frame.timestep, previous,
                )
            previous = frame.timestep
            if sort_by_id and "id" in frame.columns:
                order = np.argsort(frame.column("id"), kind="stable")
                frame = DumpFrame(
                    timestep=frame.timestep, box=frame.box,
                    columns=frame.columns, data=frame.data[order],
                    time=frame.time, units=frame.units,
                )
            yield frame
            index += 1


def read_dump_frame(
    path: str | Path, index: int = 0, *, sort_by_id: bool = True
) -> DumpFrame:
    """Read one frame of a LAMMPS dump file.

    .. versionadded:: 0.6.0

    Parameters
    ----------
    path : str or pathlib.Path
        As for :func:`iter_dump`.
    index : int
        Position of the frame in the file, from 0.
    sort_by_id : bool
        As for :func:`iter_dump`.

    Returns
    -------
    DumpFrame

    Raises
    ------
    IndexError
        If `index` is negative or the file has no frame `index`.
    BrokenDumpError
        If a frame up to and including frame `index` is broken.

    Notes
    -----
    Every frame before `index` is parsed: time is
    :math:`O((index + 1) n_a n_c)`, memory that of one frame.

    Examples
    --------
    >>> import tempfile
    >>> frame = '''ITEM: TIMESTEP
    ... {step}
    ... ITEM: NUMBER OF ATOMS
    ... 1
    ... ITEM: BOX BOUNDS pp pp pp
    ... 0 10
    ... 0 10
    ... 0 10
    ... ITEM: ATOMS id type x y z
    ... 1 1 {x} 5.0 5.0
    ... '''
    >>> text = frame.format(step=0, x=1.0) + frame.format(step=500, x=1.5)
    >>> with tempfile.TemporaryDirectory() as tmp:
    ...     path = Path(tmp) / "bead.lammpstrj"
    ...     _ = path.write_text(text)
    ...     second = read_dump_frame(path, 1)
    >>> second.timestep, second.positions()
    (500, array([[1.5, 5. , 5. ]]))
    """
    if index < 0:
        raise IndexError(f"'index' must be 0 or more, got {index}")
    count = 0
    with closing(iter_dump(path, sort_by_id=sort_by_id)) as frames:
        for count, frame in enumerate(frames, start=1):
            if count - 1 == index:
                return frame
    raise IndexError(f"{path} has {count} frames; there is no frame {index}")


class _Lines:
    """The lines of an open dump, with the number of the last line read."""

    def __init__(self, stream: IO[str], path: Path) -> None:
        self._lines = iter(stream)
        self.path = path
        self.number = 0

    def read(self) -> str | None:
        """Return the next line, or None at the end of the file."""
        line = next(self._lines, None)
        if line is not None:
            self.number += 1
        return line

    def take(self, count: int) -> list[str]:
        """Return up to `count` lines; fewer only at the end of the file."""
        taken = list(islice(self._lines, count))
        self.number += len(taken)
        return taken


def _broken(
    lines: _Lines, index: int, timestep: int | None, message: str
) -> BrokenDumpError:
    """Return an error naming the file, frame, timestep, and line."""
    step = "" if timestep is None else f" (timestep {timestep})"
    return BrokenDumpError(
        f"{lines.path}: frame {index}{step}, line {lines.number}: {message}"
    )


def _read_frame(
    lines: _Lines, index: int, units: str | None
) -> DumpFrame | None:
    """Read one frame, or return None at the end of the file."""
    line = lines.read()
    while line is not None and not line.strip():
        line = lines.read()
    if line is None:
        return None
    time: float | None = None
    while True:
        item = _item(lines, index, line)
        if item == "TIMESTEP":
            break
        if item == "UNITS":
            units = _first_field(lines, index, None, "UNITS")
        elif item == "TIME":
            time = _read_float(lines, index, None, "TIME")
        else:
            raise _broken(lines, index, None,
                          f"expected ITEM: TIMESTEP, got {line.strip()!r}")
        line = _next_line(lines, index, None)
    timestep = _read_int(lines, index, None, "TIMESTEP")
    _expect(lines, index, timestep, "NUMBER OF ATOMS")
    n_atoms = _read_int(lines, index, timestep, "NUMBER OF ATOMS")
    if n_atoms < 0:
        raise _broken(lines, index, timestep,
                      f"negative number of atoms {n_atoms}")
    header = _expect(lines, index, timestep, "BOX BOUNDS")
    box = _read_box(lines, index, timestep, header)
    columns = tuple(_expect(lines, index, timestep, "ATOMS").split()[2:])
    if not columns or len(set(columns)) != len(columns):
        raise _broken(lines, index, timestep,
                      f"ITEM: ATOMS names no columns or repeats one: "
                      f"{columns}")
    data = _read_atoms(lines, index, timestep, n_atoms, len(columns))
    return DumpFrame(timestep=timestep, box=box, columns=columns, data=data,
                     time=time, units=units)


def _next_line(lines: _Lines, index: int, timestep: int | None) -> str:
    """Return the next line; a frame that ends here is cut off."""
    line = lines.read()
    if line is None:
        raise _broken(lines, index, timestep,
                      "the file ends inside the frame header; "
                      "the dump may be cut off")
    return line


def _item(lines: _Lines, index: int, line: str) -> str:
    """Return the name of an ``ITEM:`` line, such as ``'BOX BOUNDS'``."""
    if not line.startswith(_ITEM):
        raise _broken(lines, index, None,
                      f"expected an ITEM: line, got {line.strip()!r}")
    words = line[len(_ITEM):].split()
    for name in ("NUMBER OF ATOMS", "BOX BOUNDS"):
        if words[:len(name.split())] == name.split():
            return name
    return words[0] if words else ""


def _expect(
    lines: _Lines, index: int, timestep: int | None, name: str
) -> str:
    """Read the ``ITEM: <name>`` line and return it."""
    line = _next_line(lines, index, timestep)
    if not line.startswith(_ITEM) or _item(lines, index, line) != name:
        raise _broken(lines, index, timestep,
                      f"expected ITEM: {name}, got {line.strip()!r}")
    return line


def _first_field(
    lines: _Lines, index: int, timestep: int | None, name: str
) -> str:
    """Return the first field of the value line of ``ITEM: <name>``."""
    words = _next_line(lines, index, timestep).split()
    if not words:
        raise _broken(lines, index, timestep,
                      f"ITEM: {name} is not followed by a value")
    return words[0]


def _read_int(
    lines: _Lines, index: int, timestep: int | None, name: str
) -> int:
    """Read the integer value of ``ITEM: <name>``."""
    field = _first_field(lines, index, timestep, name)
    try:
        return int(field)
    except ValueError:
        raise _broken(lines, index, timestep,
                      f"ITEM: {name} is not followed by an integer, got "
                      f"{field!r}") from None


def _read_float(
    lines: _Lines, index: int, timestep: int | None, name: str
) -> float:
    """Read the floating-point value of ``ITEM: <name>``."""
    field = _first_field(lines, index, timestep, name)
    try:
        return float(field)
    except ValueError:
        raise _broken(lines, index, timestep,
                      f"ITEM: {name} is not followed by a number, got "
                      f"{field!r}") from None


def _read_box(
    lines: _Lines, index: int, timestep: int, header: str
) -> Box:
    """Read the three bound lines of ``ITEM: BOX BOUNDS``.

    For a restricted triclinic box, the dump stores the bounds of the
    orthogonal box that encloses it; the box's own bounds are
    :math:`x_{lo} = x_{lo}^{bound} - \\min(0, xy, xz, xy + xz)`,
    :math:`x_{hi} = x_{hi}^{bound} - \\max(0, xy, xz, xy + xz)`,
    :math:`y_{lo} = y_{lo}^{bound} - \\min(0, yz)`,
    :math:`y_{hi} = y_{hi}^{bound} - \\max(0, yz)`, and the z bounds
    (LAMMPS ``Howto_triclinic``).
    """
    words = header.split()[3:]
    if words[:2] == ["abc", "origin"]:
        raise _broken(lines, index, timestep,
                      "general triclinic boxes (ITEM: BOX BOUNDS abc origin, "
                      "dump_modify triclinic/general) are not supported")
    triclinic = words[:3] == ["xy", "xz", "yz"]
    flags = words[3:] if triclinic else words
    if len(flags) != 3:
        raise _broken(lines, index, timestep,
                      f"ITEM: BOX BOUNDS needs three boundary flags, got "
                      f"{header.strip()!r}")
    width = 3 if triclinic else 2
    rows = []
    for _ in range(3):
        fields = _next_line(lines, index, timestep).split()
        try:
            row = [float(field) for field in fields]
        except ValueError:
            row = []
        if len(row) != width:
            raise _broken(lines, index, timestep,
                          f"a box bound line needs {width} numbers, got "
                          f"{' '.join(fields)!r}")
        rows.append(row)
    bounds = np.array(rows, dtype=np.float64)
    tilt = bounds[:, 2].copy() if triclinic else np.zeros(3)
    lo = bounds[:, 0].copy()
    hi = bounds[:, 1].copy()
    if triclinic:
        xy, xz, yz = tilt
        lo -= [min(0.0, xy, xz, xy + xz), min(0.0, yz), 0.0]
        hi -= [max(0.0, xy, xz, xy + xz), max(0.0, yz), 0.0]
    try:
        return Box(lo=lo, hi=hi, tilt=tilt,
                   boundary=(flags[0], flags[1], flags[2]),
                   triclinic=triclinic)
    except ValueError as err:
        raise _broken(lines, index, timestep, str(err)) from None


def _read_atoms(
    lines: _Lines, index: int, timestep: int, n_atoms: int, n_columns: int
) -> FloatArray:
    """Read a frame's atom lines into an ``(n_atoms, n_columns)`` array."""
    first = lines.number + 1
    rows = lines.take(n_atoms)
    if len(rows) < n_atoms:
        raise _broken(lines, index, timestep,
                      f"{len(rows)} of {n_atoms} atom lines; "
                      "the dump may be cut off")
    fields = "".join(rows).split()
    try:
        values = np.asarray(fields, dtype=np.float64)
    except ValueError:
        values = np.empty(0)
    if values.size != n_atoms * n_columns:
        for offset, row in enumerate(rows):
            if not _is_atom_line(row, n_columns):
                lines.number = first + offset
                break
        raise _broken(lines, index, timestep,
                      f"an atom line needs {n_columns} numbers, got "
                      f"{rows[lines.number - first].strip()!r}")
    return values.reshape(n_atoms, n_columns)


def _is_atom_line(row: str, n_columns: int) -> bool:
    """Return True if `row` holds `n_columns` numbers."""
    fields = row.split()
    try:
        [float(field) for field in fields]
    except ValueError:
        return False
    return len(fields) == n_columns
