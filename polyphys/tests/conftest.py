"""Shared fixtures for the PolyPhys test suite.

Every fixture here writes tier A data: small synthetic files generated in the
test's own ``tmp_path``, carrying no simulation data.
``docs/source/testing.md`` describes the three test-data tiers and when each
is used.

Fixtures
--------
make_temp_file
    An empty file with a given name.
make_lineage_tree
    A synthetic lineage tree of ``segment`` or ``whole`` files named by a
    chosen parser class.
make_lammps_data, make_lammps_dump
    A tiny LAMMPS data file in the ``atomic``, ``bond``, or ``angle`` atom
    style, and a tiny LAMMPS dump file in the ``atom`` style with unscaled
    coordinates.
private_data_dir
    The directory named by ``POLYPHYS_DATA_DIR``, for tests marked
    ``private_data``.

Markers
-------
private_data
    The test reads the maintainer's local full-size data (tier C). It is
    skipped when ``POLYPHYS_DATA_DIR`` is unset, which is always the case in
    CI, and it errors when the variable names something that is not a
    directory.
"""

import os
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import pytest

DATA_DIR_ENV = "POLYPHYS_DATA_DIR"


def pytest_configure(config: pytest.Config) -> None:
    """Register the ``private_data`` marker.

    The marker is registered here rather than in ``pyproject.toml`` so that a
    nested pytest run (``pytester``) that loads this file also knows it.
    """
    config.addinivalue_line(
        "markers",
        "private_data: reads the maintainer's local data under "
        f"{DATA_DIR_ENV}; skipped when the variable is unset",
    )


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Skip or fail a ``private_data`` test according to ``POLYPHYS_DATA_DIR``.

    Unset or empty: the test is skipped. Set to a path that is not a
    directory: the test errors in setup, so a misconfigured machine is not
    mistaken for one without data.
    """
    if item.get_closest_marker("private_data") is None:
        return
    data_dir = os.environ.get(DATA_DIR_ENV, "")
    if data_dir == "":
        pytest.skip(f"{DATA_DIR_ENV} is not set")
    if not Path(data_dir).is_dir():
        pytest.fail(
            f"{DATA_DIR_ENV}={data_dir!r} is not a directory", pytrace=False
        )


@pytest.fixture
def private_data_dir() -> Path:
    """Return the directory named by ``POLYPHYS_DATA_DIR``.

    Use it only in a test marked ``@pytest.mark.private_data``; the marker
    skips the test before this fixture runs when the variable is unset.

    Examples
    --------
    .. code-block:: python

        @pytest.mark.private_data
        def test_space_matches_reference(private_data_dir):
            space_dir = private_data_dir / "SumRuleCyl-allInOne"
    """
    return Path(os.environ[DATA_DIR_ENV])


@pytest.fixture
def make_temp_file(tmp_path):
    """
    Create an empty temporary file with a given filename and return its full
    path. Useful for testing parsers that expect a file path as input.
    """

    def _make_temp_file(filename: str):
        f = tmp_path / filename
        f.touch()
        return str(f)

    return _make_temp_file


# Synthetic lineage tree

# A value for every filename keyword of the parser classes in
# polyphys.manage.parser. They are small, so every derived volume fraction
# stays well below 1, and they carry no meaning beyond being parseable.
_KEYWORD_DEFAULTS: dict[str, float] = {
    "N": 20,
    "nm": 2,
    "ns": 20,
    "nl": 2,
    "nh": 2,
    "nc": 100,
    "am": 1.0,
    "ac": 1.0,
    "al": 2.0,
    "ml": 8.0,
    "sig": 1.0,
    "hl": 10.0,
    "l": 10.0,
    "r": 5.5,
    "lz": 25.0,
    "sd": 1.0,
    "epsilon": 5.0,
    "epss": 1.0,
    "epsl": 1.0,
    "kbmm": 1.0,
    "epshc": 1.0,
    "dt": 0.005,
    "bdump": 1000,
    "adump": 5000,
    "tdump": 1000,
}
# Keywords the tree sets itself, one per lineage level.
_INDEX_KEYWORDS = ("ens", "j")


@dataclass(frozen=True)
class LineageTree:
    """The files written by :func:`make_lineage_tree`.

    Attributes
    ----------
    root : Path
        Directory holding every file. Its name contains hyphens, as a lineage
        directory name does.
    lineage : {'segment', 'whole'}
        The lineage level of every file.
    group : str
        The particle group in every filename.
    prop : str
        The property in every filename.
    paths : tuple of Path
        Every file, ordered by ensemble, then whole, then segment.
    index : mapping of Path to tuple of int
        For each file, its ``(ensemble, whole, segment)`` indices, all
        starting at 1. ``whole`` equals the ``ens`` value in the name and
        ``segment`` the ``j`` value; ``segment`` is 0 for a whole file.
        Read-only, like the rest of the record.
    """

    root: Path
    lineage: str
    group: str
    prop: str
    paths: tuple[Path, ...]
    index: Mapping[Path, tuple[int, int, int]]


def lineage_tree_value(
    ensemble: int, whole: int, segment: int, n_rows: int
) -> np.ndarray:
    """Return the contents of one lineage-tree file.

    Row ``k`` holds ``1000 * ensemble + 100 * whole + 10 * segment + k``, so
    each value names its file and row, and sums and concatenations over a
    tree can be written down by hand.

    Parameters
    ----------
    ensemble, whole, segment : int
        The file's indices, as in :attr:`LineageTree.index`; each from 0 to 9.
    n_rows : int
        Number of rows, from 1 to 10.

    Returns
    -------
    numpy.ndarray
        Float64 array of shape ``(n_rows,)``.

    Examples
    --------
    >>> lineage_tree_value(2, 3, 1, 3)
    array([2310., 2311., 2312.])
    """
    base = 1000 * ensemble + 100 * whole + 10 * segment
    return base + np.arange(n_rows, dtype=np.float64)


def write_lineage_tree(
    root: Path,
    parser: type,
    *,
    group: Optional[str] = None,
    n_ensembles: int = 2,
    n_wholes: int = 2,
    n_segments: int = 2,
    lineage: str = "segment",
    ext: str = ".npy",
    prop: str = "gyrTMon",
    n_rows: int = 3,
    attrs: Optional[Mapping[str, float]] = None,
    vary: str = "nc",
    vary_step: float = 100,
) -> LineageTree:
    """Write a synthetic lineage tree of tiny files named for `parser`.

    The tree has ``n_ensembles`` state points that differ only in the
    keyword `vary`, ``n_wholes`` replicates (``ens1``, ``ens2``, ...) per
    state point, and, at the ``segment`` level, ``n_segments`` segments
    (``.j01``, ``.j02``, ...) per whole. Names follow the parser's ``whole``
    keyword order, with the values written in the file as given:

    - segment: ``<whole>.j<NN>.<group>.<prop><ext>``, as the probe phase
      writes them;
    - whole: ``<whole>-<group>-<prop><ext>``, as
      :func:`polyphys.manage.organizer.save_artifact` writes them.

    File contents follow :func:`lineage_tree_value`. A ``.npy`` file holds
    the array; a ``.csv`` file holds one column named `prop`.

    Parameters
    ----------
    root : Path
        Parent directory; the tree goes in a new subdirectory of it.
    parser : type
        A parser class from :mod:`polyphys.manage.parser`.
    group : str, optional
        Particle group; defaults to the parser's first group.
    n_ensembles, n_wholes, n_segments : int
        Tree sizes, each from 1 to 9.
    lineage : {'segment', 'whole'}
        Lineage level of the files.
    ext : {'.npy', '.csv'}
        File format.
    prop : str
        Property name in the filenames and the ``.csv`` column.
    n_rows : int
        Rows per file, from 1 to 10.
    attrs : mapping of str to float, optional
        Keyword values that override the defaults, keyed by the short
        keyword (``'r'``, ``'nc'``, ...). ``'ens'`` and ``'j'`` are set by
        the tree and cannot be given.
    vary : str
        The keyword that differs between ensembles.
    vary_step : float
        Ensemble ``e`` (from 1) gets ``vary`` equal to its base value plus
        ``(e - 1) * vary_step``.

    Returns
    -------
    LineageTree

    Raises
    ------
    ValueError
        If a size is out of range, `lineage` or `ext` is unknown, `attrs`
        sets ``'ens'`` or ``'j'``, or `vary` or a key of `attrs` is not a
        keyword of the parser's ``whole`` names.

    Notes
    -----
    Writes ``n_ensembles * n_wholes * n_segments`` files (``n_segments`` is
    1 at the ``whole`` level) of ``n_rows`` values each: time and disk space
    are O(n_ensembles * n_wholes * n_segments * n_rows).

    The parser's ``_genealogy_attributes`` table gives the keyword order.
    Take expected lineage names in a test from the parser applied to a
    written filename, not from the strings built here, so the test checks
    the code under test against the parser's own lineage rules.

    Examples
    --------
    .. code-block:: python

        def test_example(make_lineage_tree):
            tree = make_lineage_tree(SumRuleCyl, n_ensembles=1)
            assert len(tree.paths) == 4
    """
    sizes = {
        "n_ensembles": n_ensembles,
        "n_wholes": n_wholes,
        "n_segments": n_segments,
    }
    for name, size in sizes.items():
        if not 1 <= size <= 9:
            raise ValueError(f"'{name}' must be from 1 to 9, got {size}")
    if not 1 <= n_rows <= 10:
        raise ValueError(f"'n_rows' must be from 1 to 10, got {n_rows}")
    if lineage not in ("segment", "whole"):
        raise ValueError(
            f"'lineage' must be 'segment' or 'whole', " f"got {lineage!r}"
        )
    if ext not in (".npy", ".csv"):
        raise ValueError(f"'ext' must be '.npy' or '.csv', got {ext!r}")
    attrs = dict(attrs or {})
    reserved = sorted(set(attrs) & set(_INDEX_KEYWORDS))
    if reserved:
        raise ValueError(f"'attrs' cannot set {reserved}")

    keywords = [
        kw
        for kw in parser._genealogy_attributes["whole"].values()
        if kw not in _INDEX_KEYWORDS
    ]
    unknown = sorted((set(attrs) | {vary}) - set(keywords))
    if unknown:
        raise ValueError(
            f"{unknown} are not keywords of {parser.__name__} whole names"
        )
    values = {kw: attrs.get(kw, _KEYWORD_DEFAULTS[kw]) for kw in keywords}
    group = group or parser._groups[0]
    n_segs = n_segments if lineage == "segment" else 1
    stage = "segment" if lineage == "segment" else "wholeSim"
    tree_dir = root / f"{parser.__name__}-{group}-{stage}"
    tree_dir.mkdir(parents=True)

    paths: list[Path] = []
    index: dict[Path, tuple[int, int, int]] = {}
    for e in range(1, n_ensembles + 1):
        state = dict(values, **{vary: values[vary] + (e - 1) * vary_step})
        ensemble_long = "".join(f"{kw}{state[kw]}" for kw in keywords)
        for w in range(1, n_wholes + 1):
            whole = f"{ensemble_long}ens{w}"
            for s in range(1, n_segs + 1):
                if lineage == "segment":
                    filename = f"{whole}.j{s:02d}.{group}.{prop}{ext}"
                    seg = s
                else:
                    filename = f"{whole}-{group}-{prop}{ext}"
                    seg = 0
                path = tree_dir / filename
                data = lineage_tree_value(e, w, seg, n_rows)
                if ext == ".npy":
                    np.save(path, data)
                else:
                    pd.DataFrame({prop: data}).to_csv(path, index=False)
                paths.append(path)
                index[path] = (e, w, seg)
    return LineageTree(tree_dir, lineage, group, prop, tuple(paths),
                       MappingProxyType(index))


@pytest.fixture
def make_lineage_tree(tmp_path: Path) -> Callable[..., LineageTree]:
    """Factory: write a lineage tree under ``tmp_path``.

    Takes the arguments of :func:`write_lineage_tree` after `root`. Each call
    writes a new tree, so give trees in one test different `group`, `lineage`,
    or parser classes.
    """

    def _make(parser: type, **kwargs: Any) -> LineageTree:
        return write_lineage_tree(tmp_path, parser, **kwargs)

    return _make


# Tiny LAMMPS files
#
# Formats from the LAMMPS manual, read_data command (Atoms, Bonds, Angles,
# and Masses sections; header keywords) and dump command (dump style atom,
# BOX BOUNDS boundary flags): https://docs.lammps.org/read_data.html and
# https://docs.lammps.org/dump.html. The order of the dump's ITEM lines
# matches what MDAnalysis 2.9.0's LAMMPSDUMP reader parses
# (MDAnalysis.coordinates.LAMMPS.DumpReader).

_ATOM_STYLES = ("atomic", "bond", "angle")
_BOX = ((0.0, 10.0), (0.0, 10.0), (0.0, 10.0))


def _default_chain(n_atoms: int) -> np.ndarray:
    """Return a straight chain along x, one unit between neighbours."""
    positions = np.full((n_atoms, 3), 5.0)
    positions[:, 0] = 1.0 + np.arange(n_atoms)
    return positions


def _box_lines(box: Sequence[Sequence[float]]) -> list[str]:
    """Return the three ``lo hi`` lines of an orthogonal box."""
    return [f"{lo:.6f} {hi:.6f}" for lo, hi in box]


def write_lammps_data(
    path: Path,
    atom_style: str = "bond",
    positions: Optional[np.ndarray] = None,
    atom_types: Optional[Sequence[int]] = None,
    mol_ids: Optional[Sequence[int]] = None,
    box: Sequence[Sequence[float]] = _BOX,
) -> Path:
    """Write a tiny LAMMPS data file of one linear chain.

    The ``Atoms`` lines are ``atom-ID atom-type x y z`` for ``atomic`` and
    ``atom-ID molecule-ID atom-type x y z`` for ``bond`` and ``angle``. A
    ``bond`` or ``angle`` file bonds atom ``i`` to ``i + 1``; an ``angle``
    file also has the angle ``i, i + 1, i + 2``. Bond and angle types are 1.
    Every atom type gets mass 1.0. Coordinates are written with six decimals.

    Parameters
    ----------
    path : Path
        File to write.
    atom_style : {'atomic', 'bond', 'angle'}
        LAMMPS atom style; it is also written after ``Atoms #``.
    positions : numpy.ndarray, optional
        Shape ``(n_atoms, 3)``; defaults to 4 atoms one unit apart along x.
    atom_types : sequence of int, optional
        Per-atom types from 1; defaults to all 1.
    mol_ids : sequence of int, optional
        Per-atom molecule IDs, ignored for ``atomic``; defaults to all 1.
    box : sequence of (float, float)
        ``(lo, hi)`` for x, y, and z.

    Returns
    -------
    Path
        `path`.

    Raises
    ------
    ValueError
        If `atom_style` is unknown or `positions` is not ``(n_atoms, 3)``.

    Notes
    -----
    Time and space are O(n_atoms).

    Examples
    --------
    .. code-block:: python

        def test_example(make_lammps_data):
            path = make_lammps_data(atom_style="angle")
    """
    if atom_style not in _ATOM_STYLES:
        raise ValueError(
            f"'atom_style' must be one of {_ATOM_STYLES}, "
            f"got {atom_style!r}"
        )
    pos = _default_chain(4) if positions is None else np.asarray(positions)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(
            f"'positions' must have shape (n_atoms, 3), " f"got {pos.shape}"
        )
    n_atoms = pos.shape[0]
    types = list(atom_types) if atom_types is not None else [1] * n_atoms
    mols = list(mol_ids) if mol_ids is not None else [1] * n_atoms
    bonds = (
        [(i, i + 1) for i in range(1, n_atoms)]
        if atom_style in ("bond", "angle")
        else []
    )
    angles = (
        [(i, i + 1, i + 2) for i in range(1, n_atoms - 1)]
        if atom_style == "angle"
        else []
    )

    header = [f"{n_atoms} atoms", f"{max(types)} atom types"]
    if atom_style in ("bond", "angle"):
        header += [f"{len(bonds)} bonds", "1 bond types"]
    if atom_style == "angle":
        header += [f"{len(angles)} angles", "1 angle types"]
    box_lines = [
        f"{line} {ax}lo {ax}hi" for line, ax in zip(_box_lines(box), "xyz")
    ]
    masses = [f"{t} 1.0" for t in range(1, max(types) + 1)]
    atoms = []
    for i, (xyz, atype, mol) in enumerate(zip(pos, types, mols), start=1):
        coords = " ".join(f"{c:.6f}" for c in xyz)
        mol_col = "" if atom_style == "atomic" else f"{mol} "
        atoms.append(f"{i} {mol_col}{atype} {coords}")

    lines = [
        f"LAMMPS data file, PolyPhys synthetic fixture, "
        f"atom_style {atom_style}",
        "",
    ]
    lines += header + [""] + box_lines
    lines += ["", "Masses", ""] + masses
    lines += ["", f"Atoms # {atom_style}", ""] + atoms
    if bonds:
        lines += ["", "Bonds", ""]
        lines += [f"{k} 1 {i} {j}" for k, (i, j) in enumerate(bonds, 1)]
    if angles:
        lines += ["", "Angles", ""]
        lines += [
            f"{k} 1 {i} {j} {m}" for k, (i, j, m) in enumerate(angles, 1)
        ]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_lammps_dump(
    path: Path,
    positions: Optional[np.ndarray] = None,
    atom_types: Optional[Sequence[int]] = None,
    timesteps: Optional[Sequence[int]] = None,
    box: Sequence[Sequence[float]] = _BOX,
    boundary: str = "pp pp pp",
) -> Path:
    """Write a tiny LAMMPS dump file in the ``atom`` style.

    Each frame is ``ITEM: TIMESTEP``, ``ITEM: NUMBER OF ATOMS``,
    ``ITEM: BOX BOUNDS <boundary>`` with three ``lo hi`` lines, and
    ``ITEM: ATOMS id type x y z`` with one line per atom: the ``atom`` style
    with ``dump_modify scale no`` (unscaled, wrapped coordinates).
    Coordinates are written with six decimals.

    Parameters
    ----------
    path : Path
        File to write.
    positions : numpy.ndarray, optional
        Shape ``(n_frames, n_atoms, 3)``; defaults to 3 frames of a 4-atom
        chain along x that moves 0.1 along z per frame.
    atom_types : sequence of int, optional
        Per-atom types from 1; defaults to all 1.
    timesteps : sequence of int, optional
        One per frame; defaults to ``0, 1000, 2000, ...``.
    box : sequence of (float, float)
        ``(lo, hi)`` for x, y, and z, the same in every frame.
    boundary : str
        Three two-letter boundary flags, each letter one of ``p``, ``f``,
        ``s``, or ``m``, such as ``'ff ff pp'`` for a cylinder periodic
        along z.

    Returns
    -------
    Path
        `path`.

    Raises
    ------
    ValueError
        If `positions` is not ``(n_frames, n_atoms, 3)``, `timesteps` does
        not have one entry per frame, or `boundary` is malformed.

    Notes
    -----
    Time and space are O(n_frames * n_atoms).

    Examples
    --------
    .. code-block:: python

        def test_example(make_lammps_dump):
            path = make_lammps_dump(boundary="ff ff pp")
    """
    if positions is None:
        shift = np.zeros((3, 1, 3))
        shift[:, 0, 2] = 0.1 * np.arange(3)
        pos = _default_chain(4)[np.newaxis] + shift
    else:
        pos = np.asarray(positions)
    if pos.ndim != 3 or pos.shape[2] != 3:
        raise ValueError(
            f"'positions' must have shape "
            f"(n_frames, n_atoms, 3), got {pos.shape}"
        )
    n_frames, n_atoms = pos.shape[:2]
    steps = (
        list(timesteps)
        if timesteps is not None
        else [1000 * f for f in range(n_frames)]
    )
    if len(steps) != n_frames:
        raise ValueError(
            f"'timesteps' has {len(steps)} entries for " f"{n_frames} frames"
        )
    flags = boundary.split()
    if len(flags) != 3 or any(
        len(f) != 2 or set(f) - set("pfsm") for f in flags
    ):
        raise ValueError(
            f"'boundary' must be three flags such as "
            f"'pp pp pp', got {boundary!r}"
        )
    types = list(atom_types) if atom_types is not None else [1] * n_atoms

    lines: list[str] = []
    for step, frame in zip(steps, pos):
        lines += [
            "ITEM: TIMESTEP",
            str(step),
            "ITEM: NUMBER OF ATOMS",
            str(n_atoms),
            f"ITEM: BOX BOUNDS {boundary}",
        ] + _box_lines(box)
        lines.append("ITEM: ATOMS id type x y z")
        lines += [
            f"{i} {atype} " + " ".join(f"{c:.6f}" for c in xyz)
            for i, (xyz, atype) in enumerate(zip(frame, types), 1)
        ]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def make_lammps_data(tmp_path: Path) -> Callable[..., Path]:
    """Factory: write a tiny LAMMPS data file under ``tmp_path``.

    Takes `name` (default ``'chain.data'``) and the keyword arguments of
    :func:`write_lammps_data`.
    """

    def _make(name: str = "chain.data", **kwargs: Any) -> Path:
        return write_lammps_data(tmp_path / name, **kwargs)

    return _make


@pytest.fixture
def make_lammps_dump(tmp_path: Path) -> Callable[..., Path]:
    """Factory: write a tiny LAMMPS dump file under ``tmp_path``.

    Takes `name` (default ``'chain.lammpstrj'``) and the keyword arguments
    of :func:`write_lammps_dump`.
    """

    def _make(name: str = "chain.lammpstrj", **kwargs: Any) -> Path:
        return write_lammps_dump(tmp_path / name, **kwargs)

    return _make
