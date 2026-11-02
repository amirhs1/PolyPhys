"""Tests for :mod:`polyphys.simulate.dump`.

Tier A tests read dumps written in the test, by the shared
``make_lammps_dump`` factory or as hand-written frames; each expected value is
a number the test writes into the file, or follows from it by the box
formulas of the LAMMPS ``Howto_triclinic`` page.

Tier B tests read the dumps in ``polyphys/tests/data/lammps_dump/``, written
by LAMMPS (22 Jul 2025 - Update 4) from the inputs beside them. Their
expected values are other LAMMPS output for the same frame, such as the
full-precision ``x y z`` that LAMMPS wrote next to the scaled ``xs ys zs``,
read with :func:`numpy.loadtxt`, not with the module under test.
"""

import logging
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pytest

from polyphys.simulate.core import Box
from polyphys.simulate.dump import (
    BrokenDumpError,
    DumpFrame,
    iter_dump,
    read_dump_frame,
)


@pytest.fixture
def dump_dir() -> Path:
    """Return the directory of the LAMMPS-written dump fixtures."""
    return Path(__file__).resolve().parents[1] / "data" / "lammps_dump"


@pytest.fixture
def frame_text() -> Callable[..., str]:
    """Return a factory for the text of one hand-written dump frame.

    Each argument replaces one part of a valid two-atom frame; `n_atoms`
    defaults to the number of `atoms` lines, and `before` holds lines that
    come before ``ITEM: TIMESTEP``.
    """
    def _make(
        timestep: str = "0",
        n_atoms: str | None = None,
        box_header: str = "ITEM: BOX BOUNDS pp pp pp",
        box: Sequence[str] = ("0.0 10.0", "0.0 10.0", "0.0 10.0"),
        columns: str = "id type x y z",
        atoms: Sequence[str] = ("1 1 1.0 5.0 5.0", "2 1 2.0 5.0 5.0"),
        before: Sequence[str] = (),
    ) -> str:
        count = str(len(atoms)) if n_atoms is None else n_atoms
        lines = [*before, "ITEM: TIMESTEP", timestep,
                 "ITEM: NUMBER OF ATOMS", count, box_header, *box,
                 f"ITEM: ATOMS {columns}".rstrip(), *atoms]
        return "\n".join(lines) + "\n"
    return _make


@pytest.fixture
def write_dump(tmp_path: Path) -> Callable[..., Path]:
    """Return a factory that writes dump text to a file under ``tmp_path``."""
    def _write(text: str, name: str = "frame.lammpstrj") -> Path:
        path = tmp_path / name
        path.write_text(text)
        return path
    return _write


@pytest.fixture
def lammps_atoms() -> Callable[[Path, int, int], np.ndarray]:
    """Return a reader of one frame's atom lines through numpy.loadtxt.

    It returns the ``n_atoms`` lines after the frame's ``ITEM: ATOMS`` line
    as a float array sorted by its first column (the atom ID): an
    independent reading of LAMMPS output, used as the expected value.
    """
    def _read(path: Path, frame: int, n_atoms: int) -> np.ndarray:
        lines = path.read_text().splitlines()
        starts = [i for i, line in enumerate(lines)
                  if line.startswith("ITEM: ATOMS")]
        first = starts[frame] + 1
        rows = np.loadtxt(lines[first:first + n_atoms], ndmin=2)
        return rows[np.argsort(rows[:, 0])]
    return _read


# Reading frames -------------------------------------------------------------


def test_iter_dump_reads_every_frame_header(make_lammps_dump):
    """Each frame gives its timestep, atom count, and column names."""
    path = make_lammps_dump(timesteps=[0, 1000, 2000])

    frames = list(iter_dump(path))

    assert [frame.timestep for frame in frames] == [0, 1000, 2000]
    assert [frame.n_atoms for frame in frames] == [4, 4, 4]
    assert frames[0].columns == ("id", "type", "x", "y", "z")


def test_iter_dump_reads_atom_values(make_lammps_dump):
    """The per-atom values equal the numbers written in the file."""
    positions = np.array([[[1.25, 2.5, 3.75], [4.0, 5.0, 6.0]],
                          [[1.5, 2.5, 3.75], [4.0, 5.5, 6.0]]])
    path = make_lammps_dump(positions=positions, atom_types=[2, 1])

    frames = list(iter_dump(path))

    np.testing.assert_array_equal(frames[1].data, [[1, 2, 1.5, 2.5, 3.75],
                                                   [2, 1, 4.0, 5.5, 6.0]])


def test_iter_dump_yields_a_frame_before_reading_the_next(
        frame_text, write_dump):
    """A frame is yielded before the next one is read, broken or not."""
    text = frame_text() + frame_text(timestep="not a step")
    frames = iter_dump(write_dump(text))

    first = next(frames)

    assert first.timestep == 0


def test_iter_dump_ignores_trailing_blank_lines(make_lammps_dump):
    """Blank lines at the end of the file end the iteration cleanly."""
    path = make_lammps_dump()
    path.write_text(path.read_text() + "\n\n")

    frames = list(iter_dump(path))

    assert len(frames) == 3


@pytest.mark.parametrize(
    ("sort_by_id", "expected"),
    [
        pytest.param(True, [1, 2, 3, 4, 5, 6, 7, 8], id="sorted"),
        pytest.param(False, [1, 2, 8, 3, 7, 4, 6, 5], id="file-order"),
    ],
)
def test_iter_dump_sort_by_id_orders_rows(dump_dir, sort_by_id, expected):
    """Rows follow atom ID when sorting and LAMMPS's own order otherwise."""
    path = dump_dir / "ring.bond.lammpstrj"

    frame = next(iter_dump(path, sort_by_id=sort_by_id))

    np.testing.assert_array_equal(frame.column("id"), expected)


def test_iter_dump_keeps_file_order_without_an_id_column(
        frame_text, write_dump):
    """A frame without an id column keeps its rows in file order."""
    text = frame_text(columns="type x y z",
                      atoms=("2 3.0 5.0 5.0", "1 1.0 5.0 5.0"))

    frame = next(iter_dump(write_dump(text)))

    np.testing.assert_array_equal(frame.column("type"), [2, 1])


def test_iter_dump_reads_units_and_time(dump_dir):
    """UNITS, written in the first frame only, and TIME reach every frame."""
    path = dump_dir / "format.xu.lammpstrj"

    frames = list(iter_dump(path))

    # in.ring: timestep 0.005 and a dump every 500 steps, so 2.5 tau apart.
    assert [frame.units for frame in frames] == ["lj", "lj", "lj"]
    assert [frame.time for frame in frames] == [0.0, 2.5, 5.0]


def test_iter_dump_reads_gzip_like_plain_text(dump_dir):
    """A gzipped dump gives the same frames as its plain-text twin."""
    plain = list(iter_dump(dump_dir / "ring.bond.lammpstrj"))

    gzipped = list(iter_dump(dump_dir / "format.x.lammpstrj.gz"))

    np.testing.assert_array_equal(gzipped[2].data, plain[2].data)


# Boxes ----------------------------------------------------------------------


def test_iter_dump_reads_an_orthogonal_box(make_lammps_dump):
    """An orthogonal box gives its bounds, zero tilt, and boundary flags."""
    path = make_lammps_dump(box=((0.0, 10.0), (-1.0, 9.0), (2.0, 12.0)),
                            boundary="ff ff pp")

    box = next(iter_dump(path)).box

    np.testing.assert_array_equal(box.lo, [0.0, -1.0, 2.0])
    np.testing.assert_array_equal(box.hi, [10.0, 9.0, 12.0])
    np.testing.assert_array_equal(box.tilt, [0.0, 0.0, 0.0])
    assert box.boundary == ("ff", "ff", "pp")
    assert not box.triclinic


@pytest.mark.parametrize(
    ("bound_lines", "tilt"),
    [
        # lo (0, 0, 0), hi (6, 6, 6): x bounds 0 + min(0, 1, 0.5, 1.5) and
        # 6 + max(...) = 7.5; y bounds 0 + min(0, -0.5) = -0.5 and 6.
        pytest.param(("0.0 7.5 1.0", "-0.5 6.0 0.5", "0.0 6.0 -0.5"),
                     [1.0, 0.5, -0.5], id="xy-positive"),
        # x bounds 0 + min(0, -1, 0.5, -0.5) = -1 and 6 + 0.5; y bounds 0
        # and 6 + 0.5.
        pytest.param(("-1.0 6.5 -1.0", "0.0 6.5 0.5", "0.0 6.0 0.5"),
                     [-1.0, 0.5, 0.5], id="xy-negative"),
    ],
)
def test_iter_dump_triclinic_box_removes_the_bounding_offsets(
        frame_text, write_dump, bound_lines, tilt):
    """A triclinic header gives the box's own lo and hi, not the bounds."""
    text = frame_text(box_header="ITEM: BOX BOUNDS xy xz yz pp pp pp",
                      box=bound_lines)

    box = next(iter_dump(write_dump(text))).box

    np.testing.assert_array_equal(box.lo, [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(box.hi, [6.0, 6.0, 6.0])
    np.testing.assert_array_equal(box.tilt, tilt)
    assert box.triclinic


# Columns and coordinates ----------------------------------------------------


def test_iter_dump_maps_column_names_to_their_columns(frame_text, write_dump):
    """Columns in any order are found by name."""
    text = frame_text(columns="x id z type y",
                      atoms=("1.5 7 3.5 2 2.5",))

    frame = next(iter_dump(write_dump(text)))

    assert frame.column("id").tolist() == [7]
    assert frame.column("type").tolist() == [2]
    np.testing.assert_array_equal(frame.positions(), [[1.5, 2.5, 3.5]])


@pytest.mark.parametrize("name", ["id", "type", "ix"])
def test_dump_frame_column_returns_integer_attributes_as_int64(
        dump_dir, name):
    """IDs, types, and image flags come back as int64."""
    frame = next(iter_dump(dump_dir / "ring.bond.lammpstrj"))

    values = frame.column(name)

    assert values.dtype == np.int64


def test_dump_frame_column_returns_other_attributes_as_float64(
        frame_text, write_dump):
    """A non-integer column comes back as its float64 values."""
    frame = next(iter_dump(write_dump(frame_text())))

    values = frame.column("x")

    assert values.dtype == np.float64
    np.testing.assert_array_equal(values, [1.0, 2.0])


def test_dump_frame_column_raises_for_a_missing_name(make_lammps_dump):
    """Asking for a column the frame lacks raises KeyError."""
    frame = next(iter_dump(make_lammps_dump()))

    with pytest.raises(KeyError, match="no column 'vx'"):
        frame.column("vx")


@pytest.mark.parametrize(
    ("scaled", "reference", "atol"),
    [
        # dump atom writes %g, 6 significant digits: a scaled value below 1
        # is within 5e-7, times the box length 6.
        pytest.param("format.atom.lammpstrj", "ring.bond.lammpstrj", 3.0e-6,
                     id="orthogonal"),
        # Triclinic: x adds the errors of lx = 6, xy = 1, and xz = 0.5.
        pytest.param("triclinic.atom.lammpstrj", "triclinic.x.lammpstrj",
                     3.75e-6, id="triclinic"),
    ],
)
def test_positions_unscales_scaled_coordinates(
        dump_dir, lammps_atoms, scaled, reference, atol):
    """Scaled xs ys zs unscale to the x y z that LAMMPS wrote."""
    frame = read_dump_frame(dump_dir / scaled, 2)
    expected = lammps_atoms(dump_dir / reference, 2, 8)[:, -6:-3]

    positions = frame.positions()

    np.testing.assert_allclose(positions, expected, rtol=0, atol=atol)


@pytest.mark.parametrize(
    ("dump", "reference"),
    [
        pytest.param("ring.bond.lammpstrj", "format.xu.lammpstrj",
                     id="x-and-image-flags"),
        pytest.param("format.xsu.lammpstrj", "format.xu.lammpstrj",
                     id="xsu"),
        pytest.param("triclinic.x.lammpstrj", "triclinic.xu.lammpstrj",
                     id="triclinic-x-and-image-flags"),
        pytest.param("triclinic.xsu.lammpstrj", "triclinic.xu.lammpstrj",
                     id="triclinic-xsu"),
    ],
)
def test_positions_unwraps_to_lammps_unwrapped_coordinates(
        dump_dir, lammps_atoms, dump, reference):
    """Unwrapped positions from image flags or xsu equal LAMMPS's xu yu zu."""
    frame = read_dump_frame(dump_dir / dump, 0)
    expected = lammps_atoms(dump_dir / reference, 0, 8)[:, -3:]

    positions = frame.positions(unwrapped=True)

    # %.16g keeps 16 significant digits of coordinates below 20; LAMMPS
    # computes xu and xsu in float64 separately, so 1e-12 bounds both.
    np.testing.assert_allclose(positions, expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("unwrapped", [False, True])
def test_positions_raises_without_coordinate_columns(
        frame_text, write_dump, unwrapped):
    """A frame with no coordinate columns raises KeyError."""
    text = frame_text(columns="id type", atoms=("1 1", "2 1"))
    frame = next(iter_dump(write_dump(text)))

    with pytest.raises(KeyError, match="no coordinate columns"):
        frame.positions(unwrapped=unwrapped)


def test_iter_dump_reads_a_frame_without_atoms(frame_text, write_dump):
    """A frame of zero atoms gives an empty (0, n_columns) array."""
    text = frame_text(atoms=())

    frame = next(iter_dump(write_dump(text)))

    assert frame.data.shape == (0, 5)


# Broken files ---------------------------------------------------------------


def test_iter_dump_raises_on_a_cut_off_last_frame(make_lammps_dump):
    """A last frame short of atom lines raises, naming that frame."""
    path = make_lammps_dump(timesteps=[0, 1000, 2000])
    path.write_text("\n".join(path.read_text().splitlines()[:-2]) + "\n")

    with pytest.raises(BrokenDumpError,
                       match=r"frame 2 \(timestep 2000\), line 37: 2 of 4"):
        list(iter_dump(path))


def test_iter_dump_raises_when_a_header_is_cut_off(frame_text, write_dump):
    """A file that ends inside a frame header raises."""
    text = frame_text() + "ITEM: TIMESTEP\n100\n"

    with pytest.raises(BrokenDumpError, match="ends inside the frame header"):
        list(iter_dump(write_dump(text)))


@pytest.mark.parametrize(
    ("parts", "match"),
    [
        pytest.param({"before": ("garbage",)}, "expected an ITEM: line",
                     id="not-an-item"),
        pytest.param({"before": ("ITEM: FOO",)}, "expected ITEM: TIMESTEP",
                     id="unknown-item"),
        pytest.param({"before": ("ITEM: TIME", "soon")},
                     "TIME is not followed by a number", id="time-not-number"),
        pytest.param({"timestep": "1e3"}, "TIMESTEP is not followed by an "
                     "integer", id="timestep-not-integer"),
        pytest.param({"timestep": ""}, "TIMESTEP is not followed by a value",
                     id="timestep-empty"),
        pytest.param({"n_atoms": "two"}, "NUMBER OF ATOMS is not followed",
                     id="count-not-integer"),
        pytest.param({"n_atoms": "-1"}, "negative number of atoms",
                     id="count-negative"),
        pytest.param({"box_header": "ITEM: BOX"}, "expected ITEM: BOX BOUNDS",
                     id="no-box-bounds"),
        pytest.param({"box_header": "ITEM: BOX BOUNDS pp pp"},
                     "three boundary flags", id="two-flags"),
        pytest.param({"box": ("0.0", "0.0 10.0", "0.0 10.0")},
                     "needs 2 numbers", id="box-line-short"),
        pytest.param({"box": ("a b", "0.0 10.0", "0.0 10.0")},
                     "needs 2 numbers", id="box-line-not-numbers"),
        pytest.param({"box": ("10.0 0.0", "0.0 10.0", "0.0 10.0")},
                     "must be below", id="box-lo-above-hi"),
        pytest.param({"columns": ""}, "names no columns", id="no-columns"),
        pytest.param({"columns": "id id x y z"}, "repeats",
                     id="repeated-column"),
        pytest.param({"atoms": ("1 1 1.0 5.0", "2 1 2.0 5.0 5.0")},
                     r"line 10: an atom line needs 5 numbers, got '1 1 1.0 "
                     r"5.0'", id="atom-line-short"),
        pytest.param({"atoms": ("1 1 1.0 5.0 5.0", "2 1 2.0 5.0 abc")},
                     r"line 11: an atom line needs 5 numbers",
                     id="atom-value-not-number"),
    ],
)
def test_iter_dump_raises_on_a_malformed_frame(
        frame_text, write_dump, parts, match):
    """A frame that breaks the format raises BrokenDumpError."""
    path = write_dump(frame_text(**parts))

    with pytest.raises(BrokenDumpError, match=match):
        list(iter_dump(path))


def test_iter_dump_rejects_general_triclinic_boxes(frame_text, write_dump):
    """The general triclinic box form raises BrokenDumpError."""
    text = frame_text(box_header="ITEM: BOX BOUNDS abc origin pp pp pp",
                      box=("6 0 0 0", "1 6 0 0", "0.5 -0.5 6 0"))

    with pytest.raises(BrokenDumpError, match="general triclinic"):
        list(iter_dump(write_dump(text)))


def test_iter_dump_warns_when_a_timestep_repeats(make_lammps_dump, caplog):
    """A frame whose timestep does not increase is logged as a warning."""
    path = make_lammps_dump(timesteps=[0, 1000, 1000])

    with caplog.at_level(logging.WARNING,
                         logger="polyphys.simulate.dump"):
        frames = list(iter_dump(path))

    assert len(frames) == 3
    assert "frame 2 has timestep 1000, not after the previous frame's 1000" \
        in caplog.text


# read_dump_frame ------------------------------------------------------------


@pytest.mark.parametrize("index", [0, 2])
def test_read_dump_frame_returns_the_frame_at_index(make_lammps_dump, index):
    """read_dump_frame returns the frame at the given position."""
    path = make_lammps_dump(timesteps=[0, 1000, 2000])

    frame = read_dump_frame(path, index)

    assert frame.timestep == 1000 * index


@pytest.mark.parametrize(
    ("index", "match"),
    [
        pytest.param(3, "has 3 frames; there is no frame 3", id="past-end"),
        pytest.param(-1, "must be 0 or more", id="negative"),
    ],
)
def test_read_dump_frame_raises_for_a_missing_index(
        make_lammps_dump, index, match):
    """An index outside the file raises IndexError."""
    path = make_lammps_dump()

    with pytest.raises(IndexError, match=match):
        read_dump_frame(path, index)


# Records --------------------------------------------------------------------


def test_dump_frame_data_is_read_only(make_lammps_dump):
    """A frame's data cannot be changed in place."""
    frame = next(iter_dump(make_lammps_dump()))

    with pytest.raises(ValueError, match="read-only"):
        frame.data[0, 0] = 99.0


@pytest.mark.parametrize(
    ("columns", "data", "match"),
    [
        pytest.param(("id", "x"), np.zeros((2, 3)), "shape", id="width"),
        pytest.param(("id", "id"), np.zeros((2, 2)), "repeat", id="repeat"),
    ],
)
def test_dump_frame_rejects_inconsistent_data(columns, data, match):
    """Data whose width does not match the column names raises ValueError."""
    box = Box(lo=np.zeros(3), hi=np.ones(3))

    with pytest.raises(ValueError, match=match):
        DumpFrame(timestep=0, box=box, columns=columns, data=data)
