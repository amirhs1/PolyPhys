"""Tests for :mod:`polyphys.simulate.data`.

Tier A tests read data files written in the test, by the shared
``make_lammps_data`` factory or as hand-written text, and write records built
in the test; each expected value is a number the test writes, or the record
itself after a round trip.

Tier B tests read the data files and dumps in
``polyphys/tests/data/lammps_dump/``, written by LAMMPS (22 Jul 2025 -
Update 4) from the inputs beside them. Their expected values come from the
LAMMPS input scripts (box, tilt, masses, charges, ring topology) or from the
file text read with :func:`numpy.loadtxt`, not with the module under test.
"""

import dataclasses
import gzip
import logging
from operator import attrgetter
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np
import pytest

from polyphys.simulate.core import Box
from polyphys.simulate.data import (
    BrokenDataError,
    DataFile,
    read_data,
    with_frame,
    write_data,
)
from polyphys.simulate.dump import DumpFrame, read_dump_frame


@pytest.fixture
def data_dir() -> Path:
    """Return the directory of the LAMMPS-written fixtures."""
    return Path(__file__).resolve().parents[1] / "data" / "lammps_dump"


@pytest.fixture
def section_rows() -> Callable[[Path, str, int], np.ndarray]:
    """Return a reader of one data-file section through numpy.loadtxt.

    It returns the `count` lines after the line that starts with `keyword`
    and the blank line below it, as a float array sorted by its first
    column: an independent reading of LAMMPS output.
    """
    def _read(path: Path, keyword: str, count: int) -> np.ndarray:
        lines = path.read_text().splitlines()
        first = next(i for i, line in enumerate(lines)
                     if line.startswith(keyword)) + 2
        rows = np.loadtxt(lines[first:first + count], ndmin=2)
        return rows[np.argsort(rows[:, 0])]
    return _read


@pytest.fixture
def make_record() -> Callable[..., DataFile]:
    """Return a factory of small records, one per atom style.

    Three atoms in a box of side 10, with floats that only an exact format
    keeps (1/3, 1e-17); molecule IDs, charges, two bonds, and one angle as
    the style allows. Keyword arguments replace fields.
    """
    def _make(atom_style: str = "bond", **overrides: Any) -> DataFile:
        fields: dict[str, Any] = {
            "atom_style": atom_style,
            "box": Box(lo=np.zeros(3), hi=np.full(3, 10.0)),
            "ids": np.array([1, 2, 3]),
            "types": np.array([1, 1, 2]),
            "positions": np.array([[0.1, 1 / 3, 2.0], [1e-17, 2.5e-5, 3.0],
                                   [9.75, 5.0, 7.125]]),
            "n_atom_types": 2,
        }
        if atom_style != "atomic":
            fields.update(mol_ids=np.array([1, 1, 2]),
                          bonds=np.array([[1, 1, 2], [1, 2, 3]]),
                          n_bond_types=1)
        if atom_style in ("angle", "molecular", "full"):
            fields.update(angles=np.array([[1, 1, 2, 3]]), n_angle_types=1)
        if atom_style == "full":
            fields["charges"] = np.array([0.25, -0.5, 0.25])
        fields.update(overrides)
        return DataFile(**fields)
    return _make


@pytest.fixture
def data_text() -> Callable[..., str]:
    """Return a factory for the text of a hand-written data file.

    The default is two bonded atoms in atom style ``bond``; `header` and
    `sections` (pairs of a keyword line and its entry lines) replace it.
    """
    def _make(
        header: Sequence[str] = ("2 atoms", "1 atom types", "1 bonds",
                                 "1 bond types", "0 5 xlo xhi", "0 5 ylo yhi",
                                 "0 5 zlo zhi"),
        sections: Sequence[tuple[str, Sequence[str]]] = (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)),
        ),
    ) -> str:
        lines = ["hand-written data file", "", *header, ""]
        for keyword, rows in sections:
            lines += [keyword, "", *rows, ""]
        return "\n".join(lines)
    return _make


@pytest.fixture
def write_text(tmp_path: Path) -> Callable[..., Path]:
    """Return a factory that writes text to a file under ``tmp_path``."""
    def _write(text: str, name: str = "system.data") -> Path:
        path = tmp_path / name
        path.write_text(text)
        return path
    return _write


@pytest.fixture
def record_fields() -> Callable[[DataFile], dict[str, Any]]:
    """Return a function that lists a record's fields as plain values."""
    def _plain(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, Box):
            return {name: _plain(getattr(value, name))
                    for name in ("lo", "hi", "tilt", "boundary", "triclinic")}
        return value

    def _fields(record: DataFile) -> dict[str, Any]:
        return {field.name: _plain(getattr(record, field.name))
                for field in dataclasses.fields(record)}
    return _fields


# read_data ------------------------------------------------------------------


@pytest.mark.parametrize(
    ("atom_style", "bonds", "angles", "mol_ids"),
    [
        pytest.param("atomic", [], [], None, id="atomic"),
        pytest.param("bond", [[1, 1, 2], [1, 2, 3], [1, 3, 4]], [],
                     [1, 1, 2, 2], id="bond"),
        pytest.param("angle", [[1, 1, 2], [1, 2, 3], [1, 3, 4]],
                     [[1, 1, 2, 3], [1, 2, 3, 4]], [1, 1, 2, 2], id="angle"),
    ],
)
def test_read_data_reads_the_atoms_and_topology_written(
        make_lammps_data, record_fields, atom_style, bonds, angles, mol_ids):
    """Positions, types, molecule IDs, bonds, and angles read as written."""
    positions = np.array([[1.5, 2.0, 2.5], [2.5, 2.0, 2.5],
                          [3.5, 2.0, 2.5], [4.5, 2.0, 2.5]])
    path = make_lammps_data(atom_style=atom_style, positions=positions,
                            atom_types=[1, 2, 1, 1], mol_ids=[1, 1, 2, 2])

    fields = record_fields(read_data(path))

    assert fields["positions"] == positions.tolist()
    assert fields["types"] == [1, 2, 1, 1]
    assert fields["mol_ids"] == mol_ids
    assert fields["bonds"] == bonds
    assert fields["angles"] == angles
    assert fields["masses"] == [1.0, 1.0]


@pytest.mark.parametrize(
    ("name", "n_atoms", "x_column"),
    [
        pytest.param("atomic.data", 2, 2, id="atomic"),
        pytest.param("ring.bond.data", 8, 3, id="bond"),
        pytest.param("ring.angle.data", 8, 3, id="angle"),
        pytest.param("ring.molecular.data", 8, 3, id="molecular"),
        pytest.param("ring.full.data", 8, 4, id="full"),
        pytest.param("triclinic.data", 8, 3, id="triclinic"),
    ],
)
def test_read_data_reads_lammps_atoms_and_image_flags(
        data_dir, section_rows, name, n_atoms, x_column):
    """Atoms of a LAMMPS write_data file read back, sorted by atom ID."""
    path = data_dir / name
    rows = section_rows(path, "Atoms", n_atoms)

    data = read_data(path)

    np.testing.assert_array_equal(data.ids, rows[:, 0])
    np.testing.assert_array_equal(data.positions,
                                  rows[:, x_column:x_column + 3])
    np.testing.assert_array_equal(data.image_flags, rows[:, -3:])


def test_read_data_reads_a_gzipped_file(make_lammps_data, tmp_path):
    """A data file whose name ends in .gz reads like its plain text."""
    plain = make_lammps_data(atom_style="bond")
    gzipped = tmp_path / "chain.data.gz"
    with gzip.open(gzipped, "wt") as stream:
        stream.write(plain.read_text())

    data = read_data(gzipped)

    np.testing.assert_array_equal(data.positions, read_data(plain).positions)


def test_read_data_orders_velocities_by_atom_id(data_dir, section_rows):
    """Velocities, listed by LAMMPS in its own atom order, follow atom ID."""
    path = data_dir / "ring.bond.data"
    rows = section_rows(path, "Velocities", 8)

    data = read_data(path)

    np.testing.assert_array_equal(data.velocities, rows[:, 1:])


@pytest.mark.parametrize(
    ("name", "field", "expected"),
    [
        # in.atomic: mass 1 1.0 and mass 2 2.0.
        pytest.param("atomic.data", "masses", [1.0, 2.0], id="masses"),
        # in.ring with style full: set type 1 charge 0.25.
        pytest.param("ring.full.data", "charges", [0.25] * 8, id="charges"),
        # in.triclinic: region prism 0 6 0 6 0 6 1.0 0.5 -0.5.
        pytest.param("triclinic.data", "box.tilt", [1.0, 0.5, -0.5],
                     id="tilt"),
        pytest.param("triclinic.data", "box.hi", [6.0, 6.0, 6.0], id="hi"),
    ],
)
def test_read_data_reads_values_set_by_the_input_script(
        data_dir, name, field, expected):
    """Masses, charges, and the box equal those the LAMMPS script set."""
    data = read_data(data_dir / name)

    np.testing.assert_array_equal(attrgetter(field)(data), expected)


@pytest.mark.parametrize(
    ("name", "field", "expected"),
    [
        # ring8.mol: bond i joins beads i and i + 1, and bead 8 joins bead 1.
        pytest.param("ring.bond.data", "bonds",
                     {(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8),
                      (8, 1)}, id="bonds"),
        # ring8_angles.mol: angle i is beads i, i + 1, i + 2, cyclically.
        pytest.param("ring.angle.data", "angles",
                     {(1, 2, 3), (2, 3, 4), (3, 4, 5), (4, 5, 6), (5, 6, 7),
                      (6, 7, 8), (7, 8, 1), (8, 1, 2)}, id="angles"),
    ],
)
def test_read_data_reads_the_ring_topology(data_dir, name, field, expected):
    """The bonds and angles are those of the ring molecule LAMMPS built."""
    data = read_data(data_dir / name)

    members = {tuple(row) for row in getattr(data, field)[:, 1:].tolist()}

    assert members == expected


def test_read_data_accepts_any_header_order_spacing_and_comments(
        data_text, write_text):
    """Header lines in any order, runs of blanks, and comments are read."""
    text = data_text(
        header=("0  5\tzlo zhi", "1 bond types  # one", "0 5 xlo   xhi",
                "2 atoms", "1 bonds", "0 5 ylo yhi", "1  atom types",
                "2 extra bond per atom"),
        sections=(("Atoms  # bond", ("2 1 1  3.0 2.5 2.5  # second",
                                     "1\t1 1 2.0 2.5 2.5")),
                  ("Bonds", ("1 1 1 2",))),
    )

    data = read_data(write_text(text))

    np.testing.assert_array_equal(data.positions, [[2.0, 2.5, 2.5],
                                                   [3.0, 2.5, 2.5]])
    assert data.bonds.tolist() == [[1, 1, 2]]
    np.testing.assert_array_equal(data.box.hi, [5.0, 5.0, 5.0])


def test_read_data_takes_the_style_from_its_argument(data_text, write_text):
    """An Atoms line without a style comment uses the style passed in."""
    text = data_text(sections=(
        ("Atoms", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
        ("Bonds", ("1 1 1 2",)),
    ))

    data = read_data(write_text(text), atom_style="bond")

    assert data.atom_style == "bond"


@pytest.mark.parametrize(
    ("atoms", "atom_style", "match"),
    [
        pytest.param("Atoms", None, "names no atom style", id="no-style"),
        pytest.param("Atoms # bond", "angle", "'angle' but the file says "
                     "'bond'", id="disagree"),
        pytest.param("Atoms # sphere", None, "'sphere' is not supported",
                     id="unsupported"),
    ],
)
def test_read_data_raises_without_one_supported_atom_style(
        data_text, write_text, atoms, atom_style, match):
    """A missing, conflicting, or unsupported atom style raises."""
    text = data_text(sections=(
        (atoms, ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
        ("Bonds", ("1 1 1 2",)),
    ))

    with pytest.raises(BrokenDataError, match=match):
        read_data(write_text(text), atom_style=atom_style)


def test_read_data_skips_coefficient_sections(data_text, write_text, caplog):
    """Sections whose names end in Coeffs are skipped and logged."""
    text = data_text(sections=(
        ("Pair Coeffs # lj/cut", ("1 1.0 1.0",)),
        ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
        ("Bond Coeffs # fene", ("1 30.0 1.5 1.0 1.0",)),
        ("Bonds", ("1 1 1 2",)),
    ))

    with caplog.at_level(logging.DEBUG, logger="polyphys.simulate.data"):
        data = read_data(write_text(text))

    assert data.n_bonds == 1
    assert "skipped section 'Pair Coeffs'" in caplog.text
    assert "skipped section 'Bond Coeffs'" in caplog.text


@pytest.mark.parametrize(
    ("parts", "match"),
    [
        pytest.param({"header": ("2 atoms", "1 atom types", "1 dihedrals",
                                 "0 5 xlo xhi", "0 5 ylo yhi",
                                 "0 5 zlo zhi")},
                     "dihedrals are not supported", id="dihedrals"),
        pytest.param({"header": ("2 atoms", "1 atom types", "1 bonds",
                                 "1 bond types", "6 0 0 avec", "0 5 xlo xhi",
                                 "0 5 ylo yhi", "0 5 zlo zhi")},
                     "general triclinic", id="general-triclinic"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)), ("Ellipsoids", ("1 1.0 1.0 1.0",)))},
            "section 'Ellipsoids' is not supported", id="ellipsoids"),
    ],
)
def test_read_data_rejects_unsupported_features(
        data_text, write_text, parts, match):
    """Dihedrals, general triclinic boxes, and other sections raise."""
    path = write_text(data_text(**parts))

    with pytest.raises(BrokenDataError, match=match):
        read_data(path)


@pytest.mark.parametrize(
    ("parts", "match"),
    [
        pytest.param({"sections": (("Bonds", ("1 1 1 2",)),)},
                     "no Atoms section", id="no-atoms"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5",)),
            ("Bonds", ("1 1 1 2",)))},
            "Atoms has 1 entries, the header says 2", id="atom-count"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),)},
            "1 bonds but no Bonds", id="no-bonds-section"),
        pytest.param({"header": ("2 atoms", "1 atom types", "0 bonds",
                                 "1 bond types", "1 angles", "1 angle types",
                                 "0 5 xlo xhi", "0 5 ylo yhi", "0 5 zlo zhi"),
                      "sections": (("Atoms # angle", ("1 1 1 2.0 2.5 2.5",
                                                      "2 1 1 3.0 2.5 2.5")),)},
                     "1 angles but no Angles", id="no-angles-section"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 3",)))},
            r"system.data: a bond names atom IDs that do not exist: \[3\]",
            id="bond-to-missing-atom"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5 0 0 0",
                              "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)))},
            "the same on every line", id="image-flags-on-one-line"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)), ("Bonds", ("1 1 1 2",)))},
            "section 'Bonds' repeats", id="repeated-section"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)),
            ("Velocities", ("1 0 0 0", "3 0 0 0")))},
            "Velocities atom IDs differ", id="velocity-ids"),
        pytest.param({"sections": (
            ("Masses", ("2 1.0",)),
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2",)))},
            "each atom type from 1 to 1 once", id="masses"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
            ("Bonds", ("1 1 1 2.5",)))},
            "'Bonds' must hold integers", id="bond-not-integer"),
        pytest.param({"header": ("2 atoms", "1 atom types", "1 bonds",
                                 "1 bond types", "0 5 xlo xhi",
                                 "0 5 ylo yhi")},
                     r"no \['zlo'\] bounds", id="no-z-bounds"),
        pytest.param({"header": ("2 atoms", "1 atom types", "1 bonds",
                                 "1 bond types", "5 0 xlo xhi",
                                 "0 5 ylo yhi", "0 5 zlo zhi")},
                     "must be below", id="xlo-above-xhi"),
        pytest.param({"header": ("two atoms", "1 atom types")},
                     "invalid literal for int", id="count-not-integer"),
        pytest.param({"header": ("2 atoms", "1 atom types", "1 bonds",
                                 "1 bond types", "0 5 xlo xhi",
                                 "0 5 ylo yhi", "0 5 zlo zhi", "1.0 atoms")},
                     "invalid literal", id="count-not-whole"),
        pytest.param({"sections": (
            ("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 abc")),
            ("Bonds", ("1 1 1 2",)))},
            "not a number", id="atom-value-not-number"),
    ],
)
def test_read_data_raises_on_a_malformed_file(
        data_text, write_text, parts, match):
    """A file that breaks the format raises BrokenDataError naming it."""
    path = write_text(data_text(**parts))

    with pytest.raises(BrokenDataError, match=match):
        read_data(path)


def test_read_data_reads_an_empty_section_as_none(data_text, write_text):
    """A section keyword with no entries, as for zero bonds, reads as empty."""
    text = data_text(
        header=("2 atoms", "1 atom types", "0 bonds", "0 bond types",
                "0 5 xlo xhi", "0 5 ylo yhi", "0 5 zlo zhi"),
        sections=(("Atoms # bond", ("1 1 1 2.0 2.5 2.5", "2 1 1 3.0 2.5 2.5")),
                  ("Bonds", ())))

    data = read_data(write_text(text))

    assert data.bonds.shape == (0, 3)


@pytest.mark.parametrize(
    ("text", "match"),
    [
        pytest.param("title\n\n2 atoms\n0 5 xlo xhi\n", "no body",
                     id="header-only"),
        pytest.param("title\n\n2 atoms\n0 5 xlo xhi\n0 5 ylo yhi\n"
                     "0 5 zlo zhi\n1 2 3 4 5\n",
                     "an entry before any section", id="entry-first"),
    ],
)
def test_read_data_raises_on_a_file_without_a_body(write_text, text, match):
    """A file with no section keyword raises BrokenDataError."""
    with pytest.raises(BrokenDataError, match=match):
        read_data(write_text(text))


# write_data -----------------------------------------------------------------


def test_write_data_declares_bonds_with_their_own_keywords(
        make_record, tmp_path):
    """The header counts bonds as bonds and bond types, not as atoms."""
    record = make_record("bond")

    text = write_data(record, tmp_path / "out.data").read_text()

    assert text.splitlines()[2:6] == [
        "3 atoms", "2 atom types", "2 bonds", "1 bond types"]


def test_write_data_counts_every_bond_of_a_ring(make_record, tmp_path):
    """A ring of four atoms declares and lists its four bonds."""
    ring = make_record(
        "bond", ids=np.arange(1, 5), types=np.ones(4, dtype=int),
        mol_ids=np.ones(4, dtype=int), n_atom_types=1,
        positions=np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0],
                            [2.0, 2.0, 1.0], [1.0, 2.0, 1.0]]),
        bonds=np.array([[1, 1, 2], [1, 2, 3], [1, 3, 4], [1, 4, 1]]))

    lines = write_data(ring, tmp_path / "ring.data").read_text().splitlines()

    assert "4 bonds" in lines
    assert lines[lines.index("Bonds") + 2:] == [
        "1 1 1 2", "2 1 2 3", "3 1 3 4", "4 1 4 1"]


@pytest.mark.parametrize(
    "atom_style", ["atomic", "bond", "angle", "molecular", "full"])
def test_write_data_round_trips_exactly(
        make_record, record_fields, tmp_path, atom_style):
    """Writing and reading back returns every field unchanged."""
    record = make_record(
        atom_style, masses=np.array([1.0, 1 / 7]),
        image_flags=np.array([[0, 0, 0], [-1, 2, 0], [0, 0, 1]]),
        velocities=np.array([[0.1, -0.2, 1e-300], [0.0, 1.0, 2.0],
                             [3.5, -1 / 3, 0.0]]))

    back = read_data(write_data(record, tmp_path / "out.data"))

    assert record_fields(back) == record_fields(record)


def test_write_data_round_trips_a_ring_with_angles(
        make_record, record_fields, tmp_path):
    """A ring's closing bond and two wrap-around angles read back."""
    t = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    ring = make_record(
        "angle", ids=np.arange(1, 7), types=np.ones(6, dtype=int),
        mol_ids=np.ones(6, dtype=int), n_atom_types=1,
        positions=np.column_stack([5 + np.cos(t), 5 + np.sin(t),
                                   np.full(6, 5.0)]),
        bonds=np.array([[1, 1, 2], [1, 2, 3], [1, 3, 4], [1, 4, 5],
                        [1, 5, 6], [1, 6, 1]]),
        angles=np.array([[1, 1, 2, 3], [1, 2, 3, 4], [1, 3, 4, 5],
                         [1, 4, 5, 6], [1, 5, 6, 1], [1, 6, 1, 2]]))

    back = read_data(write_data(ring, tmp_path / "ring.data"))

    assert record_fields(back) == record_fields(ring)


@pytest.mark.parametrize(
    ("image_flags", "n_fields"),
    [
        pytest.param(np.array([[1, 0, -1], [0, 0, 0], [0, 0, 0]]), 9,
                     id="given"),
        pytest.param(None, 6, id="none"),
    ],
)
def test_write_data_writes_image_flags_only_when_given(
        make_record, tmp_path, image_flags, n_fields):
    """Atom lines end with image flags only when the record has them."""
    record = make_record("bond", image_flags=image_flags)

    lines = write_data(record, tmp_path / "out.data").read_text().splitlines()

    assert len(lines[lines.index("Atoms # bond") + 2].split()) == n_fields


@pytest.mark.parametrize(
    ("velocities", "written"),
    [
        pytest.param(np.zeros((3, 3)), True, id="given"),
        pytest.param(None, False, id="none"),
    ],
)
def test_write_data_writes_velocities_only_when_given(
        make_record, tmp_path, velocities, written):
    """The Velocities section appears only when the record has velocities."""
    record = make_record("bond", velocities=velocities)

    lines = write_data(record, tmp_path / "out.data").read_text().splitlines()

    assert ("Velocities" in lines) is written


def test_write_data_writes_no_masses_when_none(make_record, tmp_path):
    """Without masses there is no Masses section; the type count stays."""
    record = make_record("bond", types=np.ones(3, dtype=int), masses=None)

    lines = write_data(record, tmp_path / "out.data").read_text().splitlines()

    assert "Masses" not in lines
    assert "2 atom types" in lines


def test_write_data_writes_and_reads_tilt_factors(make_record, tmp_path):
    """A triclinic record gets an xy xz yz line that reads back exactly."""
    record = make_record("bond", box=Box(
        lo=np.zeros(3), hi=np.full(3, 10.0), tilt=np.array([1.0, 0.5, -0.5]),
        triclinic=True))

    path = write_data(record, tmp_path / "out.data")

    assert "1.0 0.5 -0.5 xy xz yz" in path.read_text().splitlines()
    np.testing.assert_array_equal(read_data(path).box.tilt,
                                  [1.0, 0.5, -0.5])


def test_write_data_rejects_a_title_with_a_line_break(make_record, tmp_path):
    """A title of more than one line raises ValueError."""
    with pytest.raises(ValueError, match="one line"):
        write_data(make_record(), tmp_path / "out.data", title="a\nb")


# DataFile -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("atom_style", "overrides", "match"),
    [
        pytest.param("sphere", {}, "'atom_style' must be one of",
                     id="unknown-style"),
        pytest.param("bond", {"positions": np.zeros((2, 3))},
                     r"'positions' must have shape \(3, 3\)",
                     id="positions-shape"),
        pytest.param("bond", {"bonds": np.array([[1, 1, 4]])},
                     r"names atom IDs that do not exist: \[4\]",
                     id="bond-to-missing-atom"),
        pytest.param("angle", {"angles": np.array([[1, 1, 2, 9]])},
                     r"names atom IDs that do not exist: \[9\]",
                     id="angle-to-missing-atom"),
        pytest.param("bond", {"ids": np.array([1, 1, 2])}, "unique",
                     id="repeated-id"),
        pytest.param("bond", {"ids": np.array([0, 1, 2]),
                              "bonds": np.array([[1, 1, 2]])},
                     "1 or more", id="id-zero"),
        pytest.param("bond", {"types": np.array([1, 2, 3])},
                     "atom types must be from 1 to 2", id="type-range"),
        pytest.param("bond", {"types": np.array([1.0, 1.5, 2.0])},
                     "'types' must hold integers", id="type-not-integer"),
        pytest.param("bond", {"types": np.array(["a", "b", "c"])},
                     "'types' must hold integers", id="type-not-number"),
        pytest.param("bond", {"mol_ids": None}, "needs 'mol_ids'",
                     id="no-mol-ids"),
        pytest.param("atomic", {"mol_ids": np.ones(3, dtype=int)},
                     "has no 'mol_ids'", id="mol-ids-in-atomic"),
        pytest.param("full", {"charges": None}, "needs 'charges'",
                     id="no-charges"),
        pytest.param("atomic", {"bonds": np.array([[1, 1, 2]]),
                                "n_bond_types": 1}, "has no bonds",
                     id="bonds-in-atomic"),
        pytest.param("bond", {"angles": np.array([[1, 1, 2, 3]]),
                              "n_angle_types": 1}, "has no angles",
                     id="angles-in-bond"),
        pytest.param("bond", {"n_bond_types": 0},
                     "'n_bond_types' must be 1 or more", id="no-bond-types"),
        pytest.param("bond", {"masses": np.array([1.0])},
                     r"'masses' must have shape \(2,\)", id="masses-shape"),
    ],
)
def test_data_file_rejects_an_invalid_record(
        make_record, atom_style, overrides, match):
    """A record that LAMMPS could not read raises ValueError."""
    with pytest.raises(ValueError, match=match):
        make_record(atom_style, **overrides)


def test_data_file_rejects_a_box_that_is_not_a_box(make_record):
    """A box given as bounds instead of a Box raises TypeError."""
    with pytest.raises(TypeError, match="'box' must be a Box"):
        make_record("bond", box=(np.zeros(3), np.full(3, 10.0)))


@pytest.mark.parametrize(
    ("image_flags", "expected_x"),
    [
        # Box of side 10: each image flag along x moves the atom by 10.
        pytest.param(np.array([[-1, 0, 0], [0, 0, 0], [2, 0, 0]]),
                     [-9.9, 0.0, 29.75], id="image-flags"),
        pytest.param(None, [0.1, 1e-17, 9.75], id="none"),
    ],
)
def test_data_file_unwraps_positions_by_image_flags(
        make_record, image_flags, expected_x):
    """Unwrapped x moves by the box length per image flag, or not at all."""
    record = make_record("bond", image_flags=image_flags)

    unwrapped = record.unwrapped_positions()

    # One addition of a multiple of 10 to each x: float64 round-off.
    np.testing.assert_allclose(unwrapped[:, 0], expected_x, rtol=0,
                               atol=1e-12)


def test_data_file_accepts_empty_bond_lists(make_record):
    """An empty list of bonds is a record with no bonds."""
    record = make_record("bond", bonds=[], n_bond_types=0)

    assert record.bonds.shape == (0, 3)


def test_data_file_arrays_are_read_only(make_record):
    """A record's arrays cannot be changed in place."""
    record = make_record()

    with pytest.raises(ValueError, match="read-only"):
        record.positions[0, 0] = 99.0


# with_frame -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("template", "dump", "atol"),
    [
        # Same step, %.16g against write_data's %.16e: float64 round-off.
        pytest.param("ring.bond.data", "ring.bond.lammpstrj", 1e-12,
                     id="unscaled"),
        # dump atom's %g: 5e-7 of a scaled value, times the box length 6.
        pytest.param("ring.bond.data", "format.atom.lammpstrj", 3.0e-6,
                     id="scaled"),
        # Triclinic x adds the errors of lx = 6, xy = 1, and xz = 0.5.
        pytest.param("triclinic.data", "triclinic.atom.lammpstrj", 3.75e-6,
                     id="triclinic-scaled"),
    ],
)
def test_with_frame_takes_positions_by_atom_id(
        data_dir, section_rows, template, dump, atol):
    """The frame LAMMPS dumped at write_data's step gives its positions."""
    base = read_data(data_dir / template)
    frame = read_dump_frame(data_dir / dump, 2, sort_by_id=False)
    expected = section_rows(data_dir / template, "Atoms", 8)

    moved = with_frame(base, frame)

    np.testing.assert_allclose(moved.positions, expected[:, 3:6],
                               rtol=0, atol=atol)
    np.testing.assert_array_equal(moved.image_flags, expected[:, -3:])


def test_with_frame_takes_the_frame_box(make_lammps_data, make_lammps_dump):
    """The box bounds come from the frame, not from the template."""
    base = read_data(make_lammps_data())
    frame = read_dump_frame(
        make_lammps_dump(box=((-1.0, 9.0), (0.0, 12.0), (2.0, 7.0))), 0)

    moved = with_frame(base, frame)

    np.testing.assert_array_equal(moved.box.lo, [-1.0, 0.0, 2.0])
    np.testing.assert_array_equal(moved.box.hi, [9.0, 12.0, 7.0])


def test_with_frame_takes_the_frame_tilt(data_dir):
    """A triclinic frame gives its tilt factors to an orthogonal template."""
    base = read_data(data_dir / "ring.bond.data")
    frame = read_dump_frame(data_dir / "triclinic.atom.lammpstrj", 0)

    moved = with_frame(base, frame)

    # in.triclinic: region prism 0 6 0 6 0 6 1.0 0.5 -0.5.
    np.testing.assert_array_equal(moved.box.tilt, [1.0, 0.5, -0.5])


def test_with_frame_matches_a_template_listed_out_of_id_order(make_record):
    """Atoms of a template not sorted by ID still get their own positions."""
    record = make_record("bond", ids=np.array([3, 1, 2]),
                         bonds=np.array([[1, 3, 1], [1, 1, 2]]))
    box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
    frame = DumpFrame(timestep=0, box=box, columns=("id", "x", "y", "z"),
                      data=np.array([[1, 1.0, 0.0, 0.0], [2, 2.0, 0.0, 0.0],
                                     [3, 3.0, 0.0, 0.0]]))

    moved = with_frame(record, frame)

    np.testing.assert_array_equal(moved.positions[:, 0], [3.0, 1.0, 2.0])


def test_with_frame_keeps_the_template_topology(data_dir):
    """Masses, bonds, molecule IDs, and types come from the template."""
    base = read_data(data_dir / "ring.angle.data")
    frame = read_dump_frame(data_dir / "ring.angle.lammpstrj", 0)

    moved = with_frame(base, frame)

    np.testing.assert_array_equal(moved.bonds, base.bonds)
    np.testing.assert_array_equal(moved.angles, base.angles)
    np.testing.assert_array_equal(moved.masses, base.masses)
    np.testing.assert_array_equal(moved.mol_ids, base.mol_ids)


def test_with_frame_takes_frame_velocities_and_drops_the_template_ones(
        make_record):
    """Velocities come from the frame when it has them, else none."""
    record = make_record("bond", velocities=np.ones((3, 3)))
    box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
    frame = DumpFrame(
        timestep=0, box=box, columns=("id", "x", "y", "z", "vx", "vy", "vz"),
        data=np.array([[3, 1.0, 1.0, 1.0, 0.3, 0.0, 0.0],
                       [1, 2.0, 2.0, 2.0, 0.1, 0.0, 0.0],
                       [2, 3.0, 3.0, 3.0, 0.2, 0.0, 0.0]]))

    moved = with_frame(record, frame)

    np.testing.assert_array_equal(moved.velocities[:, 0], [0.1, 0.2, 0.3])
    assert moved.image_flags is None


@pytest.mark.parametrize(
    ("columns", "rows", "match"),
    [
        pytest.param(("x", "y", "z"), [[1.0, 1.0, 1.0]] * 3,
                     "no 'id' column", id="no-ids"),
        pytest.param(("id", "x", "y", "z"),
                     [[1, 1.0, 1.0, 1.0], [2, 2.0, 2.0, 2.0]],
                     "atom IDs differ", id="fewer-atoms"),
        pytest.param(("id", "type", "x", "y", "z"),
                     [[1, 2, 1.0, 1.0, 1.0], [2, 1, 2.0, 2.0, 2.0],
                      [3, 2, 3.0, 3.0, 3.0]],
                     "atom types differ", id="other-types"),
    ],
)
def test_with_frame_rejects_a_frame_of_another_system(
        make_record, columns, rows, match):
    """A frame whose atoms do not match the template's raises ValueError."""
    box = Box(lo=np.zeros(3), hi=np.full(3, 10.0))
    frame = DumpFrame(timestep=0, box=box, columns=columns,
                      data=np.array(rows))

    with pytest.raises(ValueError, match=match):
        with_frame(make_record("bond"), frame)
