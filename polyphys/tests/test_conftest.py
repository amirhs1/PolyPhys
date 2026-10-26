"""Tests for the shared fixtures and the ``private_data`` marker in
``polyphys/tests/conftest.py``.

Expected file counts and indices come from the arguments each test passes.
Expected file contents come from the formula documented in
``lineage_tree_value``. Expected LAMMPS section layouts come from the LAMMPS
manual's read_data and dump pages: ``atom-ID atom-type x y z`` for
``atomic``, ``atom-ID molecule-ID atom-type x y z`` for ``bond`` and
``angle``, and one bond and one angle per consecutive pair and triple of a
linear chain. The marker tests run a nested pytest session with
``pytester``.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from polyphys.manage import parser as parsers


def _sections(text: str) -> dict[str, list[str]]:
    """Split a LAMMPS data file into its header and its named sections."""
    blocks = [b.strip().splitlines() for b in text.split("\n\n")]
    sections = {"header": blocks[1]}
    for title, body in zip(blocks[2:], blocks[3:]):
        if len(title) == 1 and not title[0][0].isdigit():
            sections[title[0].split("#")[0].strip()] = body
    return sections


def _header_counts(header: list[str]) -> dict[str, int]:
    """Map each header keyword, such as 'bonds', to its count."""
    return {line.split(maxsplit=1)[1]: int(line.split()[0])
            for line in header}


def _dump_frames(text: str) -> list[list[str]]:
    """Split a LAMMPS dump into frames, each a list of its lines."""
    frames = text.split("ITEM: TIMESTEP\n")[1:]
    return [frame.strip().splitlines() for frame in frames]


@pytest.fixture(
    params=[
        "TwoMonDepCub", "SumRuleCyl", "SumRuleCubHeteroRing",
        "SumRuleCubHeteroLinear", "TransFociCyl", "TransFociCub", "HnsCub",
        "HnsCyl",
    ],
)
def parser_class(request: pytest.FixtureRequest) -> type:
    """Each parser class in polyphys.manage.parser, one per test run."""
    return getattr(parsers, request.param)


@pytest.fixture
def private_data_project(pytester: pytest.Pytester) -> pytest.Pytester:
    """A nested pytest project with this conftest and one private_data test."""
    pytester.makeconftest(Path(__file__).with_name("conftest.py").read_text())
    pytester.makepyfile(
        test_private="""
        import pytest

        @pytest.mark.private_data
        def test_reads_data(private_data_dir):
            assert private_data_dir.is_dir()
        """
    )
    return pytester


def test_make_lineage_tree_writes_one_file_per_segment(make_lineage_tree):
    """A segment tree holds n_ensembles * n_wholes * n_segments files."""
    tree = make_lineage_tree(
        parsers.SumRuleCyl, n_ensembles=2, n_wholes=3, n_segments=2
    )

    assert len(tree.paths) == 12
    assert all(path.is_file() for path in tree.paths)


def test_make_lineage_tree_names_parse_without_missing_keywords(
    make_lineage_tree, capsys, parser_class
):
    """Every segment filename yields every keyword of the parser's grammar."""
    tree = make_lineage_tree(parser_class)

    parsed = [parser_class(path, "segment", tree.group, ispath=True)
              for path in tree.paths]

    assert len(parsed) == len(tree.paths)
    # The parsers report a missing keyword on stdout rather than raising.
    assert capsys.readouterr().out == ""


def test_make_lineage_tree_groups_into_wholes_and_ensembles(
    make_lineage_tree, parser_class
):
    """Segments group into one whole per replicate, one ensemble per state."""
    tree = make_lineage_tree(
        parser_class, n_ensembles=2, n_wholes=3, n_segments=2
    )

    parsed = [parser_class(path, "segment", tree.group, ispath=True)
              for path in tree.paths]

    assert len({p.whole for p in parsed}) == 6
    assert len({p.ensemble_long for p in parsed}) == 2


def test_make_lineage_tree_index_is_read_only(make_lineage_tree):
    """A tree's file index cannot be changed after the tree is written."""
    tree = make_lineage_tree(parsers.SumRuleCyl, n_ensembles=1)

    with pytest.raises(TypeError, match="does not support item assignment"):
        tree.index[tree.paths[0]] = (9, 9, 9)


def test_make_lineage_tree_numbers_wholes_and_segments_from_one(
    make_lineage_tree,
):
    """The ens and j values in the names count wholes and segments from 1."""
    tree = make_lineage_tree(
        parsers.SumRuleCyl, n_ensembles=1, n_wholes=3, n_segments=2
    )

    parsed = [parsers.SumRuleCyl(path, "segment", tree.group, ispath=True)
              for path in tree.paths]

    assert [(p.ensemble_id, p.segment_id) for p in parsed] == [
        (1, 1), (1, 2), (2, 1), (2, 2), (3, 1), (3, 2)
    ]


@pytest.mark.parametrize(
    ("lineage", "ext", "load", "key", "expected"),
    [
        pytest.param("segment", ".npy", np.load, (2, 3, 1),
                     [2310.0, 2311.0, 2312.0], id="segment-npy"),
        pytest.param("segment", ".csv",
                     lambda p: pd.read_csv(p)["gyrTMon"].to_numpy(),
                     (2, 3, 1), [2310.0, 2311.0, 2312.0], id="segment-csv"),
        pytest.param("whole", ".npy", np.load, (2, 3, 0),
                     [2300.0, 2301.0, 2302.0], id="whole-npy"),
    ],
)
def test_make_lineage_tree_file_holds_documented_values(
    make_lineage_tree, lineage, ext, load, key, expected
):
    """A file holds 1000*ensemble + 100*whole + 10*segment + row."""
    tree = make_lineage_tree(
        parsers.SumRuleCyl, n_ensembles=2, n_wholes=3, lineage=lineage,
        ext=ext, prop="gyrTMon", n_rows=3,
    )
    path = next(p for p, index in tree.index.items() if index == key)

    values = load(path)

    # Integers below 2**53 are exact in float64 and in the CSV text.
    np.testing.assert_array_equal(values, expected)


@pytest.mark.parametrize(
    ("atom_style", "expected"),
    [
        pytest.param("atomic", {"atoms": 4, "bonds": 0, "angles": 0},
                     id="atomic"),
        pytest.param("bond", {"atoms": 4, "bonds": 3, "angles": 0},
                     id="bond"),
        pytest.param("angle", {"atoms": 4, "bonds": 3, "angles": 2},
                     id="angle"),
    ],
)
def test_make_lammps_data_header_counts_match_sections(
    make_lammps_data, atom_style, expected
):
    """A 4-atom chain's header counts equal its section line counts."""
    path = make_lammps_data(atom_style=atom_style)

    sections = _sections(path.read_text())

    counts = _header_counts(sections["header"])
    assert {key: counts.get(key, 0) for key in expected} == expected
    assert {
        "atoms": len(sections["Atoms"]),
        "bonds": len(sections.get("Bonds", [])),
        "angles": len(sections.get("Angles", [])),
    } == expected


@pytest.mark.parametrize(
    ("atom_style", "n_columns"),
    [
        pytest.param("atomic", 5, id="atomic"),
        pytest.param("bond", 6, id="bond"),
        pytest.param("angle", 6, id="angle"),
    ],
)
def test_make_lammps_data_atoms_columns_follow_style(
    make_lammps_data, atom_style, n_columns
):
    """Atoms lines have the column count LAMMPS reads for the atom style."""
    path = make_lammps_data(atom_style=atom_style)

    atoms = _sections(path.read_text())["Atoms"]

    assert {len(line.split()) for line in atoms} == {n_columns}


def test_make_lammps_dump_round_trips_positions(make_lammps_dump):
    """Positions read back from the dump match those written."""
    positions = np.linspace(0.0, 9.5, 36).reshape(3, 4, 3)
    path = make_lammps_dump(positions=positions, timesteps=[0, 500, 1000])

    frames = _dump_frames(path.read_text())

    assert [frame[0] for frame in frames] == ["0", "500", "1000"]
    assert [int(frame[2]) for frame in frames] == [4, 4, 4]
    assert {frame[7] for frame in frames} == {"ITEM: ATOMS id type x y z"}
    read = np.array([[line.split()[2:] for line in frame[8:]]
                     for frame in frames], dtype=float)
    # Six decimals are written, so a value can be off by half of 1e-6.
    np.testing.assert_allclose(read, positions, rtol=0, atol=5e-7)


@pytest.mark.parametrize(
    "boundary",
    [pytest.param("pp pp pp", id="periodic"),
     pytest.param("ff ff pp", id="cylinder")],
)
def test_make_lammps_dump_writes_boundary_flags(make_lammps_dump, boundary):
    """Every BOX BOUNDS line carries the requested boundary flags."""
    path = make_lammps_dump(boundary=boundary)

    frames = _dump_frames(path.read_text())

    assert {frame[3] for frame in frames} == {
        f"ITEM: BOX BOUNDS {boundary}"
    }


def test_private_data_skips_when_data_dir_unset(
    private_data_project, monkeypatch
):
    """A private_data test is skipped when POLYPHYS_DATA_DIR is unset."""
    monkeypatch.delenv("POLYPHYS_DATA_DIR", raising=False)

    result = private_data_project.runpytest_subprocess()

    result.assert_outcomes(skipped=1)


def test_private_data_runs_when_data_dir_exists(
    private_data_project, monkeypatch, tmp_path
):
    """A private_data test runs when POLYPHYS_DATA_DIR names a directory."""
    monkeypatch.setenv("POLYPHYS_DATA_DIR", str(tmp_path))

    result = private_data_project.runpytest_subprocess()

    result.assert_outcomes(passed=1)


def test_private_data_errors_when_data_dir_missing(
    private_data_project, monkeypatch, tmp_path
):
    """A private_data test errors when POLYPHYS_DATA_DIR is no directory."""
    monkeypatch.setenv("POLYPHYS_DATA_DIR", str(tmp_path / "missing"))

    result = private_data_project.runpytest_subprocess()

    result.assert_outcomes(errors=1)


def test_strict_markers_rejects_misspelt_private_data(pytester):
    """A misspelt private_data marker fails collection."""
    pytester.makeconftest(Path(__file__).with_name("conftest.py").read_text())
    pytester.makepyfile(
        test_typo="""
        import pytest

        @pytest.mark.private_dta
        def test_reads_data():
            pass
        """
    )

    result = pytester.runpytest_subprocess("--strict-markers")

    result.assert_outcomes(errors=1)
