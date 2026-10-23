"""Tests for the whole-level builders, ensemble reducers, and space
builders in :mod:`polyphys.manage.organizer`.

The whole-level tests write their input under a directory whose name
contains a hyphen, as every lineage directory name does. Expected whole
names come from the ``SumRuleCyl`` parser applied to the bare filename, so
each test checks that the builder's result does not depend on the directory
the file sits in.

The space tests build a small ``analysis`` phase on disk with ``tmp_path``.
Expected attribute values come from the ``SumRuleCyl`` parser applied to the
same filenames, and the expected ``t_index`` and ``time`` values come from the
``adump`` and ``dt`` fields encoded in those filenames.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from polyphys.manage.organizer import (
    ens_avg_from_bin_edge,
    ens_avg_from_df,
    ens_from_bin_edge,
    ensemble,
    ensemble_avg,
    space_hists,
    space_tseries,
    whole_from_dist_mat_t,
    whole_from_file,
    whole_from_segments,
)
from polyphys.manage.parser import SumRuleCyl

GROUP = "all"
DIVISOR = 0.025
# Two state points that differ only in the number of crowders, so their bulk
# crowder volume fractions differ. Both dump every 5000 steps with dt=0.005.
ENSEMBLES = [
    "N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000",
    "N200epsilon5.0r10.5lz25.0sig2.0nc2000dt0.005bdump2000adump5000",
]
ADUMP = 5000
DT = 0.005
N_ROWS = 3


def _ens_avg_dir(tmp_path: Path) -> Path:
    """Create and return the ``ensAvg`` directory of a synthetic space."""
    ens_avg_dir = tmp_path / "analysis" / f"N200-{GROUP}-ensAvg"
    ens_avg_dir.mkdir(parents=True)
    return ens_avg_dir


def _write_tseries_space(tmp_path: Path, prop: str) -> Path:
    """Write two ensemble-averaged time series and one decoy file."""
    ens_avg_dir = _ens_avg_dir(tmp_path)
    for offset, ens in enumerate(ENSEMBLES):
        mean = np.arange(N_ROWS, dtype=float) + 10.0 * offset
        pd.DataFrame({
            f"{prop}-mean": mean,
            f"{prop}-var": np.full(N_ROWS, 0.1),
            f"{prop}-sem": np.full(N_ROWS, 0.05),
        }).to_csv(ens_avg_dir / f"{ens}-{GROUP}-{prop}-ensAvg.csv",
                  index=False)
    # A file of another property in the same directory must be ignored.
    pd.DataFrame({"other-mean": [1.0]}).to_csv(
        ens_avg_dir / f"{ENSEMBLES[0]}-{GROUP}-other-ensAvg.csv", index=False
    )
    return ens_avg_dir


def _write_hist_space(tmp_path: Path, prop: str) -> Path:
    """Write two ensemble-averaged histograms and one decoy file."""
    ens_avg_dir = _ens_avg_dir(tmp_path)
    for offset, ens in enumerate(ENSEMBLES):
        pd.DataFrame({
            f"{prop}-mean": np.array([1.0, 2.0, 1.0]) * (offset + 1),
            f"{prop}-var": np.full(N_ROWS, 0.1),
            f"{prop}-sem": np.full(N_ROWS, 0.05),
        }).to_csv(ens_avg_dir / f"{ens}-{GROUP}-{prop}-ensAvg.csv",
                  index=False)
    pd.DataFrame({"other-mean": [1.0]}).to_csv(
        ens_avg_dir / f"{ENSEMBLES[0]}-{GROUP}-other-ensAvg.csv", index=False
    )
    return ens_avg_dir


def _parsers() -> list[SumRuleCyl]:
    """Return a parser instance for each synthetic ensemble, in order."""
    return [SumRuleCyl(ens, "ensemble_long", GROUP) for ens in ENSEMBLES]


def _assert_crowder_columns(space: pd.DataFrame, attr: str) -> None:
    """Check the crowder column and its rounding against the parsers."""
    expected = np.repeat([getattr(p, attr) for p in _parsers()], N_ROWS)
    np.testing.assert_allclose(space[attr].to_numpy(), expected)
    rounded = space["phi_c_bulk_round"].to_numpy()
    # Rounding to the nearest multiple of DIVISOR moves a value by at most
    # half a step and lands on an integer number of steps.
    assert np.all(np.abs(rounded - expected) <= DIVISOR / 2 + 1e-12)
    np.testing.assert_allclose(
        rounded / DIVISOR, np.round(rounded / DIVISOR), atol=1e-9
    )
    # The two state points have different crowder fractions.
    assert len(np.unique(rounded)) == len(ENSEMBLES)


def test_space_tseries_concatenates_ensembles(tmp_path: Path) -> None:
    """space_tseries reads every '-<prop>-ensAvg.csv' file and no other."""
    ens_avg_dir = _write_tseries_space(tmp_path, "gyrTMon")
    space = space_tseries(
        str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
        GROUP, "adump",
    )
    assert len(space) == len(ENSEMBLES) * N_ROWS
    assert list(space.index) == list(range(len(ENSEMBLES) * N_ROWS))
    assert "other-mean" not in space.columns
    # Files are sorted by name, so the nc1000 rows come first.
    np.testing.assert_allclose(
        space["gyrTMon-mean"].to_numpy(), [0.0, 1.0, 2.0, 10.0, 11.0, 12.0]
    )
    np.testing.assert_allclose(
        space["dcyl"].to_numpy(),
        np.repeat([p.dcyl for p in _parsers()], N_ROWS),
    )


def test_space_tseries_scales_t_index_and_time(tmp_path: Path) -> None:
    """t_index is frame index times adump; time is t_index times dt."""
    ens_avg_dir = _write_tseries_space(tmp_path, "gyrTMon")
    space = space_tseries(
        str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl, "N*", [], GROUP,
        "adump",
    )
    frames = np.tile(np.arange(N_ROWS), len(ENSEMBLES))
    np.testing.assert_array_equal(space["t_index"].to_numpy(), frames * ADUMP)
    np.testing.assert_allclose(space["time"].to_numpy(), frames * ADUMP * DT)


def test_space_tseries_adds_phi_bulk_c_by_default(tmp_path: Path) -> None:
    """The crowder fraction comes from phi_bulk_c without being listed."""
    ens_avg_dir = _write_tseries_space(tmp_path, "gyrTMon")
    space = space_tseries(
        str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
        GROUP, "adump", divisor=DIVISOR,
    )
    assert "phi_c_bulk" not in space.columns
    _assert_crowder_columns(space, "phi_bulk_c")


def test_space_tseries_reads_named_crowder_attr(tmp_path: Path) -> None:
    """crowder_attr selects the parser attribute that is rounded."""
    ens_avg_dir = _write_tseries_space(tmp_path, "gyrTMon")
    space = space_tseries(
        str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl, "N*", [], GROUP,
        "adump", divisor=DIVISOR, crowder_attr="phi_bulk_m",
    )
    assert "phi_bulk_c" not in space.columns
    np.testing.assert_allclose(
        space["phi_bulk_m"].to_numpy(),
        np.repeat([p.phi_bulk_m for p in _parsers()], N_ROWS),
    )


def test_space_tseries_saves_space_csv(tmp_path: Path) -> None:
    """is_save writes the returned frame to the analysis space directory."""
    ens_avg_dir = _write_tseries_space(tmp_path, "gyrTMon")
    space = space_tseries(
        str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
        GROUP, "adump", is_save=True,
    )
    saved = (tmp_path / "analysis" / f"N200-{GROUP}-space"
             / f"N200-{GROUP}-gyrTMon-space.csv")
    assert saved.is_file()
    pd.testing.assert_frame_equal(pd.read_csv(saved), space)


def test_space_hists_concatenates_ensembles(tmp_path: Path) -> None:
    """space_hists reads every '-<prop>-ensAvg.csv' file and no other."""
    ens_avg_dir = _write_hist_space(tmp_path, "rhoMon")
    bin_center = np.array([0.5, 1.5, 2.5])
    space = space_hists(
        str(ens_avg_dir) + "/", "rhoMon", SumRuleCyl, "N*", ["dcyl"], GROUP,
        bin_center=bin_center, normalize=True, divisor=DIVISOR,
    )
    assert len(space) == len(ENSEMBLES) * N_ROWS
    assert "other-mean" not in space.columns
    np.testing.assert_allclose(
        space["bin_center"].to_numpy(), np.tile(bin_center, len(ENSEMBLES))
    )
    np.testing.assert_allclose(
        space["bin_center-norm"].to_numpy(),
        np.tile(bin_center / bin_center.max(), len(ENSEMBLES)),
    )
    # Each ensemble's normalized histogram sums to one.
    norm = space["rhoMon-norm"].to_numpy().reshape(len(ENSEMBLES), N_ROWS)
    np.testing.assert_allclose(norm.sum(axis=1), 1.0)
    np.testing.assert_allclose(
        space["dcyl"].to_numpy(),
        np.repeat([p.dcyl for p in _parsers()], N_ROWS),
    )
    assert "phi_c_bulk" not in space.columns
    _assert_crowder_columns(space, "phi_bulk_c")


def test_space_hists_saves_space_csv(tmp_path: Path) -> None:
    """is_save writes the returned frame to the analysis space directory."""
    ens_avg_dir = _write_hist_space(tmp_path, "rhoMon")
    space = space_hists(
        str(ens_avg_dir) + "/", "rhoMon", SumRuleCyl, "N*", [], GROUP,
        is_save=True,
    )
    saved = (tmp_path / "analysis" / f"N200-{GROUP}-space"
             / f"N200-{GROUP}-rhoMon-space.csv")
    assert saved.is_file()
    pd.testing.assert_frame_equal(pd.read_csv(saved), space)


def test_ens_from_bin_edge_returns_shared_edges() -> None:
    """Identical whole edges pass through as a copy of the shared array."""
    edges = np.array([0.0, 1.0, 2.0])
    name, out = ens_from_bin_edge(("ens1", {"w1": edges, "w2": edges.copy()}))
    assert name == "ens1"
    np.testing.assert_array_equal(out, edges)
    assert not np.shares_memory(out, edges)


@pytest.mark.parametrize(
    "other",
    [np.array([0.0, 1.0, 3.0]), np.array([0.0, 1.0, 2.0, 3.0])],
    ids=["different-values", "different-lengths"],
)
def test_ens_from_bin_edge_raises_on_mismatch(other: np.ndarray) -> None:
    """Differing whole edges raise ValueError instead of being merged."""
    with pytest.raises(ValueError, match="Bin edges are not identical"):
        ens_from_bin_edge(
            ("ens1", {"w1": np.array([0.0, 1.0, 2.0]), "w2": other})
        )


def test_ens_avg_from_bin_edge_copies_1d_edges() -> None:
    """A single 1D edge array, as ensemble_avg passes it, is copied."""
    edges = np.array([0.0, 0.5, 1.0])
    out = ens_avg_from_bin_edge(edges)
    np.testing.assert_array_equal(out, edges)
    assert not np.shares_memory(out, edges)


def test_ens_avg_from_bin_edge_checks_sequence() -> None:
    """A sequence of equal edges reduces to one; a mismatch raises."""
    edges = np.array([0.0, 0.5, 1.0])
    np.testing.assert_array_equal(
        ens_avg_from_bin_edge([edges, edges.copy()]), edges
    )
    with pytest.raises(ValueError, match="Bin edges are not identical"):
        ens_avg_from_bin_edge([edges, np.array([0.0, 0.5, 2.0])])


def test_ensemble_avg_bin_edge_keeps_edges() -> None:
    """ensemble_avg with 'bin_edge' returns each ensemble's edges intact."""
    edges = np.array([0.0, 0.5, 1.0, 1.5])
    out = ensemble_avg("rhoMon", {"ens1": edges}, GROUP, "bin_edge")
    np.testing.assert_array_equal(out["ens1"], edges)


def test_ens_avg_from_df_leaves_input_unchanged() -> None:
    """ens_avg_from_df returns a new frame and leaves its input unchanged."""
    ens_data = pd.DataFrame({
        "bin_center": [1.0, 2.0],
        "w1": [4.0, 5.0],
        "w2": [5.0, 7.0],
    })
    before = ens_data.copy()
    out = ens_avg_from_df("density", ens_data, ["bin_center"])
    assert out is not ens_data
    pd.testing.assert_frame_equal(ens_data, before)
    assert list(out.columns) == [
        "bin_center", "density-mean", "density-var", "density-sem"
    ]


def test_ens_avg_from_df_statistics() -> None:
    """Two wholes a, b give mean (a+b)/2, var (a-b)^2/2, sem |a-b|/2."""
    a = np.array([4.0, 5.0])
    b = np.array([5.0, 7.0])
    ens_data = pd.DataFrame({"bin_center": [1.0, 2.0], "w1": a, "w2": b})
    out = ens_avg_from_df("density", ens_data, ["bin_center"])
    np.testing.assert_allclose(out["density-mean"], (a + b) / 2)
    np.testing.assert_allclose(out["density-var"], (a - b) ** 2 / 2)
    np.testing.assert_allclose(out["density-sem"], np.abs(a - b) / 2)
    np.testing.assert_array_equal(out["bin_center"], [1.0, 2.0])


@pytest.fixture
def hist_wholes() -> dict[str, np.ndarray]:
    """Return two replicate histogram wholes of one state point."""
    base = (
        "N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000adump"
        "10000ens"
    )
    return {
        f"{base}1": np.array([1.0, 2.0, 3.0]),
        f"{base}2": np.array([3.0, 4.0, 5.0]),
    }


@pytest.fixture
def hist_edges(hist_wholes: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Return the same non-uniform bin edges for every whole in hist_wholes.

    The edges are non-uniform so a wrong centre formula cannot pass by
    accident.
    """
    return {name: np.array([0.0, 1.0, 3.0, 6.0]) for name in hist_wholes}


def test_ensemble_vector_bin_center_is_edge_midpoints(
    hist_wholes: dict[str, np.ndarray],
    hist_edges: dict[str, np.ndarray],
) -> None:
    """With whole_edges, a vector ensemble's bin_center is the midpoints."""
    ens = ensemble(
        "rhoMon", hist_wholes, SumRuleCyl, "bug", "vector",
        whole_edges=hist_edges,
    )
    (ens_df,) = ens.values()
    # Analytic midpoints of the edges [0, 1, 3, 6]
    np.testing.assert_allclose(ens_df["bin_center"], [0.5, 2.0, 4.5])
    for name, values in hist_wholes.items():
        np.testing.assert_array_equal(ens_df[name], values)


def test_ensemble_without_whole_edges_has_no_bin_center(
    hist_wholes: dict[str, np.ndarray],
) -> None:
    """Without whole_edges, a vector ensemble has one column per whole only."""
    ens = ensemble("rhoMon", hist_wholes, SumRuleCyl, "bug", "vector")
    (ens_df,) = ens.values()
    assert list(ens_df.columns) == list(hist_wholes)


def test_ensemble_avg_keeps_bin_center_from_ensemble(
    hist_wholes: dict[str, np.ndarray],
    hist_edges: dict[str, np.ndarray],
) -> None:
    """ensemble_avg passes an ensemble's bin_center through unchanged."""
    ens = ensemble(
        "rhoMon", hist_wholes, SumRuleCyl, "bug", "vector",
        whole_edges=hist_edges,
    )
    (ens_df,) = ens.values()
    avg = ensemble_avg("rhoMon", ens, "bug", "dataframe")
    (avg_df,) = avg.values()
    np.testing.assert_array_equal(avg_df["bin_center"], ens_df["bin_center"])
    assert list(avg_df.columns) == [
        "bin_center", "rhoMon-mean", "rhoMon-var", "rhoMon-sem"
    ]


@pytest.mark.parametrize("whole_type", ["matrix", "dataframe", "bin_edge"])
def test_ensemble_rejects_whole_edges_for_non_vector(
    whole_type: str,
    hist_wholes: dict[str, np.ndarray],
    hist_edges: dict[str, np.ndarray],
) -> None:
    """whole_edges with any whole_type other than 'vector' raises."""
    with pytest.raises(ValueError, match="'whole_edges' is only valid"):
        ensemble(
            "rhoMon", hist_wholes, SumRuleCyl, "bug", whole_type,
            whole_edges=hist_edges,
        )


# Whole-level builders. These use their own particle group.
WHOLE_GROUP = "bug"
WHOLE = (
    "N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000adump"
    "10000ens2"
)


def _hyphenated_dir(tmp_path: Path) -> Path:
    """Create and return a directory named like a lineage directory."""
    directory = tmp_path / "N2000-bug-wholeSim"
    directory.mkdir()
    return directory


def _expected_whole(filename: str, lineage: str) -> str:
    """Return the whole name the parser gives for a bare filename."""
    return getattr(SumRuleCyl(filename, lineage, WHOLE_GROUP), "whole")


def test_whole_from_segments_ignores_hyphenated_directory(tmp_path):
    """Segments under a hyphenated directory join into their parsed whole."""
    directory = _hyphenated_dir(tmp_path)
    paths = []
    for seg_id in (1, 2):
        path = directory / f"{WHOLE}.j0{seg_id}.{WHOLE_GROUP}.gyrTMon.npy"
        np.save(path, np.full(2, float(seg_id)))
        paths.append(str(path))

    wholes = whole_from_segments(
        "gyrTMon", paths, SumRuleCyl, WHOLE_GROUP, "tseries"
    )

    expected = _expected_whole(Path(paths[0]).stem, "segment")
    assert list(wholes) == [expected]
    np.testing.assert_array_equal(wholes[expected], [1.0, 1.0, 2.0, 2.0])


def test_whole_from_file_ignores_hyphenated_directory(tmp_path):
    """A whole under a hyphenated directory loads under its parsed name."""
    path = _hyphenated_dir(tmp_path) / f"{WHOLE}.{WHOLE_GROUP}.gyrTMon.npy"
    np.save(path, np.arange(3.0))

    wholes = whole_from_file([str(path)], SumRuleCyl, WHOLE_GROUP)

    expected = _expected_whole(path.stem, "whole")
    assert list(wholes) == [expected]
    np.testing.assert_array_equal(wholes[expected], np.arange(3.0))


def test_whole_from_dist_mat_t_ignores_hyphenated_directory(tmp_path):
    """A distance matrix under a hyphenated directory keys its parsed name."""
    filename = f"{WHOLE}.{WHOLE_GROUP}.distMatTFoci.npy"
    path = _hyphenated_dir(tmp_path) / filename
    np.save(path, np.zeros((2, 2)))
    frame = pd.DataFrame({"x": [0.0]})

    def func(whole_path, whole_info):
        return frame, frame, frame

    freqs, rdfs, tseries = whole_from_dist_mat_t(
        [(str(path),)], SumRuleCyl, WHOLE_GROUP, func
    )

    expected = _expected_whole(path.stem, "whole")
    assert list(freqs) == list(rdfs) == list(tseries) == [expected]
