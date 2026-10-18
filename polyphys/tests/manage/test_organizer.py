"""Tests for the space builders and ensemble reducers in
:mod:`polyphys.manage.organizer`.

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
    ensemble_avg,
    space_hists,
    space_tseries,
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
    """space_tseries finds every '-<prop>-ensAvg.csv' file (issue #44)."""
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
    """space_hists finds every '-<prop>-ensAvg.csv' file (issue #44)."""
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
    """Differing whole edges raise ValueError, not merge (issue #45)."""
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
    """ens_avg_from_df returns a new frame and keeps the wholes (issue #45)."""
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
