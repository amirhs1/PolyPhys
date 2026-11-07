"""Tests for :mod:`polyphys.analyze.pipeline`.

Most inputs are lineage trees from the ``make_lineage_tree`` factory in
``conftest.py``, written into one probe directory per test, one property per
subfolder. A tree's values are ``1000 e + 100 w + 10 s + k``
(``lineage_tree_value``), so two wholes of one ensemble differ by exactly
100 in every row and the mean, variance, and standard error across them are
``w1 + 50``, ``5000``, and ``50`` by hand. Bin edges, distance matrices, and
two-dimensional wholes are written over the tree's files by hand. The hooks
are also tested on their own with a ``HookContext``; the density kernel is a
stub that records its calls. The equilibrium values of the wholes are read
from ensemble CSV files written by hand. Two tests reduce the output of a
probe run in the test: ``ChainSize`` on chains written by
``make_chain_files``, which skip without MDAnalysis and take their expected
names from the parser and their expected means from the probe's own
results. The golden comparison with the research code's output is written
blind and runs only with the maintainer's data.
"""
import dataclasses
import importlib
import logging
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.clusters import foci_pair_info
from polyphys.analyze.correlations import acf
from polyphys.analyze.stats import sem_blocked
from polyphys.analyze.pipeline import (
    Autocorrelation,
    Derived,
    FociPairDistances,
    HookContext,
    LocalDensities,
    TimeAverage,
    equilibrium_wholes,
    half_contour,
    reduce_space,
    reduce_stamps,
    space_hists_from_spec,
    space_tseries_from_spec,
)
from polyphys.manage import organizer
from polyphys.tests.synthetic import ChainCyl, RingCub
from polyphys.spec import (
    AnalysisTable,
    ProjectSpec,
    edge,
    histogram,
    histogram_2d,
    matrix_tseries,
    tseries,
)
from polyphys.tests.conftest import lineage_tree_value, write_lineage_tree

# A Pandas4Warning marks a deprecated pandas call; fail on it.
pytestmark = pytest.mark.filterwarnings("error::pandas.errors.Pandas4Warning")

# The glob that lists every file of a probe directory built by write_tree:
# <probe>/<property>/<Parser>-<group>-<stage>/<file>.
HIERARCHY = "*/*/*"


# Fixtures: specs, a probe directory, and tree writers


@pytest.fixture
def cyl_spec():
    """A ChainCyl spec whose crowder size is ``dcrowd``."""
    return ProjectSpec(
        parser=ChainCyl, group="bug", divisor=0.025,
        crowder_attr="phi_bulk_c",
        physical_attrs=("space",),
        species_sizes={"Mon": "dmon", "Crd": "dcrowd"},
        tseries_properties=(),
    )


@pytest.fixture
def foci_spec():
    """A RingCub spec: a ring of 20 small and 2 large monomers."""
    return ProjectSpec(
        parser=RingCub, group="bug", divisor=0.025,
        crowder_attr="phi_bulk_c",
        physical_attrs=("space",),
        species_sizes={"Mon": "dmon_small", "Foci": "dmon_large"},
        tseries_properties=(),
    )


@pytest.fixture
def probe_dir(tmp_path):
    """A probe directory, ``<root>/probes/N-bug-probe``, next to which
    ``reduce_space`` writes ``<root>/analysis/``."""
    probe = tmp_path / "probes" / "N-bug-probe"
    probe.mkdir(parents=True)
    return probe


@pytest.fixture
def write_tree(probe_dir):
    """Factory: write a lineage tree of one property into the probe
    directory and return it."""
    def _write(parser, prop, **kwargs):
        return write_lineage_tree(probe_dir / prop, parser, prop=prop,
                                  **kwargs)
    return _write


@pytest.fixture
def write_edges():
    """Factory: write one edge file per file of a histogram tree, with the
    same edges for every file unless ``per_whole`` maps a whole index to
    its own edges."""
    def _write(tree, edge_prop, edges, per_whole=None):
        for path, (_, whole, _) in tree.index.items():
            chosen = edges if per_whole is None else per_whole[whole]
            np.save(path.with_name(path.name.replace(tree.prop, edge_prop)),
                    chosen)
    return _write


class RecordingKernel:
    """A density kernel that records its calls and returns twice and three
    times the counts."""

    def __init__(self):
        self.calls = []

    def __call__(self, counts, edges, info, direction, diameter):
        self.calls.append((info.whole, counts, edges, direction, diameter))
        return 2.0 * counts, 3.0 * counts


def _parsed(tree, parser, attr):
    """Map each tree file's (ensemble, whole) indices to a parsed name."""
    return {
        index[:2]: getattr(parser(str(path), tree.lineage, tree.group,
                                  ispath=True), attr)
        for path, index in tree.index.items()
    }


def _analysis(probe_dir, stage):
    """The analysis directory of a stage next to the probe directory."""
    group = probe_dir.name.split("-")[1]
    return probe_dir.parent.parent / "analysis" / f"N-{group}-{stage}"


# Records


def test_derived_rejects_edges_on_non_vector():
    """Edges belong to a vector histogram only."""
    with pytest.raises(ValueError, match="edges are only valid"):
        Derived("xyHistCrd", {"w": np.zeros((2, 2))}, "matrix",
                edges={"w": np.arange(3.0)})


def test_derived_rejects_unknown_whole_type():
    """A derived artifact needs one of the four whole types."""
    with pytest.raises(ValueError, match="whole_type must be"):
        Derived("gyrMon", {"w": np.zeros(2)}, "frame")


def test_dataframe_output_must_hold_dataframes(probe_dir, cyl_spec,
                                               write_tree):
    """A hook that declares ``'dataframe'`` wholes but returns arrays is
    an error, not a silent reinterpretation."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)

    def hook(wholes, context):
        return (Derived("gyrMon", wholes, "dataframe"),)

    table = AnalysisTable("bug", "whole", False,
                          (tseries("gyrT", "Mon", post_hook=hook),))
    with pytest.raises(ValueError, match="must be a pandas.DataFrame"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


# Vector properties through the lineage


def test_two_wholes_give_exact_mean_var_sem(probe_dir, cyl_spec,
                                            write_tree):
    """Two wholes 100 apart in every row average to w1 + 50 with variance
    5000 and standard error 50, by the ddof = 1 formulas."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1,
               n_wholes=2, n_rows=3)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          save=False)

    (average,) = result["gyrTMon"].values()
    assert_allclose(average["gyrTMon-mean"],
                    lineage_tree_value(1, 1, 0, 3) + 50)
    assert_allclose(average["gyrTMon-var"], 5000.0)
    assert_allclose(average["gyrTMon-sem"], 50.0)


def test_ensemble_csv_has_one_column_per_whole_named_by_parser(
        probe_dir, cyl_spec, write_tree):
    """The ensemble file holds one column per whole, under the whole names
    the parser reads from the written files."""
    tree = write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=3)
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    ensemble = pd.read_csv(
        _analysis(probe_dir, "ens") / f"{ens_name}-bug-gyrTMon.csv")
    assert list(ensemble.columns) == [names[(1, 1)], names[(1, 2)]]
    assert_allclose(ensemble[names[(1, 2)]], lineage_tree_value(1, 2, 0, 3))


def test_segment_tseries_concatenate_in_segment_order(probe_dir, cyl_spec,
                                                      write_tree):
    """A segment-lineage time series joins its segments j01 then j02 into
    the whole file."""
    tree = write_tree(ChainCyl, "gyrTMon", n_ensembles=1, n_wholes=1,
                      n_segments=2, n_rows=2)
    (whole_name,) = set(_parsed(tree, ChainCyl, "whole").values())
    table = AnalysisTable("bug", "segment", False,
                          (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    whole = np.load(
        _analysis(probe_dir, "wholeSim") / f"{whole_name}-bug-gyrTMon.npy")
    assert_array_equal(whole, np.concatenate([lineage_tree_value(1, 1, 1, 2),
                                              lineage_tree_value(1, 1, 2, 2)]))


def test_segments_are_ordered_naturally(probe_dir, cyl_spec):
    """Segments j9 and j10 join in numeric order, not in text order, which
    puts j10 first."""
    stem = "N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000ens1"
    folder = probe_dir / "gyrTMon" / "ChainCyl-bug-segment"
    folder.mkdir(parents=True)
    for seg_id in (9, 10):
        np.save(folder / f"{stem}.j{seg_id}.bug.gyrTMon.npy",
                np.array([float(seg_id)]))
    table = AnalysisTable("bug", "segment", False,
                          (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    whole = np.load(
        _analysis(probe_dir, "wholeSim") / f"{stem}-bug-gyrTMon.npy")
    assert_array_equal(whole, [9.0, 10.0])


def test_segment_histograms_sum_and_ensembles_carry_bin_centers(
        probe_dir, cyl_spec, write_tree, write_edges):
    """A segment-lineage histogram sums its segments element-wise, and its
    ensemble gets the midpoints of the shared edges as ``bin_center``."""
    tree = write_tree(ChainCyl, "rHistMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=3)
    write_edges(tree, "rEdgeMon", np.array([0.0, 1.0, 2.0, 4.0]))
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    table = AnalysisTable("bug", "segment", False, (histogram("r", "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    whole = np.load(
        _analysis(probe_dir, "wholeSim") / f"{names[(1, 1)]}-bug-rHistMon.npy")
    assert_array_equal(whole, lineage_tree_value(1, 1, 1, 3)
                       + lineage_tree_value(1, 1, 2, 3))
    ensemble = pd.read_csv(
        _analysis(probe_dir, "ens") / f"{ens_name}-bug-rHistMon.csv")
    assert_array_equal(ensemble["bin_center"], [0.5, 1.5, 3.0])
    average = result["rHistMon"][ens_name]
    assert list(average.columns) == ["bin_center", "rHistMon-mean",
                                     "rHistMon-var", "rHistMon-sem"]
    assert_allclose(average["rHistMon-mean"], whole + 100)


def test_mismatched_edges_between_wholes_of_one_ensemble_raise(
        probe_dir, cyl_spec, write_tree, write_edges):
    """Two wholes of one ensemble with different edges raise, naming the
    ensemble and both wholes."""
    tree = write_tree(ChainCyl, "rHistMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=3)
    write_edges(tree, "rEdgeMon", None, per_whole={
        1: np.array([0.0, 1.0, 2.0, 3.0]), 2: np.array([0.0, 1.0, 2.0, 3.5])})
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    table = AnalysisTable("bug", "whole", False, (histogram("r", "Mon"),))

    with pytest.raises(ValueError, match="different bin edges") as excinfo:
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)

    message = str(excinfo.value)
    assert ens_name in message
    assert names[(1, 1)] in message and names[(1, 2)] in message


def test_missing_edge_file_for_a_whole_raises(probe_dir, cyl_spec,
                                              write_tree, write_edges):
    """A histogram whose edge files do not cover every whole raises."""
    tree = write_tree(ChainCyl, "rHistMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=3)
    write_edges(tree, "rEdgeMon", np.arange(4.0))
    (first_edge,) = sorted(
        (probe_dir / "rHistMon").rglob("*ens1-bug-rEdgeMon.npy"))
    first_edge.unlink()
    table = AnalysisTable("bug", "whole", False, (histogram("r", "Mon"),))

    with pytest.raises(ValueError, match="no 'rEdgeMon' file"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


def test_property_without_files_raises(probe_dir, cyl_spec, write_tree):
    """A property that matches no file raises instead of writing nothing
    (the research code's ``Crdtwo`` row was silent)."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False,
                          (histogram("r", "Crdtwo"),))

    with pytest.raises(FileNotFoundError, match="rHistCrdtwo"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


def test_missing_property_raises_before_any_directory_is_made(
        probe_dir, cyl_spec, write_tree):
    """A property that matches no file, such as a typo, raises before the
    analysis directories are made or another property is reduced."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False,
                          (tseries("gyrT", "Mon"), tseries("gyrT", "Mno")))

    with pytest.raises(FileNotFoundError, match="gyrTMno"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    assert not (probe_dir.parent.parent / "analysis").exists()


def test_repeated_output_of_an_undeclared_hook_raises(probe_dir, cyl_spec,
                                                      write_tree):
    """Two outputs of one name from hooks that declare no outputs raise at
    the second instead of overwriting the first."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    write_tree(ChainCyl, "fsdTMon", lineage="whole", n_ensembles=1)

    def hook(wholes, context):
        return (Derived("sizeMon", wholes, "vector"),)

    table = AnalysisTable("bug", "whole", False, (
        tseries("gyrT", "Mon", post_hook=hook),
        tseries("fsdT", "Mon", post_hook=hook),
    ))

    with pytest.raises(ValueError,
                       match="output 'sizeMon' is named more than once"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


def test_default_hierarchy_is_the_spec_dir_file_layout(probe_dir, cyl_spec,
                                                       write_tree):
    """Without ``hierarchy`` the files are listed under
    ``<spec.hierarchy>/<spec.hierarchy>``, the research code's layout."""
    tree = write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    # Move the tree's folder under a name the ensemble_long prefix matches.
    tree.paths[0].parent.rename(probe_dir / "N20ac1.0")
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, save=False)

    assert cyl_spec.hierarchy == "N*"
    assert len(result["gyrTMon"]) == 1


def test_empty_glob_raises(probe_dir, cyl_spec):
    """A hierarchy that matches no file raises."""
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))
    with pytest.raises(FileNotFoundError, match="No file matches"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


@pytest.mark.parametrize("prop, n_rows", [
    pytest.param("gyrTMon", 3, id="gyrTMon"),
    pytest.param("xgyrTMon", 5, id="xgyrTMon"),
])
def test_file_selection_matches_the_last_name_token_exactly(
        probe_dir, cyl_spec, write_tree, prop, n_rows):
    """With ``gyrTMon`` and ``xgyrTMon`` files side by side, each property
    reads only its own files, told apart by their row counts."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1,
               n_rows=3)
    write_tree(ChainCyl, "xgyrTMon", lineage="whole", n_ensembles=1,
               n_rows=5)
    name = prop[:-len("Mon")]
    table = AnalysisTable("bug", "whole", False, (tseries(name, "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          save=False)

    (average,) = result[prop].values()
    assert len(average) == n_rows


def test_outputs_land_in_the_analysis_stage_folders(probe_dir, cyl_spec,
                                                    write_tree):
    """A segment-lineage run writes the wholeSim, ens, and ensAvg folders
    next to the probe phase, with the organizer's file names."""
    tree = write_tree(ChainCyl, "gyrTMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=2)
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    table = AnalysisTable("bug", "segment", False,
                          (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    analysis = probe_dir.parent.parent / "analysis"
    assert sorted(path.name for path in analysis.iterdir()) == [
        "N-bug-ens", "N-bug-ensAvg", "N-bug-wholeSim"]
    assert sorted(path.name for path in (analysis / "N-bug-wholeSim")
                  .iterdir()) == sorted(
        f"{name}-bug-gyrTMon.npy" for name in names.values())
    assert [path.name for path in (analysis / "N-bug-ens").iterdir()] == [
        f"{ens_name}-bug-gyrTMon.csv"]
    assert [path.name for path in (analysis / "N-bug-ensAvg").iterdir()] == [
        f"{ens_name}-bug-gyrTMon-ensAvg.csv"]


def test_whole_lineage_writes_no_wholesim_folder(probe_dir, cyl_spec,
                                                 write_tree):
    """A whole-lineage run writes the ens and ensAvg folders only."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    analysis = probe_dir.parent.parent / "analysis"
    assert sorted(path.name for path in analysis.iterdir()) == [
        "N-bug-ens", "N-bug-ensAvg"]


def test_save_false_writes_nothing(probe_dir, cyl_spec, write_tree):
    """With ``save=False`` no analysis folder appears."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY, save=False)

    assert not (probe_dir.parent.parent / "analysis").exists()


def test_reduce_space_logs_each_property(probe_dir, cyl_spec, write_tree,
                                         caplog):
    """Progress is logged at INFO, one record per property."""
    write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    with caplog.at_level(logging.INFO, logger="polyphys.analyze.pipeline"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)

    assert any("'gyrTMon'" in record.message for record in caplog.records)


# Matrix and edge properties


def test_matrix_outputs_average_element_wise(probe_dir, cyl_spec,
                                             write_tree):
    """A two-dimensional histogram's ensemble average is taken element by
    element across the wholes, as one array per statistic."""
    tree = write_tree(ChainCyl, "xyHistCrd", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=1)
    matrices = {1: np.array([[1.0, 2.0], [3.0, 4.0]]),
                2: np.array([[3.0, 2.0], [7.0, 0.0]])}
    for path, (_, whole, _) in tree.index.items():
        np.save(path, matrices[whole])
    table = AnalysisTable("all", "whole", False, (histogram_2d("xy", "Crd"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          save=False)

    (average,) = result["xyHistCrd"].values()
    assert sorted(average) == ["xyHistCrd-mean", "xyHistCrd-sem",
                               "xyHistCrd-var"]
    assert_allclose(average["xyHistCrd-mean"], [[2.0, 2.0], [5.0, 2.0]])
    assert_allclose(average["xyHistCrd-var"], [[2.0, 0.0], [8.0, 8.0]])
    assert_allclose(average["xyHistCrd-sem"], [[1.0, 0.0], [2.0, 2.0]])


def test_edge_property_passes_shared_edges_through(probe_dir, cyl_spec,
                                                   write_tree):
    """An edge property's ensemble average is the edges themselves."""
    tree = write_tree(ChainCyl, "xEdge", group="all", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=1)
    edges = np.array([-1.0, 0.0, 1.0])
    for path in tree.paths:
        np.save(path, edges)
    table = AnalysisTable("all", "whole", False, (edge("x"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          save=False)

    (average,) = result["xEdge"].values()
    assert_array_equal(average, edges)


# Hooks on their own


def test_time_average_drops_first_column_and_names_the_output():
    """``TimeAverage`` averages over the frames, drops column 0 on request,
    and names the output after the property's species."""
    hook = TimeAverage("clustersHist", drop_first_column=True)
    context = HookContext(matrix_tseries("clustersHistT", "Foci"))
    whole = np.array([[0, 1, 2, 3], [0, 3, 4, 5], [0, 5, 6, 7]])

    (out,) = hook({"w": whole}, context)

    assert (out.name, out.whole_type) == ("clustersHistFoci", "vector")
    assert_allclose(out.wholes["w"], [3.0, 4.0, 5.0])


def test_time_average_keeps_every_column_by_default():
    """Without the flag, column 0 is averaged like the others."""
    (out,) = TimeAverage("bondsHist")(
        {"w": np.array([[1.0, 2.0], [3.0, 4.0]])},
        HookContext(matrix_tseries("bondsHistT", "Foci")))
    assert_allclose(out.wholes["w"], [2.0, 3.0])


def test_autocorrelation_passes_series_through_and_adds_acf(probe_dir,
                                                            cyl_spec,
                                                            write_tree):
    """``Autocorrelation`` writes the series under its own name and the
    autocorrelation of each whole, equal to ``correlations.acf``, under
    ``<property>-acf``."""
    tree = write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=1)
    rng = np.random.default_rng(7)
    series = {}
    for path, (_, whole, _) in tree.index.items():
        values = np.empty(64)
        values[0] = rng.standard_normal()
        for t in range(1, 64):
            values[t] = 0.6 * values[t - 1] + rng.standard_normal()
        np.save(path, values)
        series[whole] = values
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    table = AnalysisTable("bug", "whole", False, (
        tseries("gyrT", "Mon", post_hook=Autocorrelation(nlags=5)),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    assert set(result) == {"gyrTMon", "gyrTMon-acf"}
    ensemble = pd.read_csv(
        _analysis(probe_dir, "ens") / f"{ens_name}-bug-gyrTMon-acf.csv")
    for whole, values in series.items():
        assert_allclose(ensemble[names[(1, whole)]],
                        acf(values, 5, alpha=None).values)
    assert (_analysis(probe_dir, "ensAvg")
            / f"{ens_name}-bug-gyrTMon-acf-ensAvg.csv").exists()


def test_local_densities_calls_the_kernel_once_per_whole(probe_dir, cyl_spec,
                                                         write_tree,
                                                         write_edges):
    """``LocalDensities`` calls the kernel with each whole's counts, edges,
    direction, and the species' diameter from the spec, and reduces the
    histogram and both profiles on the same bins."""
    tree = write_tree(ChainCyl, "rHistCrd", group="all", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=3)
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    write_edges(tree, "rEdgeCrd", edges)
    names = _parsed(tree, ChainCyl, "whole")
    kernel = RecordingKernel()
    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=LocalDensities()),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          density=kernel, save=False)

    assert set(result) == {"rHistCrd", "rRhoCrd", "rPhiCrd"}
    assert sorted(call[0] for call in kernel.calls) == sorted(names.values())
    for whole_name, counts, call_edges, direction, diameter in kernel.calls:
        index = next(idx for idx, name in names.items() if name == whole_name)
        assert_array_equal(counts, lineage_tree_value(*index, 0, 3))
        assert_array_equal(call_edges, edges)
        assert direction == "r"
        assert diameter == 1.0   # ac1.0 in the written names
    (rho,) = result["rRhoCrd"].values()
    (phi,) = result["rPhiCrd"].values()
    assert_allclose(rho["rRhoCrd-mean"],
                    2.0 * (lineage_tree_value(1, 1, 0, 3) + 50))
    assert_allclose(phi["rPhiCrd-mean"],
                    3.0 * (lineage_tree_value(1, 1, 0, 3) + 50))
    assert_array_equal(phi["bin_center"], [0.5, 1.5, 2.5])


def test_local_densities_without_a_kernel_raises(probe_dir, cyl_spec,
                                                 write_tree, write_edges):
    """A ``LocalDensities`` property needs the kernel passed to
    ``reduce_space``."""
    tree = write_tree(ChainCyl, "rHistCrd", group="all", lineage="whole",
                      n_ensembles=1)
    write_edges(tree, "rEdgeCrd", np.arange(4.0))
    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=LocalDensities()),))

    with pytest.raises(ValueError, match="density kernel"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


def test_local_densities_needs_a_species_with_a_size(cyl_spec):
    """A species the spec gives no size for is an error, not a guess."""
    context = HookContext(histogram("r", "Dna"), spec=cyl_spec, group="all",
                          edges={"w": np.arange(3.0)},
                          density=RecordingKernel())
    with pytest.raises(ValueError, match="'Dna' has no size"):
        LocalDensities()({"w": np.ones(2)}, context)


def test_local_densities_needs_the_spec_and_the_edges(cyl_spec):
    """Without the study's spec or the histogram's edges the hook cannot
    find the species' size or the bins."""
    context = HookContext(histogram("r", "Crd"), spec=cyl_spec, group="all",
                          density=RecordingKernel())
    with pytest.raises(ValueError, match="needs the study's spec and"):
        LocalDensities()({"w": np.ones(2)}, context)


def test_local_densities_rejects_a_non_histogram_property(cyl_spec):
    """``LocalDensities`` is for ``histogram(...)`` properties only."""
    context = HookContext(tseries("gyrT", "Mon"), spec=cyl_spec,
                          density=RecordingKernel())
    with pytest.raises(ValueError, match="needs a histogram"):
        LocalDensities()({"w": np.ones(2)}, context)


# The error of the per-whole normalized profile


def test_derived_rejects_a_normalizer_on_non_vector():
    """A per-whole normalizer belongs to a vector profile only."""
    with pytest.raises(ValueError, match="normalizer is only valid"):
        Derived("xyHistCrd", {"w": np.zeros((2, 2))}, "matrix",
                normalizer="first")


def test_derived_rejects_an_unknown_normalizer():
    """A per-whole normalizer is 'first', 'last', or 'max'."""
    with pytest.raises(ValueError, match="normalizer must be"):
        Derived("rRhoCrd", {"w": np.zeros(2)}, "vector", normalizer="mean")


@pytest.mark.parametrize("direction, expected", [
    pytest.param("r", "first", id="cylinder-r"),
    pytest.param("z", "max", id="cylinder-z"),
])
def test_local_densities_take_the_sum_rule_default_normalizer(
    cyl_spec, direction, expected
):
    """The two density profiles carry the ``PROFILE_NORMALIZERS`` entry of
    the spec's geometry and the histogram's direction; the histogram does
    not."""
    whole = ("N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000"
             "ens1")
    context = HookContext(
        histogram(direction, "Crd"), spec=cyl_spec, group="all",
        infos={whole: ChainCyl(whole, "whole", "all")},
        edges={whole: np.arange(3.0)}, density=RecordingKernel())
    hist, rho, phi = LocalDensities()({whole: np.ones(2)}, context)
    assert hist.normalizer is None
    assert (rho.normalizer, phi.normalizer) == (expected, expected)


@pytest.fixture
def write_density_wholes(probe_dir):
    """Factory: write the histogram and edge files of two wholes of one
    ChainCyl state point, group 'all', in the probe directory, and
    return the histogram's property name."""
    stem = "N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000"

    def _write(direction, counts_1, counts_2, edges=(0.0, 1.0, 2.0, 3.0)):
        for k, counts in enumerate((counts_1, counts_2), start=1):
            np.save(probe_dir / f"{stem}ens{k}-all-{direction}HistCrd.npy",
                    np.asarray(counts, dtype=float))
            np.save(probe_dir / f"{stem}ens{k}-all-{direction}EdgeCrd.npy",
                    np.asarray(edges))
        return f"{direction}HistCrd"

    return _write


@pytest.mark.parametrize("direction, counts_1, counts_2, expected", [
    # 'first': the wholes 2 x [1, 2, 4] and 2 x [2, 2, 2], each over its
    # own bin 0, are [1, 2, 4] and [1, 1, 1]; the sem of two values with
    # ddof = 1 is half their difference.
    pytest.param("r", [1.0, 2.0, 4.0], [2.0, 2.0, 2.0], [0.0, 0.5, 1.5],
                 id="cylinder-r-first"),
    # 'max': the ensemble mean 2 x [1.5, 2, 3] peaks in bin 2, so the
    # wholes are [1/4, 1/2, 1] and [1, 1, 1].
    pytest.param("z", [1.0, 2.0, 4.0], [2.0, 2.0, 2.0], [0.375, 0.25, 0.0],
                 id="cylinder-z-max"),
])
def test_density_profile_error_is_the_spread_over_wholes_of_the_ratio(
    probe_dir, cyl_spec, write_density_wholes, direction, counts_1,
    counts_2, expected
):
    """``<prop>-norm-sem`` is the standard error over the wholes of each
    whole's profile divided by its own value in the bin the ensemble's
    normalizer names; the ensemble average itself is not normalized."""
    write_density_wholes(direction, counts_1, counts_2)
    table = AnalysisTable("all", "whole", False, (
        histogram(direction, "Crd", post_hook=LocalDensities()),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*",
                          density=RecordingKernel(), save=False)

    (rho,) = result[f"{direction}RhoCrd"].values()
    # Exact: ratios and differences of small dyadic floats
    assert_allclose(rho[f"{direction}RhoCrd-norm-sem"], expected, rtol=1e-15,
                    atol=0)
    assert_allclose(rho[f"{direction}RhoCrd-mean"],
                    np.add(counts_1, counts_2))
    assert f"{direction}HistCrd-norm-sem" not in next(
        iter(result[f"{direction}HistCrd"].values())).columns


def test_density_profile_error_is_scale_free(probe_dir, cyl_spec,
                                             write_density_wholes):
    """Property: the number density and the volume fraction, which differ
    by a constant factor, carry the same ``-norm-sem``."""
    write_density_wholes("r", [1.0, 2.0, 4.0], [2.0, 2.0, 2.0])
    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=LocalDensities()),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*",
                          density=RecordingKernel(), save=False)

    (rho,) = result["rRhoCrd"].values()
    (phi,) = result["rPhiCrd"].values()
    # float64 round-off of one division
    assert_allclose(phi["rPhiCrd-norm-sem"], rho["rRhoCrd-norm-sem"],
                    rtol=1e-14, atol=0)


def test_density_profile_error_is_nan_and_warns_on_a_zero_reference(
    probe_dir, cyl_spec, write_density_wholes
):
    """A whole that is zero in the reference bin gives NaN in every bin of
    ``-norm-sem``, with a warning naming the ensemble and the bin."""
    write_density_wholes("r", [0.0, 2.0, 4.0], [2.0, 2.0, 2.0])
    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=LocalDensities()),))

    with pytest.warns(UserWarning, match="zero in its reference bin 0"):
        result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*",
                              density=RecordingKernel(), save=False)

    (rho,) = result["rRhoCrd"].values()
    assert rho["rRhoCrd-norm-sem"].isna().all()


def test_density_profile_error_is_saved_in_the_ensavg_csv(
    probe_dir, cyl_spec, write_density_wholes
):
    """The ``-norm-sem`` column is in the ensAvg file the pipeline writes,
    after the organizer's three statistics."""
    write_density_wholes("r", [1.0, 2.0, 4.0], [2.0, 2.0, 2.0])
    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=LocalDensities()),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy="*",
                 density=RecordingKernel())

    (path,) = (probe_dir.parent.parent / "analysis" / "N-all-ensAvg").glob(
        "*-all-rRhoCrd-ensAvg.csv")
    saved = pd.read_csv(path)
    assert list(saved.columns) == ["bin_center", "rRhoCrd-mean",
                                   "rRhoCrd-var", "rRhoCrd-sem",
                                   "rRhoCrd-norm-sem"]
    assert_allclose(saved["rRhoCrd-norm-sem"], [0.0, 0.5, 1.5], rtol=1e-15,
                    atol=0)


def test_density_profile_error_takes_the_last_bin_as_reference(
    probe_dir, cyl_spec, write_density_wholes
):
    """A ``'last'`` normalizer divides each whole by its own value in the
    final bin, not in the first bin or the bin where the mean peaks."""
    write_density_wholes("r", [4.0, 2.0, 1.0], [2.0, 2.0, 2.0])

    def last_bin(wholes, context):
        return (Derived(context.prop.prop, wholes, "vector", context.edges,
                        "last"),)

    table = AnalysisTable("all", "whole", False, (
        histogram("r", "Crd", post_hook=last_bin),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*",
                          save=False)

    (hist,) = result["rHistCrd"].values()
    # The wholes over their own bin 2 are [4, 2, 1] and [1, 1, 1]; the sem
    # of two values with ddof = 1 is half their difference. Bin 0, which
    # 'first' and 'max' would pick, gives [0, 0.25, 0.375] instead. Exact:
    # ratios and differences of small dyadic floats.
    assert_allclose(hist["rHistCrd-norm-sem"], [1.5, 0.5, 0.0], rtol=1e-15,
                    atol=0)


# The foci pair distances


@pytest.fixture
def foci_context(foci_spec):
    """A context for ``distMatTFoci`` of the RingCub ring."""
    return HookContext(matrix_tseries("distMatT", "Foci"), spec=foci_spec)


def _distance_frames(genomic_pos, separations):
    """Stack symmetric distance matrices with `genomic_pos` on the diagonal
    of frame 0; `separations` lists the upper triangle of each frame."""
    n_foci = len(genomic_pos)
    first, second = np.triu_indices(n_foci, 1)
    frames = []
    for values in separations:
        frame = np.zeros((n_foci, n_foci))
        frame[first, second] = values
        frame[second, first] = values
        frames.append(frame)
    frames[0][np.diag_indices(n_foci)] = genomic_pos
    return np.stack(frames)


def test_foci_pair_distances_counts_and_tags(foci_context):
    """Three foci on a ring of 10 over two frames: the tags are the pairs
    of ``foci_pair_info`` with genomic distances in bonds, the counts are
    the hand tallies, and the separations are the frames' values."""
    whole = _distance_frames([0, 4, 8], [[1.5, 2.5, 0.5], [1.5, 2.9, 2.5]])
    hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)

    hist, density, series = hook({"w": whole}, foci_context)

    pairs = foci_pair_info([0, 4, 8], 10, "ring")
    tags = [f"loc{i}loc{j}genDist{d}" for i, j, d in pairs]
    assert tags == ["loc0loc4genDist4", "loc0loc8genDist2",
                    "loc4loc8genDist4"]
    assert (hist.name, density.name, series.name) == (
        "pairDistHistFoci", "pairDistRdfFoci", "pairDistTFoci")
    frame = hist.wholes["w"]
    assert list(frame.columns) == tags + ["bin_center"]
    assert frame["loc0loc4genDist4"].to_list() == [0, 2, 0]
    assert frame["loc0loc8genDist2"].to_list() == [0, 0, 2]
    assert frame["loc4loc8genDist4"].to_list() == [1, 0, 1]
    assert_array_equal(frame["bin_center"], [0.5, 1.5, 2.5])
    assert series.wholes["w"]["loc4loc8genDist4"].to_list() == [0.5, 2.5]


def test_foci_pair_density_integrates_to_one_per_pair(foci_context):
    """Each pair's density satisfies ``sum(P_k dV_k) = 1`` over the shells."""
    whole = _distance_frames([0, 4, 8], [[1.5, 2.5, 0.5], [1.5, 2.9, 2.5]])
    hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)

    _, density, _ = hook({"w": whole}, foci_context)

    frame = density.wholes["w"]
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    shells = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    for tag in frame.columns.drop("bin_center"):
        assert_allclose(np.sum(frame[tag] * shells), 1.0)


def test_foci_pair_histograms_conserve_the_frames(foci_context):
    """Every pair has one separation per frame, so each histogram column
    sums to ``n_frames`` when the grid covers every separation."""
    rng = np.random.default_rng(3)
    n_frames = 40
    whole = _distance_frames(
        [0, 4, 8], rng.uniform(0.0, 7.5, size=(n_frames, 3)))
    hook = FociPairDistances(bin_size=0.5, r_max=half_contour(10, 1.5),
                             n_monomers=10)

    hist, _, _ = hook({"w": whole}, foci_context)

    frame = hist.wholes["w"]
    assert_array_equal(frame.drop(columns="bin_center").sum(axis=0),
                       n_frames)


def test_foci_pair_distances_reject_a_separation_above_r_max(foci_context):
    """A separation above the grid raises and names the whole and the
    value, since a dropped frame would break the normalization."""
    whole = _distance_frames([0, 4, 8], [[1.5, 2.5, 0.5], [1.5, 3.5, 2.5]])
    hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)

    with pytest.raises(ValueError, match="whole 'w1'.*3.5.*above r_max"):
        hook({"w1": whole}, foci_context)


def test_foci_pair_distances_accept_a_separation_at_r_max(foci_context):
    """A separation exactly at an `r_max` that is a whole number of bins
    is on the grid, although ``ceil(r_max / bin_size) * bin_size`` falls
    an ulp below it: 15.9 is 530 bins of 0.03 in decimal arithmetic."""
    whole = _distance_frames([0, 4, 8], [[15.9, 1.0, 1.0]])
    hook = FociPairDistances(bin_size=0.03, r_max=15.9, n_monomers=10)

    hist, _, _ = hook({"w": whole}, foci_context)

    frame = hist.wholes["w"]
    assert len(frame) == 530
    assert frame["loc0loc4genDist4"].iloc[-1] == 1


def test_foci_pair_distances_need_the_spec_for_the_topology():
    """Without the study's spec the hook cannot know the ring topology."""
    hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)
    with pytest.raises(ValueError, match="needs the study's spec"):
        hook({"w": _distance_frames([0, 4, 8], [[1.0, 1.0, 1.0]])},
             HookContext(matrix_tseries("distMatT", "Foci")))


def test_foci_pair_distances_reject_a_non_matrix_whole(foci_context):
    """A whole that is not a ``(n_frames, n_foci, n_foci)`` array raises."""
    hook = FociPairDistances(bin_size=1.0, r_max=3.0, n_monomers=10)
    with pytest.raises(ValueError, match="distance matrix time series"):
        hook({"w": np.ones((2, 3))}, foci_context)


def test_foci_pair_distances_resolve_callables_on_the_whole(foci_context):
    """A number-valued function of the whole's parser instance is read for
    each whole."""
    foci_context = HookContext(
        foci_context.prop, spec=foci_context.spec,
        infos={"w": RingCub(
            "al2.0nl2ml8.0ns20ac1.0nc100l10.0dt0.005bdump1000adump5000ens1",
            "whole", "bug")})
    whole = _distance_frames([0, 11], [[2.5], [2.5]])
    hook = FociPairDistances(
        bin_size=lambda info: 0.5 * info.dmon_small,
        r_max=lambda info: half_contour(
            info.nmon_small + info.nmon_large, 1.5 * info.dmon_small),
        n_monomers=lambda info: info.nmon_small + info.nmon_large)

    hist, _, _ = hook({"w": whole}, foci_context)

    frame = hist.wholes["w"]
    assert list(frame.columns) == ["loc0loc11genDist11", "bin_center"]
    # bin_size 0.5 to r_max 11 * 1.5 = 16.5: 33 bins; both frames in bin 5.
    assert len(frame) == 33
    assert frame["loc0loc11genDist11"].to_list() == [0] * 5 + [2] + [0] * 27


def test_foci_pair_distances_need_a_parser_instance_for_a_callable(
        foci_context):
    """A callable attribute with no parser instance in the context raises."""
    hook = FociPairDistances(bin_size=lambda info: 1.0, r_max=3.0,
                             n_monomers=10)
    with pytest.raises(ValueError, match="no parser instance"):
        hook({"w": _distance_frames([0, 4, 8], [[1.0, 1.0, 1.0]])},
             foci_context)


@pytest.mark.parametrize("overrides, what", [
    pytest.param({"bin_size": lambda info: -1.0}, "bin_size",
                 id="bin_size-negative"),
    pytest.param({"bin_size": lambda info: 0.0}, "bin_size",
                 id="bin_size-zero"),
    pytest.param({"r_max": lambda info: float("nan")}, "r_max",
                 id="r_max-nan"),
    pytest.param({"n_monomers": lambda info: "10"}, "n_monomers",
                 id="n_monomers-string"),
])
def test_foci_pair_distances_reject_a_bad_callable_result(
        foci_context, overrides, what):
    """A callable quantity that does not give a finite positive number
    raises, naming the quantity and the whole."""
    foci_context = HookContext(foci_context.prop, spec=foci_context.spec,
                               infos={"w": object()})
    hook = FociPairDistances(
        **({"bin_size": 1.0, "r_max": 3.0, "n_monomers": 10} | overrides))
    with pytest.raises(ValueError,
                       match=f"{what} of whole 'w' must resolve to a finite"):
        hook({"w": _distance_frames([0, 4, 8], [[1.0, 1.0, 1.0]])},
             foci_context)


@pytest.mark.parametrize("n_monomers, bond_max, match", [
    pytest.param(1, 1.5, "at least 2", id="one-monomer"),
    pytest.param(10, 0.0, "must be positive", id="zero-bond"),
])
def test_half_contour_rejects_bad_input(n_monomers, bond_max, match):
    """``half_contour`` needs at least two monomers and a positive bond."""
    with pytest.raises(ValueError, match=match):
        half_contour(n_monomers, bond_max)


def test_half_contour_is_half_the_bonds_times_the_bond():
    """A ring of 11 monomers has at most 5 bonds between two monomers."""
    assert half_contour(11, 2.0) == 10.0


def test_dataframe_wholes_of_unequal_length_raise(probe_dir, cyl_spec,
                                                  write_tree):
    """Two 'dataframe' wholes of one ensemble with 3 and 2 rows raise,
    naming the ensemble and the lengths, instead of averaging over the
    first whole's rows."""
    tree = write_tree(ChainCyl, "gyrTMon", lineage="whole", n_ensembles=1,
                      n_wholes=2)
    names = _parsed(tree, ChainCyl, "whole")
    (ens_name,) = set(_parsed(tree, ChainCyl, "ensemble_long").values())
    lengths = {names[(1, 1)]: 3, names[(1, 2)]: 2}

    def hook(wholes, context):
        frames = {name: pd.DataFrame({"x": np.arange(float(lengths[name]))})
                  for name in wholes}
        return (Derived("gyrMon", frames, "dataframe"),)

    table = AnalysisTable("bug", "whole", False,
                          (tseries("gyrT", "Mon", post_hook=hook),))
    with pytest.raises(ValueError, match=re.escape(
            f"The wholes of ensemble '{ens_name}' have unequal lengths: "
            f"{names[(1, 1)]} 3, {names[(1, 2)]} 2.")):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)


def test_dataframe_wholes_reduce_to_long_ensembles_and_row_statistics(
        probe_dir, foci_spec, write_tree):
    """Through ``reduce_space``, the pair histograms of two wholes give an
    ensemble file with every whole's rows and an ensemble average with the
    mean, variance, and standard error of each pair's counts, row by row."""
    tree = write_tree(RingCub, "distMatTFoci", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=1)
    separations = {1: [[1.5], [1.5]], 2: [[2.5], [2.5]]}
    for path, (_, whole, _) in tree.index.items():
        np.save(path, _distance_frames([0, 11], separations[whole]))
    names = _parsed(tree, RingCub, "whole")
    (ens_name,) = set(_parsed(tree, RingCub, "ensemble_long").values())
    table = AnalysisTable("bug", "whole", False, (
        matrix_tseries("distMatT", "Foci", post_hook=FociPairDistances(
            bin_size=1.0, r_max=3.0,
            n_monomers=lambda info: info.nmon_small + info.nmon_large)),))

    result = reduce_space(probe_dir, foci_spec, table, hierarchy=HIERARCHY)

    assert set(result) == {"pairDistHistFoci", "pairDistRdfFoci",
                           "pairDistTFoci"}
    ensemble = pd.read_csv(
        _analysis(probe_dir, "ens") / f"{ens_name}-bug-pairDistHistFoci.csv")
    assert list(ensemble.columns) == ["whole", "row", "loc0loc11genDist11",
                                      "bin_center"]
    assert ensemble["whole"].to_list() == [names[(1, 1)]] * 3 + [
        names[(1, 2)]] * 3
    assert ensemble["loc0loc11genDist11"].to_list() == [0, 2, 0, 0, 0, 2]
    average = result["pairDistHistFoci"][ens_name]
    assert list(average.columns) == [
        "bin_center", "loc0loc11genDist11-mean", "loc0loc11genDist11-var",
        "loc0loc11genDist11-sem"]
    assert_allclose(average["loc0loc11genDist11-mean"], [0.0, 1.0, 1.0])
    assert_allclose(average["loc0loc11genDist11-var"], [0.0, 2.0, 2.0])
    assert_allclose(average["loc0loc11genDist11-sem"], [0.0, 1.0, 1.0])
    saved = pd.read_csv(_analysis(probe_dir, "ensAvg")
                        / f"{ens_name}-bug-pairDistHistFoci-ensAvg.csv")
    pd.testing.assert_frame_equal(saved, average)
    series = pd.read_csv(
        _analysis(probe_dir, "ens") / f"{ens_name}-bug-pairDistTFoci.csv")
    assert series["loc0loc11genDist11"].to_list() == [1.5, 1.5, 2.5, 2.5]


# Stamps


def _write_stamps(probe_dir, parser_name, group, whole_names, n_segments):
    """Write one segment stamp file per segment of each whole, with the
    columns the organizer's stamp reducers read."""
    folder = probe_dir / "stamps" / f"{parser_name}-{group}-segment"
    folder.mkdir(parents=True)
    for ens_id, whole in enumerate(whole_names, start=1):
        for seg_id in range(1, n_segments + 1):
            pd.DataFrame({
                "space": ["N20ac1.0"], "group": [group],
                "lineage_name": ["segment"],
                "segment": [f"{whole}.j{seg_id:02d}"], "whole": [whole],
                "ensemble_long": [whole[:-len(f"ens{ens_id}")]],
                "ensemble_id": [ens_id], "segment_id": [seg_id],
                "n_frames": [10 * seg_id],
            }).to_csv(folder / f"{whole}.j{seg_id:02d}.{group}.stamps.csv",
                      index=False)


def test_reduce_stamps_rolls_segments_up_to_wholes_and_ensembles(probe_dir):
    """Segment stamps give one row per whole, with the segments counted and
    the frames summed, then one row per ensemble with the wholes counted."""
    stem = "N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000"
    wholes = [f"{stem}ens1", f"{stem}ens2"]
    _write_stamps(probe_dir, "ChainCyl", "bug", wholes, n_segments=2)
    files = [(str(path),) for path in sorted(probe_dir.rglob("*.stamps.csv"))]

    whole_stamps, ensemble_stamps = reduce_stamps(files, "segment")

    assert whole_stamps["whole"].to_list() == wholes
    assert whole_stamps["n_segments"].to_list() == [2, 2]
    assert whole_stamps["n_frames"].to_list() == [30, 30]
    assert ensemble_stamps["ensemble_long"].to_list() == [stem]
    assert ensemble_stamps["n_ensembles"].to_list() == [2]


def test_reduce_stamps_combines_whole_stamps_directly(probe_dir):
    """Whole-lineage stamps combine into one row per whole, with their
    segment columns dropped, then roll up to one row per ensemble."""
    stem = "N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000"
    wholes = [f"{stem}ens1", f"{stem}ens2"]
    _write_stamps(probe_dir, "ChainCyl", "bug", wholes, n_segments=1)
    files = [(str(path),) for path in sorted(probe_dir.rglob("*.stamps.csv"))]

    with pytest.warns(UserWarning, match="dropped"):
        whole_stamps, ensemble_stamps = reduce_stamps(files, "whole")

    assert whole_stamps["whole"].to_list() == wholes
    assert "segment_id" not in whole_stamps.columns
    assert ensemble_stamps["n_ensembles"].to_list() == [2]


def test_reduce_space_writes_the_stamp_files_of_each_stage(probe_dir,
                                                           cyl_spec):
    """With ``stamps=True`` the segment, whole, and ensemble stamps land in
    the wholeSim, ens, and ensAvg folders."""
    stem = "N20r5.5lz25.0ac1.0nc100dt0.005bdump1000adump5000"
    _write_stamps(probe_dir, "ChainCyl", "bug", [f"{stem}ens1"],
                  n_segments=2)
    table = AnalysisTable("bug", "segment", True, ())

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    assert (_analysis(probe_dir, "wholeSim")
            / "N20ac1.0-bug-stamps-segment.csv").exists()
    assert (_analysis(probe_dir, "ens") / "N20ac1.0-bug-stamps-whole.csv"
            ).exists()
    assert (_analysis(probe_dir, "ensAvg")
            / "N20ac1.0-bug-stamps-ensemble_long.csv").exists()


# End to end


def test_segment_tree_end_to_end(probe_dir, cyl_spec, write_tree,
                                 write_edges):
    """One segment-lineage run of a time series with its autocorrelation, a
    histogram with density profiles, and a time-averaged histogram gives
    the hand-composed averages of the segment, whole, and ensemble steps."""
    gyr = write_tree(ChainCyl, "gyrTMon", n_ensembles=1, n_wholes=2,
                     n_segments=2, n_rows=3)
    hist = write_tree(ChainCyl, "rHistMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=3)
    write_edges(hist, "rEdgeMon", np.arange(4.0))
    clusters = write_tree(ChainCyl, "clustersHistTMon", n_ensembles=1,
                          n_wholes=2, n_segments=2, n_rows=2)
    for path, (e, w, s) in clusters.index.items():
        # Two frames per segment, three columns: the tree value plus 0, 1, 2.
        np.save(path, lineage_tree_value(e, w, s, 2)[:, None]
                + np.arange(3.0))
    kernel = RecordingKernel()
    table = AnalysisTable("bug", "segment", False, (
        tseries("gyrT", "Mon", post_hook=Autocorrelation(nlags=2)),
        histogram("r", "Mon", post_hook=LocalDensities()),
        matrix_tseries("clustersHistT", "Mon", post_hook=TimeAverage(
            "clustersHist", drop_first_column=True)),
    ))
    (ens_name,) = set(_parsed(gyr, ChainCyl, "ensemble_long").values())

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                          density=kernel)

    assert set(result) == {"gyrTMon", "gyrTMon-acf", "rHistMon", "rRhoMon",
                           "rPhiMon", "clustersHistMon"}
    whole_1 = np.concatenate([lineage_tree_value(1, 1, 1, 3),
                              lineage_tree_value(1, 1, 2, 3)])
    assert_allclose(result["gyrTMon"][ens_name]["gyrTMon-mean"],
                    whole_1 + 50)
    hist_1 = lineage_tree_value(1, 1, 1, 3) + lineage_tree_value(1, 1, 2, 3)
    assert_allclose(result["rHistMon"][ens_name]["rHistMon-mean"],
                    hist_1 + 100)
    assert_allclose(result["rRhoMon"][ens_name]["rRhoMon-mean"],
                    2.0 * (hist_1 + 100))
    # The time average over the four frames of whole 1 is the mean of its
    # two segments' values, plus the column offsets 1 and 2; whole 2 is 100
    # higher, so the ensemble mean is 50 higher than whole 1.
    frames_1 = np.concatenate([lineage_tree_value(1, 1, 1, 2),
                               lineage_tree_value(1, 1, 2, 2)])
    clusters_avg = result["clustersHistMon"][ens_name]
    assert_allclose(clusters_avg["clustersHistMon-mean"],
                    frames_1.mean() + np.array([1.0, 2.0]) + 50)
    # The wholeSim folder holds the combined wholes and every derived whole.
    (names,) = [set(_parsed(hist, ChainCyl, "whole").values())]
    assert sorted(path.name for path in _analysis(probe_dir, "wholeSim")
                  .iterdir()) == sorted(
        f"{name}-bug-{prop}.npy" for name in names for prop in (
            "gyrTMon", "gyrTMon-acf", "rHistMon", "rEdgeMon", "rRhoMon",
            "rPhiMon", "clustersHistTMon", "clustersHistMon"))


# A probe's own output, reduced


@pytest.fixture
def probe_space(probe_dir, make_chain_files):
    """Run ``ChainSize`` on two segments of each of two replicates of one
    ChainCyl state point, saving its ``rfloryT`` series and stamps into
    one folder per segment of the probe directory. Each chain is straight
    along z with bond length ``e + s / 10 + f / 100`` for replicate ``e``,
    segment ``s``, and frame ``f``, so every value differs. Returns the
    probe's ``rfloryT`` results by (replicate, segment). The dumps carry no
    time step, so the tests that use it ignore the dump reader's warning
    that the frames have no ``dt``."""
    pytest.importorskip("MDAnalysis")
    from polyphys.probe import ChainSize, make_universe, write_stamps

    stem = "N4r2.5lz10ac1nc0dt0.005bdump1adump1"
    results = {}
    for ens in (1, 2):
        for seg in (1, 2):
            name = f"{stem}ens{ens}.j{seg:02d}"
            bonds = ens + seg / 10 + np.arange(3) / 100
            frames = np.full((3, 4, 3), 10.0)
            frames[:, :, 2] = 5.0 + bonds[:, None] * np.arange(4)
            data, dump = make_chain_files(frames, name=name)
            parser = ChainCyl(name, "segment", "bug")
            folder = probe_dir / name
            folder.mkdir()
            universe = make_universe(data, dump)
            analysis = ChainSize(universe.select_atoms("resid 1"),
                                 ("rfloryT",)).run()
            analysis.save(folder, f"{parser.name}-{parser.group}")
            write_stamps(parser, analysis.n_frames, folder)
            results[(ens, seg)] = analysis.results["rfloryT"]
    return results


@pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")
def test_probe_outputs_reduce_to_the_stage_files(probe_dir, cyl_spec,
                                                 probe_space):
    """The files a probe writes parse in the pipeline, which writes each
    whole of the probed segments to wholeSim and the ensemble and its
    average to ens and ensAvg, with the stamps of each level, under the
    names the parser gives the probed segments."""
    parsers = [ChainCyl(next((probe_dir / f).glob("*.npy")).name, "segment",
                        "bug", ispath=True)
               for f in sorted(path.name for path in probe_dir.iterdir())]
    (space,) = {parser.space for parser in parsers}
    (ensemble,) = {parser.ensemble_long for parser in parsers}
    wholes = sorted({parser.whole for parser in parsers})
    table = AnalysisTable("bug", "segment", True, (tseries("rfloryT", "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*/*")

    assert list(result["rfloryTMon"]) == [ensemble]
    assert sorted(path.name for path in _analysis(probe_dir, "wholeSim")
                  .iterdir()) == sorted(
        [f"{space}-bug-stamps-segment.csv"]
        + [f"{whole}-bug-rfloryTMon.npy" for whole in wholes])
    assert sorted(path.name for path in _analysis(probe_dir, "ens")
                  .iterdir()) == sorted([
                      f"{space}-bug-stamps-whole.csv",
                      f"{ensemble}-bug-rfloryTMon.csv"])
    assert sorted(path.name for path in _analysis(probe_dir, "ensAvg")
                  .iterdir()) == sorted([
                      f"{space}-bug-stamps-ensemble_long.csv",
                      f"{ensemble}-bug-rfloryTMon-ensAvg.csv"])


@pytest.mark.filterwarnings(
    "ignore:Reader has no dt information:UserWarning")
def test_probe_outputs_reduce_to_the_mean_of_the_wholes(probe_dir, cyl_spec,
                                                        probe_space):
    """The ensemble average of a probed series is, frame by frame, the mean
    over the replicates of each replicate's segments joined in order: a
    property of the reduction, computed from the probe's own results."""
    whole_1 = np.concatenate([probe_space[(1, 1)], probe_space[(1, 2)]])
    whole_2 = np.concatenate([probe_space[(2, 1)], probe_space[(2, 2)]])
    table = AnalysisTable("bug", "segment", False,
                          (tseries("rfloryT", "Mon"),))

    result = reduce_space(probe_dir, cyl_spec, table, hierarchy="*/*",
                          save=False)

    (average,) = result["rfloryTMon"].values()
    # Both sides average the same float64 copies of the probe's values, so
    # only float64 round-off separates them.
    assert_allclose(average["rfloryTMon-mean"], (whole_1 + whole_2) / 2,
                    rtol=1e-12)


# The equilibrium values of the wholes


@pytest.fixture
def ens_dir(tmp_path):
    """An ``ens`` stage directory, ``<root>/analysis/N-bug-ens``."""
    folder = tmp_path / "analysis" / "N-bug-ens"
    folder.mkdir(parents=True)
    return folder


@pytest.fixture
def write_ensemble(ens_dir):
    """Factory: write the ensemble CSV of one property, one column per
    whole named ``<stem>ens<i>``, and return the whole names."""
    def _write(prop, columns, stem=None):
        stem = stem or ("N20r5.5lz25.0ac1.0nc100"
                        "dt0.005bdump1000adump5000")
        frame = pd.DataFrame({f"{stem}ens{i}": column
                              for i, column in enumerate(columns, start=1)})
        frame.to_csv(ens_dir / f"{stem}-bug-{prop}.csv", index=False)
        return list(frame.columns)
    return _write


@pytest.fixture
def alternating():
    """A series of 256 frames alternating between 0 and 2: its time
    average is 1, its variance 256/255, and its blocked error 0, since
    every block of two frames averages to 1."""
    return np.tile([0.0, 2.0], 128)


def test_equilibrium_wholes_reads_each_ensemble_file_once(
        ens_dir, write_ensemble, alternating, monkeypatch):
    """Two properties over two ensembles read four files, each once."""
    wholes = (write_ensemble("gyrTMon", [alternating, alternating + 1])
              + write_ensemble("gyrTMon", [alternating], stem="N20ac1.0nc0"))
    write_ensemble("shapeTMon", [alternating, alternating])
    write_ensemble("shapeTMon", [alternating], stem="N20ac1.0nc0")
    stamps = pd.DataFrame({"whole": wholes})
    reads = []
    read_csv = pd.read_csv

    def counting_read_csv(path, *args, **kwargs):
        reads.append(Path(path).name)
        return read_csv(path, *args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", counting_read_csv)

    equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon", "shapeTMon": "shapeMon"},
                       stamps)

    assert len(reads) == 4
    assert len(set(reads)) == 4


@pytest.mark.parametrize("method", [
    pytest.param("lee", id="lee"),
    pytest.param("sem_err", id="sem_err"),
])
def test_equilibrium_wholes_mean_and_var_are_exact_and_sem_is_blocked(
        ens_dir, write_ensemble, method):
    """The mean and variance of each whole are the hand values of an
    alternating series and numpy's of a random one, and the error is
    ``sem_blocked`` of that whole's column under the chosen rule."""
    n_frames = 16384
    rng = np.random.default_rng(7)
    series = [np.tile([0.0, 2.0], n_frames // 2),
              rng.standard_normal(n_frames),
              np.repeat(rng.standard_normal(n_frames // 8), 8)]
    wholes = write_ensemble("gyrTMon", series)
    stamps = pd.DataFrame({"whole": wholes})

    table = equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps,
                               method=method)

    # The alternating series: deviations of exactly 1, so var = n/(n - 1);
    # the sums are exact in float64, one rounding in the division.
    assert table["gyrMon-mean"][0] == pytest.approx(1.0, rel=1e-12)
    assert table["gyrMon-var"][0] == pytest.approx(
        n_frames / (n_frames - 1), rel=1e-12)
    assert table["gyrMon-sem"][0] == 0.0
    for row in (1, 2):
        # The same numpy reductions on the same values: round-off only.
        assert table["gyrMon-mean"][row] == pytest.approx(
            series[row].mean(), rel=1e-12)
        assert table["gyrMon-var"][row] == pytest.approx(
            series[row].var(ddof=1), rel=1e-12)
        assert table["gyrMon-sem"][row] == pytest.approx(
            sem_blocked(series[row], method=method), rel=1e-12)


def test_equilibrium_wholes_error_grows_with_the_correlation_time(
        ens_dir, write_ensemble):
    """A series that holds each value for 8 frames has a blocked error
    about ``sqrt(8)`` times the per-frame one; an independent series has
    about the same."""
    rng = np.random.default_rng(7)
    independent = rng.standard_normal(16384)
    held = np.repeat(rng.standard_normal(2048), 8)
    wholes = write_ensemble("gyrTMon", [independent, held])
    stamps = pd.DataFrame({"whole": wholes})

    table = equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps)

    naive = [series.std(ddof=1) / np.sqrt(series.size)
             for series in (independent, held)]
    # The bounds are those of the sem_blocked doctest on the same signal:
    # a statistical band around sqrt(8) = 2.83 and around 1.
    assert 2.3 < table["gyrMon-sem"][1] / naive[1] < 3.4
    assert 0.8 < table["gyrMon-sem"][0] / naive[0] < 1.25


def test_equilibrium_wholes_gives_a_constant_series_zero_error(
        ens_dir, write_ensemble, alternating):
    """A whole whose series never varies has an exact average and no
    error, instead of the blocking analysis rejecting it."""
    wholes = write_ensemble("gyrTMon", [np.full(256, 3.0), alternating])
    stamps = pd.DataFrame({"whole": wholes})

    table = equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps)

    assert table["gyrMon-mean"][0] == 3.0
    assert table["gyrMon-var"][0] == 0.0
    assert table["gyrMon-sem"][0] == 0.0


def test_equilibrium_wholes_keeps_the_stamps_order_and_columns(
        ens_dir, write_ensemble, alternating):
    """The table holds the stamps' columns, then the three statistics of
    each property in the order given, with the rows in the stamps' order,
    not the files'."""
    wholes = write_ensemble("gyrTMon", [alternating, alternating + 1])
    write_ensemble("shapeTMon", [alternating + 2, alternating + 3])
    stamps = pd.DataFrame({"whole": wholes[::-1], "n_frames": [256, 256]})

    table = equilibrium_wholes(
        ens_dir, {"shapeTMon": "shapeMon", "gyrTMon": "gyrMon"}, stamps)

    assert list(table.columns) == [
        "whole", "n_frames", "shapeMon-mean", "shapeMon-var", "shapeMon-sem",
        "gyrMon-mean", "gyrMon-var", "gyrMon-sem"]
    assert table["whole"].to_list() == wholes[::-1]
    assert_allclose(table["gyrMon-mean"], [2.0, 1.0])
    assert_allclose(table["shapeMon-mean"], [4.0, 3.0])


def test_equilibrium_wholes_appends_the_topology_to_match_the_stamps(
        ens_dir, write_ensemble, alternating):
    """Stamps keyed by the whole's file name, ``...ens1.ring``, match the
    ensemble columns once the topology is appended, and the table keeps
    the stamps' names."""
    wholes = write_ensemble("gyrTMon", [alternating, alternating + 1])
    stamps = pd.DataFrame({"whole": [f"{whole}.ring" for whole in wholes]})

    table = equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps,
                               topology="ring")

    assert table["whole"].to_list() == [f"{whole}.ring" for whole in wholes]
    assert_allclose(table["gyrMon-mean"], [1.0, 2.0])


@pytest.mark.parametrize("stamps_wholes, match", [
    pytest.param(slice(0, 1), r"In the files only: \[.*ens2'\]",
                 id="stamps-lack-a-whole"),
    pytest.param(slice(0, 3), r"In the stamps only: \[.*ens3'\]",
                 id="files-lack-a-whole"),
])
def test_equilibrium_wholes_rejects_wholes_on_one_side_only(
        ens_dir, write_ensemble, alternating, stamps_wholes, match):
    """A whole in the ensemble files but not in the stamps, or in the
    stamps but not in the files, raises and is named, instead of a row
    being dropped or left empty."""
    wholes = write_ensemble("gyrTMon", [alternating, alternating + 1])
    stem = wholes[0][:-len("ens1")]
    stamps = pd.DataFrame({"whole": [f"{stem}ens{i}" for i in (1, 2, 3)]
                           [stamps_wholes]})

    with pytest.raises(ValueError, match=match):
        equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps)


def test_equilibrium_wholes_rejects_repeated_output_names(
        ens_dir, write_ensemble, alternating):
    """Two properties mapped to one output name would overwrite each
    other's columns, so they are rejected."""
    wholes = write_ensemble("gyrTMon", [alternating])
    stamps = pd.DataFrame({"whole": wholes})
    with pytest.raises(ValueError, match=r"more than once: \['gyrMon'\]"):
        equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon",
                                     "shapeTMon": "gyrMon"}, stamps)


@pytest.mark.parametrize("stamps, match", [
    pytest.param(pd.DataFrame({"n_frames": [256]}), "needs a 'whole'",
                 id="no-whole-column"),
    pytest.param(pd.DataFrame({"whole": ["w1", "w1"]}), "more than once",
                 id="repeated-whole"),
])
def test_equilibrium_wholes_rejects_bad_stamps(ens_dir, write_ensemble,
                                               alternating, stamps, match):
    """The stamps need a ``whole`` column that lists each whole once."""
    write_ensemble("gyrTMon", [alternating])
    with pytest.raises(ValueError, match=match):
        equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"}, stamps)


def test_equilibrium_wholes_property_without_files_raises(
        ens_dir, write_ensemble, alternating):
    """A property with no ensemble file raises and is named."""
    wholes = write_ensemble("gyrTMon", [alternating])
    stamps = pd.DataFrame({"whole": wholes})
    with pytest.raises(FileNotFoundError, match="'shapeTMon'"):
        equilibrium_wholes(ens_dir, {"shapeTMon": "shapeMon"}, stamps)


def test_equilibrium_wholes_empty_directory_raises(ens_dir):
    """An ``ens`` directory without CSV files raises."""
    with pytest.raises(FileNotFoundError, match="No file matches"):
        equilibrium_wholes(ens_dir, {"gyrTMon": "gyrMon"},
                           pd.DataFrame({"whole": ["w1"]}))


# Golden comparison with the research code, written blind


def _study_module():  # pragma: no cover - tier C helpers
    """The TransFociCub spec module. The study folders are not an installed
    package: import them from the repository root, as projects/conftest.py
    does for their own tests."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    return importlib.import_module("projects.trans_foci_cub.spec")


def _research_folder(private_data_dir, phase, suffix):  # pragma: no cover
    """The one ``TransFociCub/<phase>/*-bug-<suffix>`` folder of the research
    code's output under the data directory, or a skip."""
    parent = private_data_dir / "TransFociCub" / phase
    folders = sorted(parent.glob(f"*-bug-{suffix}"))
    if len(folders) != 1:
        pytest.skip(f"no single *-bug-{suffix} folder under {parent}")
    return folders[0]


def _bug_table(module):  # pragma: no cover
    """The study's analysis table of the ``bug`` group."""
    (table,) = [t for t in module.ANALYSIS if t.group == "bug"]
    return table


def _cap_autocorrelation_lags(table, n_frames):  # pragma: no cover
    """The table with each ``Autocorrelation`` hook's lags capped at
    ``n_frames - 1``. The research code's statsmodels ``acf`` clamped lags
    beyond the series silently, and PolyPhys's rejects them, so a fixture
    of short wholes needs the cap; the spec itself is unchanged."""
    return dataclasses.replace(table, properties=tuple(
        dataclasses.replace(prop, post_hook=Autocorrelation(
            min(prop.post_hook.nlags, n_frames - 1)))
        if isinstance(prop.post_hook, Autocorrelation) else prop
        for prop in table.properties))


def _assert_ens_avg_matches(produced, reference):  # pragma: no cover
    """Every ``*-ensAvg.csv`` of `produced` with a counterpart in `reference`
    equals it, to float64 round-off, and at least one does. The
    pair-distance artifacts and the ACF band are left out: they changed on
    purpose (criteria section 4), and the ACFs run on capped lags."""
    compared = 0
    for path in sorted(produced.glob("*-ensAvg.csv")):
        if "pairDist" in path.name or "-acf" in path.name:
            continue
        counterpart = reference / path.name
        if not counterpart.exists():
            continue
        pd.testing.assert_frame_equal(
            pd.read_csv(path), pd.read_csv(counterpart), rtol=1e-10,
            atol=1e-12, check_like=True)
        compared += 1
    assert compared > 0, f"no ensAvg file of {produced} has a counterpart"


@pytest.mark.private_data
def test_ens_avg_matches_the_research_code(  # pragma: no cover - needs data
        private_data_dir, tmp_path):
    """The ensemble averages of one TransFociCub space, reduced from the
    research code's probe output, match the research code's, for every
    property it reduced the same way.

    Expects ``<data>/TransFociCub/probes/<space>-bug-probe/`` (the probe
    phase of one space) and
    ``<data>/TransFociCub/analysis/<prefix>-bug-ensAvg/`` (the research
    code's ensemble averages). The probe folder is linked into a temporary
    tree so that nothing is written under the data directory.
    """
    module = _study_module()
    probe = _research_folder(private_data_dir, "probes", "probe")
    reference = _research_folder(private_data_dir, "analysis", "ensAvg")
    linked = tmp_path / "probes" / probe.name
    linked.parent.mkdir()
    linked.symlink_to(probe, target_is_directory=True)
    # The frame count from one whole's gyrTMon header, not its values.
    n_frames = np.load(next(probe.rglob("*-gyrTMon.npy")),
                       mmap_mode="r").shape[0]
    table = _cap_autocorrelation_lags(_bug_table(module), n_frames)

    reduce_space(linked, module.SPEC, table)

    produced = tmp_path / "analysis" / f"{probe.name.split('-')[0]}-bug-ensAvg"
    _assert_ens_avg_matches(produced, reference)


def _uncrowded(whole, wholes):  # pragma: no cover
    """The whole without crowders of the family and replicate of `whole`:
    the name that differs only in the crowder count and box side,
    ``nc<n>l<side>``, and has ``nc0``."""
    family = re.sub(r"nc\d+l[\d.]+", "", whole)
    (match,) = [other for other in wholes if "nc0l" in other
                and re.sub(r"nc\d+l[\d.]+", "", other) == family]
    return match


def _probe_tier_b_wholes(module, wholes, probe_dir):  # pragma: no cover
    """Probe each whole's tier B polymer dump (``polyphys/tests/data/probe/``)
    with the study's ``ProbeSpec``, opened with the data file of its whole
    without crowders, and save the chain-size observables, the foci pair
    quantities, and the stamps under the research names in
    ``probe_dir/<whole>.ring/``."""
    from polyphys.probe import (
        ChainSize, PairDistances, make_universe, write_stamps)

    data = Path(__file__).resolve().parents[1] / "data" / "probe"
    for whole in wholes:
        topology = _uncrowded(whole, wholes)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", "Reader has no dt information", UserWarning)
            universe = make_universe(
                data / topology / f"{topology}.ring.all.data.gz",
                data / whole / f"{whole}.ring.bug.last50.lammpstrj.gz")
        parser = module.SPEC.parser(f"{whole}.ring", "whole",
                                    module.SPEC.group)
        folder = probe_dir / f"{whole}.ring"
        folder.mkdir(parents=True)
        name = f"{parser.name}-{parser.group}"
        ChainSize.from_spec(universe, module.PROBE).run().save(folder, name)
        PairDistances.from_spec(universe, module.PROBE, parser).run().save(
            folder, name)
        write_stamps(parser, universe.trajectory.n_frames, folder)


@pytest.mark.private_data
def test_probe_and_reduce_match_the_research_ens_avg(  # pragma: no cover
        private_data_dir, tmp_path):
    """Probing the TransFociCub polymer dumps with the study's ``ProbeSpec``,
    saving under the research names, and reducing them with ``reduce_space``
    reproduces the research code's ensemble averages end to end. The wholes
    are those of the research probe folder; the dumps are their tier B
    copies, the same 50 frames the research probe read."""
    pytest.importorskip("MDAnalysis")
    module = _study_module()
    research = _research_folder(private_data_dir, "probes", "probe")
    reference = _research_folder(private_data_dir, "analysis", "ensAvg")
    wholes = sorted(path.name.removesuffix(".ring")
                    for path in research.iterdir() if path.is_dir())
    probe_dir = tmp_path / "probes" / research.name
    _probe_tier_b_wholes(module, wholes, probe_dir)
    # The tier B polymer dumps hold the wholes' last 50 frames.
    table = _cap_autocorrelation_lags(_bug_table(module), 50)

    reduce_space(probe_dir, module.SPEC, table)

    produced = (tmp_path / "analysis"
                / f"{research.name.split('-')[0]}-bug-ensAvg")
    _assert_ens_avg_matches(produced, reference)


@pytest.mark.private_data
@pytest.mark.parametrize("prop", ["gyrTMon", "shapeTMon", "asphericityTMon"])
def test_space_tseries_from_spec_matches_the_research_code(  # pragma: no cover
        private_data_dir, prop):
    """The space table of a TransFociCub time series, built from the
    research code's ensemble averages, equals the research code's space
    table: the same rows by state point and time index, the same
    statistics, time, and attributes, under the spec's renaming of the
    crowder volume fraction and the documented mapping of the short names
    (Usage, "Lineage names")."""
    spec = _study_module().SPEC
    ens_avg = _research_folder(private_data_dir, "analysis", "ensAvg")
    space_dir = _research_folder(private_data_dir, "analysis", "space")
    theirs = pd.read_csv(
        space_dir / f"{space_dir.name.split('-')[0]}-bug-{prop}-space.csv")
    theirs = theirs.rename(columns={"phi_c_bulk": "phi_bulk_c",
                                    "phi_c_bulk_round": "phi_bulk_c_round"})
    # A research short name, parsed at its lineage and rebuilt, gives the
    # current form.
    theirs = theirs.assign(**{
        lineage: theirs[lineage].map(
            lambda text, lineage=lineage: getattr(
                spec.parser(text, lineage, spec.group), lineage))
        for lineage in ("ensemble", "space")})
    key = ["ensemble_long", "t_index"]

    ours = space_tseries_from_spec(ens_avg, spec, prop)

    ours = ours.sort_values(key).reset_index(drop=True)
    theirs = theirs.sort_values(key).reset_index(drop=True)
    assert set(theirs.columns) <= set(ours.columns)
    # The statistics are float64 values read from the same CSV files; the
    # integer attributes compare exactly. The research table wrote some
    # attributes as integers that PolyPhys holds as floats, hence no dtype
    # check.
    pd.testing.assert_frame_equal(ours[theirs.columns], theirs, rtol=1e-10,
                                  atol=1e-12, check_dtype=False)


# The space stage from the spec


@pytest.fixture
def ens_avg_db(tmp_path):
    """An ``ensAvg`` directory with the gyrTMon time series and the rHistMon
    histogram of one ChainCyl state point with n = 1000 crowders of
    a = 2 in a cylinder of r = 10.5 and lz = 25."""
    database = tmp_path / "analysis" / "N200-bug-ensAvg"
    database.mkdir(parents=True)
    state = "N200r10.5lz25.0ac2.0nc1000dt0.005bdump2000adump5000"
    pd.DataFrame({"gyrTMon-mean": [1.0, 2.0]}).to_csv(
        database / f"{state}-bug-gyrTMon-ensAvg.csv", index=False)
    pd.DataFrame({"rHistMon-mean": [1.0, 3.0]}).to_csv(
        database / f"{state}-bug-rHistMon-ensAvg.csv", index=False)
    return database


@pytest.mark.parametrize("divisor, rounded", [
    pytest.param(0.001, 0.329, id="step-0.001"),
    pytest.param(0.025, 0.325, id="step-0.025"),
])
def test_space_tseries_from_spec_rounds_with_the_spec_divisor(
        ens_avg_db, cyl_spec, divisor, rounded):
    """The space table rounds phi_bulk_c to the spec's divisor."""
    spec = dataclasses.replace(cyl_spec, divisor=divisor)

    space = space_tseries_from_spec(ens_avg_db, spec, "gyrTMon")

    # phi_c = n pi a^3 / 6 / (pi Z (D - a)^2 / 4) = 2 n a^3 / (3 Z (D - a)^2)
    # with n = 1000, a = 2, Z = 2 lz = 50, D = 2 r - 1 = 20:
    # 16000 / 48600 = 0.32922, which rounds to 0.329 on a grid of 0.001 and
    # to 0.325 on a grid of 0.025.
    assert space["phi_bulk_c_round"].tolist() == [rounded, rounded]


def test_space_tseries_from_spec_matches_the_organizer_with_spec_fields(
        ens_avg_db, cyl_spec):
    """The spec-taking time-series table equals the organizer's table built
    from the spec's fields one by one."""
    expected = organizer.space_tseries(
        f"{ens_avg_db}/", "gyrTMon", cyl_spec.parser, cyl_spec.hierarchy,
        list(cyl_spec.physical_attrs), cyl_spec.group,
        divisor=cyl_spec.divisor, crowder_attr=cyl_spec.crowder_attr)

    space = space_tseries_from_spec(ens_avg_db, cyl_spec, "gyrTMon")

    pd.testing.assert_frame_equal(space, expected)


def test_space_hists_from_spec_matches_the_organizer_with_spec_fields(
        ens_avg_db, cyl_spec):
    """The spec-taking histogram table equals the organizer's table built
    from the spec's fields one by one, with the bin centers and the
    normalization passed through."""
    centers = np.array([0.5, 1.5])
    expected = organizer.space_hists(
        f"{ens_avg_db}/", "rHistMon", cyl_spec.parser, cyl_spec.hierarchy,
        list(cyl_spec.physical_attrs), cyl_spec.group, bin_center=centers,
        normalize=True, divisor=cyl_spec.divisor,
        crowder_attr=cyl_spec.crowder_attr)

    space = space_hists_from_spec(ens_avg_db, cyl_spec, "rHistMon",
                                  bin_center=centers, normalize=True)

    pd.testing.assert_frame_equal(space, expected)
