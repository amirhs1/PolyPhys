"""Tests for :mod:`polyphys.analyze.pipeline`.

Most inputs are lineage trees from the ``make_lineage_tree`` factory in
``conftest.py``, written into one probe directory per test, one property per
subfolder. A tree's values are ``1000 e + 100 w + 10 s + k``
(``lineage_tree_value``), so two wholes of one ensemble differ by exactly
100 in every row and the mean, variance, and standard error across them are
``w1 + 50``, ``5000``, and ``50`` by hand. Bin edges, distance matrices, and
two-dimensional wholes are written over the tree's files by hand. The hooks
are also tested on their own with a ``HookContext``; the density kernel is a
stub that records its calls. The golden comparison with the research code's
output is written blind and runs only with the maintainer's data.
"""
import importlib
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.clusters import foci_pair_info
from polyphys.analyze.correlations import acf
from polyphys.analyze.pipeline import (
    AnalysisTable,
    Autocorrelation,
    Derived,
    FociPairDistances,
    HookContext,
    LocalDensities,
    PropertySpec,
    TimeAverage,
    edge,
    half_contour,
    histogram,
    histogram_2d,
    matrix_tseries,
    reduce_space,
    reduce_stamps,
    tseries,
)
from polyphys.manage.parser import SumRuleCyl, TransFociCub
from polyphys.manage.projects import ProjectSpec
from polyphys.tests.conftest import lineage_tree_value, write_lineage_tree

# A Pandas4Warning marks a deprecated pandas call; fail on it.
pytestmark = pytest.mark.filterwarnings("error::pandas.errors.Pandas4Warning")

# The glob that lists every file of a probe directory built by write_tree:
# <probe>/<property>/<Parser>-<group>-<stage>/<file>.
HIERARCHY = "*/*/*"


# Fixtures: specs, a probe directory, and tree writers


@pytest.fixture
def cyl_spec():
    """A SumRuleCyl spec whose crowder size is ``dcrowd``."""
    return ProjectSpec(
        parser=SumRuleCyl, group="bug", divisor=0.025,
        physical_attrs=("space",),
        species_sizes={"Mon": "dmon", "Crd": "dcrowd"},
        tseries_properties=(),
    )


@pytest.fixture
def foci_spec():
    """A TransFociCub spec: a ring of 20 small and 2 large monomers."""
    return ProjectSpec(
        parser=TransFociCub, group="bug", divisor=0.025,
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


@pytest.mark.parametrize("kwargs, match", [
    pytest.param(dict(name="", species="Mon", relation="tseries",
                      whole_type="vector"), "must not be empty",
                 id="empty-name"),
    pytest.param(dict(name="gyrT", species="Mon", relation="series",
                      whole_type="vector"), "relation must be",
                 id="unknown-relation"),
    pytest.param(dict(name="gyrT", species="Mon", relation="tseries",
                      whole_type="frame"), "whole_type must be",
                 id="unknown-whole-type"),
    pytest.param(dict(name="xEdge", species="", relation="bin_edge",
                      whole_type="vector"), "do not match",
                 id="edge-relation-vector-whole"),
    pytest.param(dict(name="xEdge", species="", relation="tseries",
                      whole_type="bin_edge"), "do not match",
                 id="tseries-relation-edge-whole"),
    pytest.param(dict(name="gyrT", species="Mon", relation="tseries",
                      whole_type="vector", edge_prop="gyrEdgeMon"),
                 "edge_prop is only valid", id="edges-on-tseries"),
    pytest.param(dict(name="xyHist", species="Crd", relation="histogram",
                      whole_type="matrix", edge_prop="xyEdgeCrd"),
                 "edge_prop is only valid", id="edges-on-matrix"),
])
def test_property_spec_rejects_inconsistent_kinds(kwargs, match):
    """A property whose relation, shape, and edges disagree is rejected."""
    with pytest.raises(ValueError, match=match):
        PropertySpec(**kwargs)


def test_constructors_build_the_kinds_the_studies_use():
    """Each named constructor sets the property, relation, shape, and edges
    of its kind."""
    assert (tseries("gyrT", "Mon").prop, tseries("gyrT", "Mon").relation,
            tseries("gyrT", "Mon").whole_type) == ("gyrTMon", "tseries",
                                                   "vector")
    hist = histogram("z", "Crd")
    assert (hist.prop, hist.relation, hist.edge_prop) == (
        "zHistCrd", "histogram", "zEdgeCrd")
    assert (histogram_2d("xz", "Dna").prop,
            histogram_2d("xz", "Dna").whole_type) == ("xzHistDna", "matrix")
    assert (edge("y").prop, edge("y").relation, edge("y").has_edges) == (
        "yEdge", "bin_edge", False)
    assert matrix_tseries("principalT", "Mon").whole_type == "matrix"


def test_analysis_table_rejects_bad_lineage_and_repeated_properties():
    """A table needs a known lineage and lists each artifact once."""
    with pytest.raises(ValueError, match="lineage must be"):
        AnalysisTable("bug", "ensemble", False, ())
    with pytest.raises(ValueError, match=r"more than once: \['gyrTMon'\]"):
        AnalysisTable("bug", "whole", False,
                      (tseries("gyrT", "Mon"), tseries("gyrT", "Mon")))


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
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)

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
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1,
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
    tree = write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=3)
    names = _parsed(tree, SumRuleCyl, "whole")
    (ens_name,) = set(_parsed(tree, SumRuleCyl, "ensemble_long").values())
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
    tree = write_tree(SumRuleCyl, "gyrTMon", n_ensembles=1, n_wholes=1,
                      n_segments=2, n_rows=2)
    (whole_name,) = set(_parsed(tree, SumRuleCyl, "whole").values())
    table = AnalysisTable("bug", "segment", False,
                          (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    whole = np.load(
        _analysis(probe_dir, "wholeSim") / f"{whole_name}-bug-gyrTMon.npy")
    assert_array_equal(whole, np.concatenate([lineage_tree_value(1, 1, 1, 2),
                                              lineage_tree_value(1, 1, 2, 2)]))


def test_segments_are_ordered_naturally(probe_dir, cyl_spec):
    """Segments j2 and j10 join in numeric order, not in text order."""
    stem = "N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000ens1"
    folder = probe_dir / "gyrTMon" / "SumRuleCyl-bug-segment"
    folder.mkdir(parents=True)
    for seg_id in (2, 10):
        np.save(folder / f"{stem}.j{seg_id}.bug.gyrTMon.npy",
                np.array([float(seg_id)]))
    table = AnalysisTable("bug", "segment", False,
                          (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    whole = np.load(
        _analysis(probe_dir, "wholeSim") / f"{stem}-bug-gyrTMon.npy")
    assert_array_equal(whole, [2.0, 10.0])


def test_segment_histograms_sum_and_ensembles_carry_bin_centres(
        probe_dir, cyl_spec, write_tree, write_edges):
    """A segment-lineage histogram sums its segments element-wise, and its
    ensemble gets the midpoints of the shared edges as ``bin_center``."""
    tree = write_tree(SumRuleCyl, "rHistMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=3)
    write_edges(tree, "rEdgeMon", np.array([0.0, 1.0, 2.0, 4.0]))
    names = _parsed(tree, SumRuleCyl, "whole")
    (ens_name,) = set(_parsed(tree, SumRuleCyl, "ensemble_long").values())
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
    tree = write_tree(SumRuleCyl, "rHistMon", lineage="whole", n_ensembles=1,
                      n_wholes=2, n_rows=3)
    write_edges(tree, "rEdgeMon", None, per_whole={
        1: np.array([0.0, 1.0, 2.0, 3.0]), 2: np.array([0.0, 1.0, 2.0, 3.5])})
    names = _parsed(tree, SumRuleCyl, "whole")
    (ens_name,) = set(_parsed(tree, SumRuleCyl, "ensemble_long").values())
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
    tree = write_tree(SumRuleCyl, "rHistMon", lineage="whole", n_ensembles=1,
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
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False,
                          (histogram("r", "Crdtwo"),))

    with pytest.raises(FileNotFoundError, match="rHistCrdtwo"):
        reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY,
                     save=False)


def test_default_hierarchy_is_the_spec_dir_file_layout(probe_dir, cyl_spec,
                                                       write_tree):
    """Without ``hierarchy`` the files are listed under
    ``<spec.hierarchy>/<spec.hierarchy>``, the research code's layout."""
    tree = write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)
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
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1,
               n_rows=3)
    write_tree(SumRuleCyl, "xgyrTMon", lineage="whole", n_ensembles=1,
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
    tree = write_tree(SumRuleCyl, "gyrTMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=2)
    names = _parsed(tree, SumRuleCyl, "whole")
    (ens_name,) = set(_parsed(tree, SumRuleCyl, "ensemble_long").values())
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
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    analysis = probe_dir.parent.parent / "analysis"
    assert sorted(path.name for path in analysis.iterdir()) == [
        "N-bug-ens", "N-bug-ensAvg"]


def test_save_false_writes_nothing(probe_dir, cyl_spec, write_tree):
    """With ``save=False`` no analysis folder appears."""
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)
    table = AnalysisTable("bug", "whole", False, (tseries("gyrT", "Mon"),))

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY, save=False)

    assert not (probe_dir.parent.parent / "analysis").exists()


def test_reduce_space_logs_each_property(probe_dir, cyl_spec, write_tree,
                                         caplog):
    """Progress is logged at INFO, one record per property."""
    write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1)
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
    tree = write_tree(SumRuleCyl, "xyHistCrd", lineage="whole",
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
    tree = write_tree(SumRuleCyl, "xEdge", group="all", lineage="whole",
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
    tree = write_tree(SumRuleCyl, "gyrTMon", lineage="whole", n_ensembles=1,
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
    names = _parsed(tree, SumRuleCyl, "whole")
    (ens_name,) = set(_parsed(tree, SumRuleCyl, "ensemble_long").values())
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
    tree = write_tree(SumRuleCyl, "rHistCrd", group="all", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=3)
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    write_edges(tree, "rEdgeCrd", edges)
    names = _parsed(tree, SumRuleCyl, "whole")
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
        assert diameter == 1.0   # sig1.0 in the written names
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
    tree = write_tree(SumRuleCyl, "rHistCrd", group="all", lineage="whole",
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


# The foci pair distances


@pytest.fixture
def foci_context(foci_spec):
    """A context for ``distMatTFoci`` of the TransFociCub ring."""
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
        infos={"w": TransFociCub(
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


def test_dataframe_wholes_reduce_to_long_ensembles_and_row_statistics(
        probe_dir, foci_spec, write_tree):
    """Through ``reduce_space``, the pair histograms of two wholes give an
    ensemble file with every whole's rows and an ensemble average with the
    mean, variance, and standard error of each pair's counts, row by row."""
    tree = write_tree(TransFociCub, "distMatTFoci", lineage="whole",
                      n_ensembles=1, n_wholes=2, n_rows=1)
    separations = {1: [[1.5], [1.5]], 2: [[2.5], [2.5]]}
    for path, (_, whole, _) in tree.index.items():
        np.save(path, _distance_frames([0, 11], separations[whole]))
    names = _parsed(tree, TransFociCub, "whole")
    (ens_name,) = set(_parsed(tree, TransFociCub, "ensemble_long").values())
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
    stem = "N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000"
    wholes = [f"{stem}ens1", f"{stem}ens2"]
    _write_stamps(probe_dir, "SumRuleCyl", "bug", wholes, n_segments=2)
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
    stem = "N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000"
    wholes = [f"{stem}ens1", f"{stem}ens2"]
    _write_stamps(probe_dir, "SumRuleCyl", "bug", wholes, n_segments=1)
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
    stem = "N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000"
    _write_stamps(probe_dir, "SumRuleCyl", "bug", [f"{stem}ens1"],
                  n_segments=2)
    table = AnalysisTable("bug", "segment", True, ())

    reduce_space(probe_dir, cyl_spec, table, hierarchy=HIERARCHY)

    assert (_analysis(probe_dir, "wholeSim")
            / "N20ac1.0-bug-segment-stamps.csv").exists()
    assert (_analysis(probe_dir, "ens") / "N20ac1.0-bug-stamps-whole.csv"
            ).exists()
    assert (_analysis(probe_dir, "ensAvg") / "N20ac1.0-bug-stamps-ensAvg.csv"
            ).exists()


# End to end


def test_segment_tree_end_to_end(probe_dir, cyl_spec, write_tree,
                                 write_edges):
    """One segment-lineage run of a time series with its autocorrelation, a
    histogram with density profiles, and a time-averaged histogram gives
    the hand-composed averages of the segment, whole, and ensemble steps."""
    gyr = write_tree(SumRuleCyl, "gyrTMon", n_ensembles=1, n_wholes=2,
                     n_segments=2, n_rows=3)
    hist = write_tree(SumRuleCyl, "rHistMon", n_ensembles=1, n_wholes=2,
                      n_segments=2, n_rows=3)
    write_edges(hist, "rEdgeMon", np.arange(4.0))
    clusters = write_tree(SumRuleCyl, "clustersHistTMon", n_ensembles=1,
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
    (ens_name,) = set(_parsed(gyr, SumRuleCyl, "ensemble_long").values())

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
    (names,) = [set(_parsed(hist, SumRuleCyl, "whole").values())]
    assert sorted(path.name for path in _analysis(probe_dir, "wholeSim")
                  .iterdir()) == sorted(
        f"{name}-bug-{prop}.npy" for name in names for prop in (
            "gyrTMon", "gyrTMon-acf", "rHistMon", "rEdgeMon", "rRhoMon",
            "rPhiMon", "clustersHistTMon", "clustersHistMon"))


# Golden comparison with the research code, written blind


@pytest.mark.private_data
def test_ens_avg_matches_the_research_code(  # pragma: no cover - needs data
        private_data_dir, tmp_path):
    """The ensemble averages of one TransFociCub space match the research
    code's, for every property the research code reduced the same way.

    Expects ``<data>/TransFociCub/probes/<space>-bug-probe/`` (the probe
    phase of one space) and
    ``<data>/TransFociCub/analysis/<prefix>-bug-ensAvg/`` (the research
    code's ensemble averages). The probe folder is linked into
    a temporary tree so that nothing is written under the data directory.
    The pair-distance artifacts and the ACF band are left out: they changed
    on purpose (criteria section 4).
    """
    # The study folders are not an installed package: import them from the
    # repository root, as projects/conftest.py does for their own tests.
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    spec_module = importlib.import_module("projects.trans_foci_cub.spec")
    tables, spec = spec_module.ANALYSIS, spec_module.SPEC

    study = private_data_dir / "TransFociCub"
    probes = sorted((study / "probes").glob("*-bug-probe"))
    reference_dirs = sorted((study / "analysis").glob("*-bug-ensAvg"))
    if not probes or not reference_dirs:
        pytest.skip(f"no TransFociCub probe and ensAvg folders under {study}")
    probe, reference = probes[0], reference_dirs[0]
    linked = tmp_path / "probes" / probe.name
    linked.parent.mkdir()
    linked.symlink_to(probe, target_is_directory=True)
    (table,) = [t for t in tables if t.group == "bug"]

    reduce_space(linked, spec, table)

    produced = tmp_path / "analysis" / f"{probe.name.split('-')[0]}-bug-ensAvg"
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
