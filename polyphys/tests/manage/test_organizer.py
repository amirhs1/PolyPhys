"""Tests for :mod:`polyphys.manage.organizer`.

Most inputs come from the ``make_lineage_tree`` factory in
``polyphys/tests/conftest.py``: tiny ``.npy`` files named for the
``SumRuleCyl`` parser, under a directory whose name contains hyphens, as
every lineage directory name does. File ``(e, w, s)`` (ensemble, whole,
segment, each from 1; ``s = 0`` for a whole file) holds
``1000 * e + 100 * w + 10 * s + k`` in row ``k``, so every expected sum,
concatenation, and replicate statistic below is written down by hand from
that rule. Expected lineage names come from the parser applied to the
written filename, never from the strings the tree builds.

The space tests build a small ``analysis`` phase on disk; their expected
``t_index`` and ``time`` values come from the ``dt`` field and the dump
period of the group's trajectory (``adump`` for ``all``, ``bdump`` for
``bug``) encoded in the filenames. The stamp tests use the column layout of the
SciOps probe phase's stamp files: genealogy columns, parser attributes, and
``n_frames``.
"""
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from polyphys.manage.organizer import (
    children_stamps,
    create_fullname,
    ens_avg_from_bin_edge,
    ens_avg_from_df,
    ens_avg_from_ndarray,
    ens_from_bin_edge,
    ens_from_df,
    ens_from_mat_t,
    ens_from_vec,
    ensemble,
    ensemble_avg,
    find_unique_properties,
    make_database,
    parents_stamps,
    save_artifact,
    space_hists,
    space_tseries,
    whole_from_dist_mat_t,
    whole_from_file,
    whole_from_segments,
)
from polyphys.manage.parser import SumRuleCyl
from polyphys.manage.utils import round_to_multiple


# Fixtures: lineage trees and the names the parser gives their files


@pytest.fixture
def segment_tree(make_lineage_tree):
    """Return a segment tree: 2 ensembles x 2 wholes x 2 segments x 3 rows."""
    return make_lineage_tree(SumRuleCyl, group="bug", n_rows=3)


@pytest.fixture
def whole_tree(make_lineage_tree):
    """Return a whole tree: 2 ensembles x 2 wholes x 3 rows."""
    return make_lineage_tree(
        SumRuleCyl, group="bug", lineage="whole", n_rows=3
    )


def _parsed_names(tree, attr):
    """Map each tree file's (ensemble, whole) indices to a parsed name."""
    return {
        index[:2]: getattr(
            SumRuleCyl(str(path), tree.lineage, tree.group, ispath=True),
            attr,
        )
        for path, index in tree.index.items()
    }


@pytest.fixture
def segment_whole_names(segment_tree):
    """Return the parsed whole name of each (ensemble, whole) of the tree."""
    return _parsed_names(segment_tree, "whole")


@pytest.fixture
def whole_names(whole_tree):
    """Return the parsed whole name of each (ensemble, whole) in whole_tree."""
    return _parsed_names(whole_tree, "whole")


@pytest.fixture
def ensemble_names(whole_tree):
    """Return the parsed ensemble_long name of each ensemble in whole_tree."""
    return {
        e: name for (e, _), name in
        _parsed_names(whole_tree, "ensemble_long").items()
    }


@pytest.fixture
def write_segments(tmp_path):
    """Factory: write arrays as segments j01, j02, ... of one whole."""
    stem = (
        "N20epsilon5.0r5.5lz25.0sig1.0nc100dt0.005bdump1000adump5000ens1"
    )

    def _write(arrays):
        directory = tmp_path / "SumRuleCyl-bug-segment"
        directory.mkdir(exist_ok=True)
        paths = []
        for seg_id, array in enumerate(arrays, start=1):
            path = directory / f"{stem}.j{seg_id:02d}.bug.binEdge.npy"
            np.save(path, np.asarray(array))
            paths.append(str(path))
        return paths

    return _write


# Naming and persistence


def test_create_fullname_joins_with_hyphens():
    """A full name is the name, group, and property joined by hyphens."""
    assert create_fullname("w", "bug", "rhoMon") == "w-bug-rhoMon"


def test_save_artifact_writes_ndarray_to_npy(tmp_path):
    """An array is saved as <filename>.npy and loads back unchanged."""
    data = np.arange(3.0)
    save_artifact("w-bug-rhoMon", data, tmp_path)
    np.testing.assert_array_equal(np.load(tmp_path / "w-bug-rhoMon.npy"), data)


def test_save_artifact_writes_dataframe_to_csv_without_index(tmp_path):
    """A DataFrame is saved as <filename>.csv without its index."""
    data = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, index=[5, 6])
    save_artifact("w-bug-rhoMon", data, tmp_path)
    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "w-bug-rhoMon.csv"),
        data.reset_index(drop=True),
    )


def test_save_artifact_writes_dict_entries_by_measure(tmp_path):
    """Each '<prop>-<measure>' entry is saved as <filename>-<measure>.npy."""
    data = {"rhoMon-mean": np.array([1.0]), "rhoMon-var": np.array([2.0])}
    save_artifact("w-bug-rhoMon", data, tmp_path)
    assert sorted(os.listdir(tmp_path)) == [
        "w-bug-rhoMon-mean.npy", "w-bug-rhoMon-var.npy"
    ]
    np.testing.assert_array_equal(
        np.load(tmp_path / "w-bug-rhoMon-var.npy"), [2.0]
    )


def test_save_artifact_rejects_unsupported_type(tmp_path):
    """Data that is not an array, a DataFrame, or a dict raises TypeError."""
    with pytest.raises(TypeError, match="Unsupported data type list"):
        save_artifact("w-bug-rhoMon", [1.0, 2.0], tmp_path)


@pytest.mark.parametrize(
    "group, directory",
    [
        pytest.param("bug", "N20-bug-ens", id="with-group"),
        pytest.param(None, "N20-ens", id="without-group"),
    ],
)
def test_make_database_creates_sibling_phase_directory(
    tmp_path, group, directory
):
    """The new directory is <phase>/<prefix>-<group>-<stage> beside the old."""
    old = tmp_path / "probes" / "N20-bug-segment"
    new = make_database(old, "analysis", "ens", group)
    assert new == str(tmp_path / "analysis" / directory) + "/"
    assert Path(new).is_dir()


def test_make_database_returns_existing_directory(tmp_path):
    """Asking twice for the same directory returns the same path."""
    old = tmp_path / "probes" / "N20-bug-segment"
    first = make_database(old, "analysis", "ens", "bug")
    assert make_database(old, "analysis", "ens", "bug") == first


def test_make_database_logs_reuse_of_existing_directory(
    tmp_path, caplog, capsys
):
    """Reusing a directory logs an INFO record naming it and prints nothing."""
    old = tmp_path / "probes" / "N20-bug-segment"
    make_database(old, "analysis", "ens", "bug")
    caplog.set_level(logging.INFO, logger="polyphys.manage.organizer")

    make_database(old, "analysis", "ens", "bug")

    assert [(record.levelname, record.getMessage()) for record in
            caplog.records] == [
        ("INFO", f"Reusing the existing directory "
                 f"{tmp_path / 'analysis' / 'N20-bug-ens'}.")]
    assert capsys.readouterr().out == ""


def test_make_database_rejects_unknown_phase(tmp_path):
    """A phase outside the phase vocabulary raises ValueError."""
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        make_database(tmp_path / "probes" / "N20-bug-segment", "nope", "ens",
                      "bug")


# Building wholes


def test_whole_from_segments_tseries_concatenates_in_segment_order(
    segment_tree, segment_whole_names
):
    """A tseries whole is its segments' rows, j01 first, then j02."""
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "tseries")
    np.testing.assert_array_equal(
        wholes[segment_whole_names[(1, 2)]],
        [1210.0, 1211.0, 1212.0, 1220.0, 1221.0, 1222.0],
    )


def test_whole_from_segments_histogram_sums_segments(
    segment_tree, segment_whole_names
):
    """A histogram whole is the element-wise sum of its segments."""
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "histogram")
    # (1210 + 1220, 1211 + 1221, 1212 + 1222)
    np.testing.assert_array_equal(
        wholes[segment_whole_names[(1, 2)]], [2430.0, 2432.0, 2434.0]
    )


def test_whole_from_segments_keys_wholes_by_parsed_name(
    segment_tree, segment_whole_names
):
    """There is one whole per parsed whole name, whatever the directory."""
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "tseries")
    assert "-" in segment_tree.root.name
    assert list(wholes) == list(segment_whole_names.values())


def test_whole_from_segments_bin_edge_keeps_shared_edges(write_segments):
    """Identical segment edges reduce to one copy of them."""
    edges = np.array([0.0, 1.0, 3.0])
    paths = write_segments([edges, edges])
    wholes = whole_from_segments("binEdge", paths, SumRuleCyl, "bug",
                                 "bin_edge")
    (whole,) = wholes.values()
    np.testing.assert_array_equal(whole, edges)


def test_whole_from_segments_bin_edge_rejects_differing_edges(write_segments):
    """Segments of one whole with different edges raise ValueError."""
    paths = write_segments([[0.0, 1.0, 3.0], [0.0, 1.0, 4.0]])
    with pytest.raises(ValueError, match="Bin edges are not identical"):
        whole_from_segments("binEdge", paths, SumRuleCyl, "bug", "bin_edge")


def test_whole_from_segments_saves_each_whole(
    tmp_path, segment_tree, segment_whole_names
):
    """save_to writes each whole to <whole>-<group>-<prop>.npy."""
    out = tmp_path / "out"
    out.mkdir()
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "tseries", save_to=out)
    name = segment_whole_names[(2, 1)]
    np.testing.assert_array_equal(
        np.load(out / f"{name}-bug-gyrTMon.npy"), wholes[name]
    )


@pytest.mark.parametrize("relation", ["histogram", "tseries", "bin_edge"])
def test_whole_from_segments_rejects_descending_segment_ids(
    segment_tree, relation
):
    """Segments of a whole given in descending segment_id order raise."""
    paths = [str(p) for p in segment_tree.paths[1::-1]]
    with pytest.raises(ValueError, match=r"segment_id .*\[2, 1\]"):
        whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug", relation)


@pytest.mark.parametrize("relation", ["histogram", "tseries", "bin_edge"])
def test_whole_from_segments_rejects_duplicate_segment(segment_tree, relation):
    """A segment passed twice raises instead of being counted twice."""
    first = str(segment_tree.paths[0])
    with pytest.raises(ValueError, match=r"segment_id .*\[1, 1\]"):
        whole_from_segments(
            "gyrTMon", [first, first], SumRuleCyl, "bug", relation
        )


def test_whole_from_segments_accepts_interleaved_wholes(segment_tree):
    """Wholes interleaved segment by segment give the same wholes."""
    ordered = [str(p) for p in segment_tree.paths]
    # All j01 files first, then all j02 files; each whole stays in order.
    interleaved = ordered[0::2] + ordered[1::2]
    expected = whole_from_segments("gyrTMon", ordered, SumRuleCyl, "bug",
                                   "tseries")
    wholes = whole_from_segments("gyrTMon", interleaved, SumRuleCyl, "bug",
                                 "tseries")
    assert list(wholes) == list(expected)
    np.testing.assert_array_equal(
        np.stack(list(wholes.values())), np.stack(list(expected.values()))
    )


def test_whole_from_segments_rejects_unknown_relation(segment_tree):
    """A relation other than histogram, tseries, or bin_edge raises."""
    paths = [str(p) for p in segment_tree.paths]
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug", "nope")


def test_whole_from_file_loads_each_whole_by_parsed_name(
    whole_tree, whole_names
):
    """Each whole file loads under its parsed whole name."""
    wholes = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    assert "-" in whole_tree.root.name
    assert list(wholes) == list(whole_names.values())
    np.testing.assert_array_equal(
        wholes[whole_names[(2, 1)]], [2100.0, 2101.0, 2102.0]
    )


def test_whole_from_file_accepts_grouped_paths(whole_tree):
    """1-tuples of paths, as sort_filenames gives, load like bare paths."""
    grouped = list(zip(map(str, whole_tree.paths)))
    bare = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    wholes = whole_from_file(grouped, SumRuleCyl, "bug")
    assert list(wholes) == list(bare)
    np.testing.assert_array_equal(
        np.stack(list(wholes.values())), np.stack(list(bare.values()))
    )


def test_whole_from_file_rejects_unsupported_path_item():
    """An item that is neither a path nor a sequence of one raises."""
    with pytest.raises(TypeError, match="Unsupported path item"):
        whole_from_file([42], SumRuleCyl, "bug")


def test_whole_from_dist_mat_t_keys_outputs_by_parsed_name(
    whole_tree, whole_names
):
    """func runs once per whole and its three outputs key by whole name."""
    calls = []

    def func(whole_path, whole_info):
        calls.append(whole_path)
        frame = pd.DataFrame({"x": [0.0]})
        return frame, frame, frame

    grouped = list(zip(map(str, whole_tree.paths)))
    freqs, rdfs, tseries = whole_from_dist_mat_t(
        grouped, SumRuleCyl, "bug", func
    )
    expected = list(whole_names.values())
    assert list(freqs) == list(rdfs) == list(tseries) == expected
    assert calls == [str(p) for p in whole_tree.paths]


# Grouping wholes into ensembles


def test_ens_from_vec_puts_one_whole_per_column():
    """Each whole becomes a column, in the order of the dictionary."""
    name, ens_df = ens_from_vec(
        ("ens1", {"w2": np.array([1.0, 2.0]), "w1": np.array([3.0, 4.0])})
    )
    assert name == "ens1"
    pd.testing.assert_frame_equal(
        ens_df, pd.DataFrame({"w2": [1.0, 2.0], "w1": [3.0, 4.0]})
    )


def test_ens_from_mat_t_stacks_wholes_along_axis_0():
    """Matrices stack into (n_wholes, n_rows, n_cols) in dictionary order."""
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[5.0, 6.0], [7.0, 8.0]])
    name, stacked = ens_from_mat_t(("ens1", {"w1": a, "w2": b}))
    assert name == "ens1"
    np.testing.assert_array_equal(stacked, np.stack([a, b]))


def test_ens_from_df_averages_wholes_row_by_row():
    """Whole DataFrames average element-wise into one: (a + b) / 2."""
    a = pd.DataFrame({"x": [1.0, 2.0], "y": [3.0, 4.0]})
    b = pd.DataFrame({"x": [3.0, 6.0], "y": [5.0, 8.0]})
    name, ens_df = ens_from_df(("ens1", {"w1": a, "w2": b}))
    assert name == "ens1"
    pd.testing.assert_frame_equal(ens_df, (a + b) / 2)


def test_ens_from_bin_edge_returns_shared_edges():
    """Identical whole edges pass through as a copy of the shared array."""
    edges = np.array([0.0, 1.0, 2.0])
    name, out = ens_from_bin_edge(("ens1", {"w1": edges, "w2": edges.copy()}))
    assert name == "ens1"
    np.testing.assert_array_equal(out, edges)
    assert not np.shares_memory(out, edges)


@pytest.mark.parametrize(
    "other",
    [
        pytest.param(np.array([0.0, 1.0, 3.0]), id="different-values"),
        pytest.param(np.array([0.0, 1.0, 2.0, 3.0]), id="different-lengths"),
    ],
)
def test_ens_from_bin_edge_raises_on_mismatch(other):
    """Differing whole edges raise ValueError instead of being merged."""
    with pytest.raises(ValueError, match="Bin edges are not identical"):
        ens_from_bin_edge(
            ("ens1", {"w1": np.array([0.0, 1.0, 2.0]), "w2": other})
        )


def test_ensemble_groups_wholes_by_state_point(
    whole_tree, whole_names, ensemble_names
):
    """Wholes of two state points form two ensembles keyed by ensemble_long."""
    wholes = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector")
    assert list(ens) == [ensemble_names[1], ensemble_names[2]]
    assert list(ens[ensemble_names[2]].columns) == [
        whole_names[(2, 1)], whole_names[(2, 2)]
    ]


def test_ensemble_takes_bin_center_from_each_ensemble_edges(
    whole_tree, whole_names, ensemble_names
):
    """Each ensemble's bin_center comes from its own wholes' edges."""
    wholes = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    edges_1 = np.array([0.0, 1.0, 3.0, 6.0])
    edges_2 = np.array([0.0, 2.0, 4.0, 8.0])
    whole_edges = {
        whole_names[(1, 1)]: edges_1,
        whole_names[(1, 2)]: edges_1,
        whole_names[(2, 1)]: edges_2,
        whole_names[(2, 2)]: edges_2,
    }
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector",
                   whole_edges=whole_edges)
    # Analytic midpoints of each ensemble's edges
    np.testing.assert_allclose(
        ens[ensemble_names[1]]["bin_center"], [0.5, 2.0, 4.5]
    )
    np.testing.assert_allclose(
        ens[ensemble_names[2]]["bin_center"], [1.0, 3.0, 6.0]
    )


@pytest.mark.parametrize(
    "whole_type, first, second, expected",
    [
        pytest.param(
            "matrix",
            np.array([[1.0, 2.0]]), np.array([[3.0, 4.0]]),
            np.array([[[1.0, 2.0]], [[3.0, 4.0]]]),
            id="matrix-stacks",
        ),
        pytest.param(
            "dataframe",
            pd.DataFrame({"x": [1.0, 2.0]}), pd.DataFrame({"x": [3.0, 6.0]}),
            np.array([[2.0], [4.0]]),
            id="dataframe-averages",
        ),
        pytest.param(
            "bin_edge",
            np.array([0.0, 1.0]), np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            id="bin-edge-shares",
        ),
    ],
)
def test_ensemble_merges_by_whole_type(
    whole_names, ensemble_names, whole_type, first, second, expected
):
    """matrix stacks, dataframe averages, and bin_edge shares the wholes."""
    wholes = {whole_names[(1, 1)]: first, whole_names[(1, 2)]: second}
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", whole_type)
    np.testing.assert_array_equal(np.asarray(ens[ensemble_names[1]]),
                                  expected)


def test_ensemble_saves_each_ensemble(
    tmp_path, whole_tree, ensemble_names
):
    """save_to writes each ensemble to <ensemble_long>-<group>-<prop>.csv."""
    out = tmp_path / "out"
    out.mkdir()
    wholes = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector",
                   save_to=out)
    name = ensemble_names[2]
    pd.testing.assert_frame_equal(
        pd.read_csv(out / f"{name}-bug-gyrTMon.csv"), ens[name]
    )


def test_ensemble_rejects_unknown_whole_type(whole_tree):
    """A whole_type outside the four accepted values raises ValueError."""
    wholes = whole_from_file(whole_tree.paths, SumRuleCyl, "bug")
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "nope")


@pytest.fixture
def hist_wholes():
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
def hist_edges(hist_wholes):
    """Return the same non-uniform bin edges for every whole in hist_wholes.

    The edges are non-uniform so a wrong centre formula cannot pass by
    accident.
    """
    return {name: np.array([0.0, 1.0, 3.0, 6.0]) for name in hist_wholes}


def test_ensemble_vector_bin_center_is_edge_midpoints(hist_wholes,
                                                      hist_edges):
    """With whole_edges, a vector ensemble's bin_center is the midpoints."""
    ens = ensemble(
        "rhoMon", hist_wholes, SumRuleCyl, "bug", "vector",
        whole_edges=hist_edges,
    )
    (ens_df,) = ens.values()
    # Analytic midpoints of the edges [0, 1, 3, 6]
    np.testing.assert_allclose(ens_df["bin_center"], [0.5, 2.0, 4.5])
    pd.testing.assert_frame_equal(
        ens_df.drop(columns="bin_center"), pd.DataFrame(hist_wholes)
    )


def test_ensemble_without_whole_edges_has_no_bin_center(hist_wholes):
    """Without whole_edges, a vector ensemble has one column per whole only."""
    ens = ensemble("rhoMon", hist_wholes, SumRuleCyl, "bug", "vector")
    (ens_df,) = ens.values()
    assert list(ens_df.columns) == list(hist_wholes)


def test_ensemble_avg_keeps_bin_center_from_ensemble(hist_wholes,
                                                     hist_edges):
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
    whole_type, hist_wholes, hist_edges
):
    """whole_edges with any whole_type other than 'vector' raises."""
    with pytest.raises(ValueError, match="'whole_edges' is only valid"):
        ensemble(
            "rhoMon", hist_wholes, SumRuleCyl, "bug", whole_type,
            whole_edges=hist_edges,
        )


# Averaging over an ensemble


def test_ens_avg_from_ndarray_statistics_along_axis_0():
    """Two wholes a, b give mean (a+b)/2, var (a-b)^2/2, sem |a-b|/2."""
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[2.0, 5.0], [3.0, 8.0]])
    out = ens_avg_from_ndarray("gyrTMon", np.stack([a, b]))
    assert list(out) == ["gyrTMon-mean", "gyrTMon-var", "gyrTMon-sem"]
    # Float64 round-off on values of order 10
    np.testing.assert_allclose(out["gyrTMon-mean"], (a + b) / 2, atol=1e-12)
    np.testing.assert_allclose(out["gyrTMon-var"], (a - b) ** 2 / 2,
                               atol=1e-12)
    np.testing.assert_allclose(out["gyrTMon-sem"], np.abs(a - b) / 2,
                               atol=1e-12)


def test_ens_avg_from_df_leaves_input_unchanged():
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


def test_ens_avg_from_df_statistics():
    """Two wholes a, b give mean (a+b)/2, var (a-b)^2/2, sem |a-b|/2."""
    a = np.array([4.0, 5.0])
    b = np.array([5.0, 7.0])
    ens_data = pd.DataFrame({"bin_center": [1.0, 2.0], "w1": a, "w2": b})
    out = ens_avg_from_df("density", ens_data, ["bin_center"])
    np.testing.assert_allclose(out["density-mean"], (a + b) / 2)
    np.testing.assert_allclose(out["density-var"], (a - b) ** 2 / 2)
    np.testing.assert_allclose(out["density-sem"], np.abs(a - b) / 2)
    np.testing.assert_array_equal(out["bin_center"], [1.0, 2.0])


def test_ens_avg_from_bin_edge_copies_1d_edges():
    """A single 1D edge array, as ensemble_avg passes it, is copied."""
    edges = np.array([0.0, 0.5, 1.0])
    out = ens_avg_from_bin_edge(edges)
    np.testing.assert_array_equal(out, edges)
    assert not np.shares_memory(out, edges)


def test_ens_avg_from_bin_edge_reduces_equal_sequence():
    """A sequence of equal edge arrays reduces to one copy of them."""
    edges = np.array([0.0, 0.5, 1.0])
    np.testing.assert_array_equal(
        ens_avg_from_bin_edge([edges, edges.copy()]), edges
    )


def test_ens_avg_from_bin_edge_rejects_differing_sequence():
    """A sequence of differing edge arrays raises ValueError."""
    edges = np.array([0.0, 0.5, 1.0])
    with pytest.raises(ValueError, match="Bin edges are not identical"):
        ens_avg_from_bin_edge([edges, np.array([0.0, 0.5, 2.0])])


def test_ensemble_avg_bin_edge_keeps_edges():
    """ensemble_avg with 'bin_edge' returns each ensemble's edges intact."""
    edges = np.array([0.0, 0.5, 1.0, 1.5])
    out = ensemble_avg("rhoMon", {"ens1": edges}, "all", "bin_edge")
    np.testing.assert_array_equal(out["ens1"], edges)


def test_ensemble_avg_ndarray_averages_each_ensemble():
    """ensemble_avg with 'ndarray' gives each ensemble's statistics dict."""
    a = np.array([[1.0, 2.0]])
    b = np.array([[3.0, 6.0]])
    out = ensemble_avg("gyrTMon", {"ens1": np.stack([a, b])}, "bug",
                       "ndarray")
    np.testing.assert_allclose(out["ens1"]["gyrTMon-mean"], [[2.0, 4.0]])


def test_ensemble_avg_rejects_unknown_ens_type():
    """An ens_type outside dataframe, ndarray, and bin_edge raises."""
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        ensemble_avg("gyrTMon", {"ens1": np.zeros(2)}, "bug", "nope")


@pytest.mark.parametrize(
    "ensembles, ens_type, filenames",
    [
        pytest.param(
            {"ens1": pd.DataFrame({"w1": [1.0], "w2": [3.0]})},
            "dataframe",
            ["ens1-bug-rhoMon-ensAvg.csv"],
            id="dataframe-csv",
        ),
        pytest.param(
            {"ens1": np.zeros((2, 1, 1))},
            "ndarray",
            [
                "ens1-bug-rhoMon-ensAvg-mean.npy",
                "ens1-bug-rhoMon-ensAvg-sem.npy",
                "ens1-bug-rhoMon-ensAvg-var.npy",
            ],
            id="ndarray-npy-per-measure",
        ),
        pytest.param(
            {"ens1": np.array([0.0, 1.0])},
            "bin_edge",
            ["ens1-bug-rhoMon-ensAvg.npy"],
            id="bin-edge-npy",
        ),
    ],
)
def test_ensemble_avg_saves_under_ens_avg_name(
    tmp_path, ensembles, ens_type, filenames
):
    """save_to writes each result under <ensemble>-<group>-<prop>-ensAvg."""
    ensemble_avg("rhoMon", ensembles, "bug", ens_type, save_to=tmp_path)
    assert sorted(os.listdir(tmp_path)) == filenames


# Stamps


@pytest.fixture
def segment_stamps():
    """Return stamps of 2 wholes x 3 segments of one state point.

    n_frames is 10, 20, 30 and gyr is 1, 2, 3 for the segments of w1, and
    gyr is 4, 5, 6 for those of w2.
    """
    return pd.DataFrame({
        "lineage_name": ["w1.j01", "w1.j02", "w1.j03",
                         "w2.j01", "w2.j02", "w2.j03"],
        "segment": ["w1.j01", "w1.j02", "w1.j03",
                    "w2.j01", "w2.j02", "w2.j03"],
        "whole": ["w1", "w1", "w1", "w2", "w2", "w2"],
        "ensemble_long": ["s1"] * 6,
        "ensemble": ["e1"] * 6,
        "space": ["N20"] * 6,
        "group": ["bug"] * 6,
        "segment_id": [1, 2, 3, 1, 2, 3],
        "ensemble_id": [1, 1, 1, 2, 2, 2],
        "gyr": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        "n_frames": [10, 20, 30, 10, 20, 30],
    })


@pytest.fixture
def whole_stamps():
    """Return stamps of 2 state points x 2 wholes, each of 60 frames.

    The second state point keeps replicates ens1 and ens3, as if ens2 had
    failed, so the number of wholes differs from the last ensemble_id.
    """
    return pd.DataFrame({
        "lineage_name": ["w1", "w2", "w3", "w4"],
        "segment": ["N/A"] * 4,
        "whole": ["w1", "w2", "w3", "w4"],
        "ensemble_long": ["s1", "s1", "s2", "s2"],
        "ensemble": ["e1", "e1", "e2", "e2"],
        "space": ["N20"] * 4,
        "group": ["bug"] * 4,
        "segment_id": [0, 0, 0, 0],
        "ensemble_id": [1, 2, 1, 3],
        "n_frames": [60, 60, 60, 60],
    })


@pytest.fixture
def write_stamps(tmp_path):
    """Factory: write each row of a stamps frame to its own CSV file."""

    def _write(stamps):
        directory = tmp_path / "stamps"
        directory.mkdir()
        paths = []
        for row in range(len(stamps)):
            path = directory / f"{stamps.loc[row, 'lineage_name']}-stamps.csv"
            stamps.iloc[[row]].to_csv(path, index=False)
            paths.append((str(path),))
        return paths

    return _write


def test_children_stamps_concatenates_in_input_order(
    write_stamps, segment_stamps
):
    """Segment stamp files combine row by row, in the order given."""
    paths = write_stamps(segment_stamps)
    combined = children_stamps(paths, "segment")
    pd.testing.assert_frame_equal(combined, segment_stamps)


def test_children_stamps_whole_drops_segment_columns(
    write_stamps, whole_stamps
):
    """Whole stamps lose 'segment' and 'segment_id', with a UserWarning."""
    paths = write_stamps(whole_stamps)
    with pytest.warns(UserWarning, match="'segment' and 'segment_id'"):
        combined = children_stamps(paths, "whole")
    pd.testing.assert_frame_equal(
        combined, whole_stamps.drop(columns=["segment", "segment_id"])
    )


def test_children_stamps_whole_without_segment_columns_is_unchanged(
    write_stamps, whole_stamps, capsys
):
    """Whole stamps without segment columns pass through, printing nothing."""
    stamps = whole_stamps.drop(columns=["segment", "segment_id"])
    paths = write_stamps(stamps)

    combined = children_stamps(paths, "whole")

    pd.testing.assert_frame_equal(combined, stamps)
    assert capsys.readouterr().out == ""


def test_children_stamps_saves_combined_stamps(
    tmp_path, write_stamps, segment_stamps
):
    """save_to writes the result as <space>-<group>-<lineage>-stamps.csv."""
    paths = write_stamps(segment_stamps)
    combined = children_stamps(paths, "segment", save_to=tmp_path)
    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "N20-bug-segment-stamps.csv"), combined
    )


def test_children_stamps_rejects_unknown_lineage():
    """A lineage other than segment or whole raises ValueError."""
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        children_stamps([], "nope")


def test_parents_stamps_segment_counts_segments_and_sums_frames(
    segment_stamps
):
    """Segment stamps give one row per whole with n_segments and n_frames."""
    parents = parents_stamps(segment_stamps, "segment")
    assert parents["whole"].tolist() == ["w1", "w2"]
    assert parents["n_segments"].tolist() == [3, 3]
    # 10 + 20 + 30 frames per whole
    assert parents["n_frames"].tolist() == [60, 60]


def test_parents_stamps_whole_counts_wholes_per_ensemble(whole_stamps):
    """Whole stamps give one row per state point with n_ensembles."""
    parents = parents_stamps(whole_stamps, "whole")
    assert parents["ensemble_long"].tolist() == ["s1", "s2"]
    assert parents["n_ensembles"].tolist() == [2, 2]
    assert parents["n_frames"].tolist() == [60, 60]
    assert {"segment", "segment_id", "whole"}.isdisjoint(parents.columns)


def test_parents_stamps_applies_given_aggregation(segment_stamps):
    """properties replaces the default 'last' aggregation of a column."""
    parents = parents_stamps(segment_stamps, "segment",
                             properties={"gyr": "mean"})
    # Means of 1, 2, 3 and of 4, 5, 6
    assert parents["gyr"].tolist() == pytest.approx([2.0, 5.0])


@pytest.mark.parametrize(
    "stamps_fixture, lineage, filename",
    [
        pytest.param("segment_stamps", "segment", "N20-bug-stamps-whole.csv",
                     id="segment-to-whole"),
        pytest.param("whole_stamps", "whole", "N20-bug-stamps-ensAvg.csv",
                     id="whole-to-ensAvg"),
    ],
)
def test_parents_stamps_saves_under_parent_stage(
    request, tmp_path, stamps_fixture, lineage, filename
):
    """save_to names the parent stamps after the parents' stage."""
    stamps = request.getfixturevalue(stamps_fixture)
    parents = parents_stamps(stamps, lineage, save_to=tmp_path)
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / filename), parents)


@pytest.mark.parametrize(
    "stamps_fixture, lineage",
    [
        pytest.param("segment_stamps", "segment", id="segment"),
        pytest.param("whole_stamps", "whole", id="whole"),
    ],
)
def test_parents_stamps_without_lineage_name_drops_child_column(
    request, stamps_fixture, lineage
):
    """Stamps without 'lineage_name' still lose the child lineage column."""
    stamps = request.getfixturevalue(stamps_fixture).drop(
        columns="lineage_name")

    parents = parents_stamps(stamps, lineage)

    assert lineage not in parents.columns


def test_parents_stamps_without_lineage_name_prints_nothing(
    segment_stamps, capsys
):
    """An absent 'lineage_name' column is skipped without printing."""
    stamps = segment_stamps.drop(columns="lineage_name")

    parents_stamps(stamps, "segment")

    assert capsys.readouterr().out == ""


def test_parents_stamps_rejects_unknown_lineage(segment_stamps):
    """A lineage other than segment or whole raises ValueError."""
    with pytest.raises(ValueError, match="'nope' is an invalid option"):
        parents_stamps(segment_stamps, "nope")


# Unique properties


@pytest.fixture
def property_files(tmp_path):
    """Return a glob over four ensemble-averaged files of one ensemble."""
    for name in ("N20-bug-gyrTMon-ensAvg.csv",
                 "N20-bug-gyrTMon-acf-ensAvg.csv",
                 "N20-bug-rhoMon-ensAvg.csv",
                 "N20-bug-transSizeTMon-ensAvg.csv"):
        (tmp_path / name).touch()
    return str(tmp_path / "*")


def test_find_unique_properties_splits_properties_and_measures(
    property_files
):
    """Properties and property-measures come back sorted and separate."""
    props, measures = find_unique_properties(property_files, 2, ["-ensAvg"])
    assert props == ["gyrTMon", "rhoMon", "transSizeTMon"]
    assert measures == ["gyrTMon-acf"]


def test_find_unique_properties_drops_named_properties(property_files):
    """Names in drop_properties are left out of the result."""
    props, _ = find_unique_properties(
        property_files, 2, ["-ensAvg"], drop_properties=["transSizeTMon"]
    )
    assert props == ["gyrTMon", "rhoMon"]


# Space level: ensemble-averaged files of two state points that differ only
# in the number of crowders. Both dump every 5000 steps with dt=0.005.


@pytest.fixture
def space_states():
    """Return the ensemble_long names of the two synthetic state points."""
    return [
        "N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000adump5000",
        "N200epsilon5.0r10.5lz25.0sig2.0nc2000dt0.005bdump2000adump5000",
    ]


@pytest.fixture
def make_ens_avg_space(tmp_path, space_states):
    """Factory: write one ensAvg CSV per state point, plus a decoy file.

    Takes the property, one '-mean' column per state point, and the group
    (default 'all'); '-var' and '-sem' are constant. Returns the ensAvg
    directory, ending with '/'.
    """

    def _make(prop, means, group="all"):
        ens_avg_dir = tmp_path / "analysis" / f"N200-{group}-ensAvg"
        ens_avg_dir.mkdir(parents=True)
        for state, mean in zip(space_states, means):
            pd.DataFrame({
                f"{prop}-mean": mean,
                f"{prop}-var": np.full(len(mean), 0.1),
                f"{prop}-sem": np.full(len(mean), 0.05),
            }).to_csv(ens_avg_dir / f"{state}-{group}-{prop}-ensAvg.csv",
                      index=False)
        # A file of another property in the same directory must be ignored.
        pd.DataFrame({"other-mean": [1.0]}).to_csv(
            ens_avg_dir / f"{space_states[0]}-{group}-other-ensAvg.csv",
            index=False,
        )
        return str(ens_avg_dir) + "/"

    return _make


@pytest.fixture
def state_attr(space_states):
    """Factory: an attribute of each state point, repeated n_rows times."""

    def _values(attr, n_rows):
        parsers = [SumRuleCyl(s, "ensemble_long", "all") for s in space_states]
        return np.repeat([getattr(p, attr) for p in parsers], n_rows)

    return _values


def test_space_tseries_concatenates_ensembles(make_ens_avg_space,
                                              state_attr):
    """space_tseries reads every '-<prop>-ensAvg.csv' file and no other."""
    database = make_ens_avg_space(
        "gyrTMon", [[0.0, 1.0, 2.0], [10.0, 11.0, 12.0]]
    )
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
                          "all")
    assert list(space.index) == list(range(6))
    assert "other-mean" not in space.columns
    # Files are sorted by name, so the nc1000 rows come first.
    np.testing.assert_allclose(
        space["gyrTMon-mean"], [0.0, 1.0, 2.0, 10.0, 11.0, 12.0]
    )
    np.testing.assert_allclose(space["dcyl"], state_attr("dcyl", 3))


def test_space_tseries_scales_t_index_and_time(make_ens_avg_space):
    """For all, t_index is frame index times adump; time is t_index * dt."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", [], "all")
    frames = np.array([0, 1, 2, 0, 1, 2])
    # adump = 5000 and dt = 0.005 in the filenames
    np.testing.assert_array_equal(space["t_index"], frames * 5000)
    np.testing.assert_allclose(space["time"], frames * 5000 * 0.005)


def test_space_tseries_scales_bug_by_bdump(make_ens_avg_space):
    """For bug, t_index is frame index times bdump, not adump."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2,
                                  group="bug")
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", [], "bug")
    frames = np.array([0, 1, 2, 0, 1, 2])
    # bdump = 2000 and dt = 0.005 in the filenames
    np.testing.assert_array_equal(space["t_index"], frames * 2000)
    np.testing.assert_allclose(space["time"], frames * 2000 * 0.005)


def test_space_tseries_rejects_dump_key_argument(make_ens_avg_space):
    """A dump keyword passed by position fails instead of binding divisor."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    with pytest.raises(TypeError, match="takes 6 positional arguments"):
        space_tseries(database, "gyrTMon", SumRuleCyl, "N*", [], "all",
                      "adump")


def test_space_tseries_adds_phi_bulk_c_by_default(make_ens_avg_space,
                                                  state_attr):
    """The crowder fraction comes from phi_bulk_c without being listed."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
                          "all")
    assert "phi_c_bulk" not in space.columns
    np.testing.assert_allclose(space["phi_bulk_c"],
                               state_attr("phi_bulk_c", 3))


def test_space_tseries_rounds_crowder_with_round_to_multiple(
    make_ens_avg_space
):
    """phi_c_bulk_round is round_to_multiple of the crowder column."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", [], "all",
                          divisor=0.05, round_to=2)
    expected = space["phi_bulk_c"].apply(round_to_multiple, args=(0.05, 2))
    np.testing.assert_array_equal(space["phi_c_bulk_round"], expected)


def test_space_tseries_reads_named_crowder_attr(make_ens_avg_space,
                                                state_attr):
    """crowder_attr selects the parser attribute that is rounded."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", [], "all",
                          crowder_attr="phi_bulk_m")
    assert "phi_bulk_c" not in space.columns
    np.testing.assert_allclose(space["phi_bulk_m"],
                               state_attr("phi_bulk_m", 3))


def test_space_tseries_saves_space_csv(tmp_path, make_ens_avg_space):
    """is_save writes the returned frame to the analysis space directory."""
    database = make_ens_avg_space("gyrTMon", [[0.0, 1.0, 2.0]] * 2)
    space = space_tseries(database, "gyrTMon", SumRuleCyl, "N*", ["dcyl"],
                          "all", is_save=True)
    saved = tmp_path / "analysis" / "N200-all-space" / \
        "N200-all-gyrTMon-space.csv"
    pd.testing.assert_frame_equal(pd.read_csv(saved), space)


def test_space_hists_concatenates_and_normalizes(make_ens_avg_space,
                                                 state_attr):
    """space_hists adds bin centres and each normalized histogram sums to 1."""
    database = make_ens_avg_space(
        "rhoMon", [[1.0, 2.0, 1.0], [2.0, 4.0, 2.0]]
    )
    bin_center = np.array([0.5, 1.5, 2.5])
    space = space_hists(database, "rhoMon", SumRuleCyl, "N*", ["dcyl"],
                        "all", bin_center=bin_center, normalize=True)
    assert len(space) == 6
    assert "other-mean" not in space.columns
    np.testing.assert_allclose(space["bin_center"], np.tile(bin_center, 2))
    np.testing.assert_allclose(space["bin_center-norm"],
                               np.tile(bin_center / 2.5, 2))
    # [1, 2, 1] / 4 and [2, 4, 2] / 8
    np.testing.assert_allclose(space["rhoMon-norm"],
                               [0.25, 0.5, 0.25, 0.25, 0.5, 0.25])
    np.testing.assert_allclose(space["dcyl"], state_attr("dcyl", 3))
    np.testing.assert_allclose(space["phi_bulk_c"],
                               state_attr("phi_bulk_c", 3))


def test_space_hists_zero_histogram_normalizes_to_zero(make_ens_avg_space):
    """An all-zero histogram warns and gets a zero normalized column."""
    database = make_ens_avg_space("rhoMon", [[0.0, 0.0], [0.0, 0.0]])
    with pytest.warns(UserWarning, match="All values are zero"):
        space = space_hists(database, "rhoMon", SumRuleCyl, "N*", [], "all",
                            normalize=True)
    np.testing.assert_array_equal(space["rhoMon-norm"], np.zeros(4))


def test_space_hists_saves_space_csv(tmp_path, make_ens_avg_space):
    """is_save writes the returned frame to the analysis space directory."""
    database = make_ens_avg_space("rhoMon", [[1.0, 2.0, 1.0]] * 2)
    space = space_hists(database, "rhoMon", SumRuleCyl, "N*", [], "all",
                        is_save=True)
    saved = tmp_path / "analysis" / "N200-all-space" / \
        "N200-all-rhoMon-space.csv"
    pd.testing.assert_frame_equal(pd.read_csv(saved), space)


# Integration along the lineage: segment -> whole -> ensemble -> ensAvg
# -> space, on one segment tree.


def test_lineage_tseries_ensemble_average_matches_tree_values(
    segment_tree, ensemble_names
):
    """Two wholes 100 apart average to 1000e + 150 + 10s + k, var 5000."""
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "tseries")
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector")
    avg = ensemble_avg("gyrTMon", ens, "bug", "dataframe")
    ens_1 = avg[ensemble_names[1]]
    np.testing.assert_allclose(
        ens_1["gyrTMon-mean"],
        [1160.0, 1161.0, 1162.0, 1170.0, 1171.0, 1172.0],
    )
    # var (ddof=1) of two values 100 apart is 100**2 / 2; sem is
    # sqrt(5000 / 2). Float64 round-off on values of order 1e3.
    np.testing.assert_allclose(ens_1["gyrTMon-var"], 5000.0, rtol=1e-12)
    np.testing.assert_allclose(ens_1["gyrTMon-sem"], 50.0, rtol=1e-12)


def test_lineage_saved_ens_avg_feeds_space_tseries(tmp_path, segment_tree):
    """ensAvg files saved by ensemble_avg are the files space_tseries reads."""
    ens_avg_dir = tmp_path / "analysis" / "N20-bug-ensAvg"
    ens_avg_dir.mkdir(parents=True)
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "tseries")
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector")
    ensemble_avg("gyrTMon", ens, "bug", "dataframe", save_to=ens_avg_dir)
    space = space_tseries(str(ens_avg_dir) + "/", "gyrTMon", SumRuleCyl,
                          "N*", [], "bug")
    # 2 ensembles x 6 frames; ensemble 2 is 1000 above ensemble 1
    np.testing.assert_allclose(
        space["gyrTMon-mean"],
        [1160.0, 1161.0, 1162.0, 1170.0, 1171.0, 1172.0,
         2160.0, 2161.0, 2162.0, 2170.0, 2171.0, 2172.0],
    )
    # bdump = 1000, the dump period of the bug trajectory, in the tree's
    # filenames
    np.testing.assert_array_equal(
        space["t_index"], np.tile(np.arange(6), 2) * 1000
    )


def test_lineage_histogram_conserves_total_counts(segment_tree):
    """Histogram counts summed over a space equal those over its segments."""
    paths = [str(p) for p in segment_tree.paths]
    wholes = whole_from_segments("gyrTMon", paths, SumRuleCyl, "bug",
                                 "histogram")
    ens = ensemble("gyrTMon", wholes, SumRuleCyl, "bug", "vector")
    avg = ensemble_avg("gyrTMon", ens, "bug", "dataframe")
    # Ensemble mean times 2 wholes, summed over bins and ensembles
    total = 2 * pd.concat(avg.values())["gyrTMon-mean"].sum()
    # Sum of 1000e + 100w + 10s + k over e, w, s in {1, 2} and k in
    # {0, 1, 2}: 24 * (1500 + 150 + 15) + 8 * (0 + 1 + 2)
    assert total == pytest.approx(39984.0, rel=1e-12)
