"""Tests for :mod:`polyphys.analyze.measurer`.

Every expected value is hand arithmetic on the positions or bounds written in
the test, or a closed form: the spherical-segment and sphere-sphere
intersection volumes from Weisstein, MathWorld, "Spherical Segment" and
"Sphere-Sphere Intersection", and the sample mean, the unbiased variance, and
the standard error of the mean for the statistics. Property tests say so in
their docstrings. Closed-form values are compared at ``rel=1e-12``, far above
the float64 round-off of the few operations behind them (about 1e-16).
"""
import warnings
from functools import partial

import pytest
import numpy as np
from hypothesis import given, strategies as st
from hypothesis.extra.numpy import arrays
from numpy.testing import assert_allclose, assert_array_equal

from polyphys.analyze.measurer import (
    apply_pbc_orthogonal,
    pair_distance,
    end_to_end,
    transverse_size,
    max_distance,
    fsd,
    simple_stats,
    sem,
    spherical_segment,
    sphere_sphere_intersection,
    create_bin_edge_and_hist,
    fixedsize_bins,
    radial_cyl_histogram,
    axial_histogram,
    radial_histogram,
    azimuth_cyl_histogram,
    planar_cartesian_histogram
)
from polyphys.manage.types import AxisT


@pytest.fixture
def positions_3d():
    """Three atom positions in 3D, for tests where their values are
    incidental."""
    return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]])


def test_apply_pbc_orthogonal_sets_lengths_and_inverses_on_named_axes():
    """Each named axis gets its box length and its inverse; the other axes
    keep their values."""
    lengths, inverses = apply_pbc_orthogonal(
        np.zeros(3), np.zeros(3), {0: 10.0, 2: 20.0}
    )
    # float64 round-off of one division; atol for the untouched zeros:
    assert_allclose(lengths, [10.0, 0.0, 20.0], rtol=1e-12, atol=1e-12)
    assert_allclose(inverses, [0.1, 0.0, 0.05], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize(
    "length,error,match",
    [
        pytest.param(0.0, ZeroDivisionError, "cannot be zero", id="zero"),
        pytest.param(-5.0, ValueError, "must be positive", id="negative"),
    ],
)
def test_apply_pbc_orthogonal_rejects_nonpositive_length(length, error, match):
    """A box length that is zero or negative is rejected."""
    with pytest.raises(error, match=match):
        apply_pbc_orthogonal(np.zeros(3), np.zeros(3), {1: length})


def test_pair_distance_is_coordinate_difference_without_pbc():
    """Without periodic boundaries the pair vector is the second position
    minus the first."""
    positions = np.array([[1.0, 2.0, 3.0], [4.0, 2.0, 3.0]])
    result = pair_distance(positions)
    # float64 round-off of one subtraction; atol for the zero components:
    assert_allclose(result, [3.0, 0.0, 0.0], rtol=1e-12, atol=1e-12)


def test_pair_distance_applies_minimum_image():
    """A separation of 8.5 along a periodic axis of length 8 folds to 0.5,
    the minimum image."""
    positions = np.array([[1.0, 2.0, 3.0], [9.5, 2.0, 3.0]])
    result = pair_distance(positions, pbc={0: 8.0})
    # float64 round-off of a subtraction and a product; atol for the zeros:
    assert_allclose(result, [0.5, 0.0, 0.0], rtol=1e-12, atol=1e-12)


def test_pair_distance_rejects_other_than_two_positions():
    """A positions array that does not hold exactly two atoms is rejected."""
    with pytest.raises(ValueError, match="only works for two atoms"):
        pair_distance(np.array([[1.0, 2.0, 3.0]]))


def test_end_to_end_of_straight_chain_is_its_length():
    """Three beads one unit apart on a line are two units end to end."""
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    result = end_to_end(positions)
    # float64 round-off of centring and a norm:
    assert result == pytest.approx(2.0, rel=1e-12)


def test_end_to_end_rejects_single_atom():
    """A chain of fewer than two atoms has no end-to-end distance."""
    with pytest.raises(ValueError, match="at least two atoms"):
        end_to_end(np.array([[0.0, 0.0, 0.0]]))


def test_transverse_size_is_twice_mean_distance_from_axis():
    """Two atoms at (0, 0) and (2, 2) in the y-z plane are each sqrt(2) from
    their centre, so the transverse size about x is 2 sqrt(2)."""
    positions = np.array([[1, 0, 0], [1, 2, 2]])
    result = transverse_size(positions, axis=0)
    # float64 round-off of centring and a norm:
    assert result == pytest.approx(2 * np.sqrt(2), rel=1e-12)


@pytest.mark.parametrize(
    "measure,positions,axis",
    [
        pytest.param(
            transverse_size, np.array([[1, 0, 0], [1, 2, 2]]), 3,
            id="transverse-size",
        ),
        pytest.param(fsd, np.array([[1, 2, 3]]), 5, id="fsd"),
    ],
)
def test_size_measures_reject_axis_beyond_positions(measure, positions, axis):
    """An axis index at or beyond the number of coordinates is rejected."""
    with pytest.raises(IndexError, match="out of bounds"):
        measure(positions, axis=axis)


def test_max_distance_is_extent_along_each_axis():
    """The extent along each axis is the spread of that coordinate."""
    positions = np.array([[0, 0, 0], [4, 2, 6]])
    result = max_distance(positions)
    # float64 round-off of centring:
    assert_allclose(result, [4.0, 2.0, 6.0], rtol=1e-12)


def test_fsd_is_spread_along_axis():
    """Feret's diameter along y is the spread of the y coordinates, 8 - 2."""
    positions = np.array([[1, 2, 3], [4, 8, 6]])
    result = fsd(positions, axis=1)
    # float64 round-off of centring:
    assert result == pytest.approx(6.0, rel=1e-12)


@pytest.mark.parametrize(
    "key,expected",
    [
        pytest.param("energy_mean", 2.0, id="mean"),
        pytest.param("energy_var", 1.0, id="unbiased-variance"),
        pytest.param("energy_sem", 1 / np.sqrt(3), id="sem"),
    ],
)
def test_simple_stats_gives_mean_variance_and_sem(key, expected):
    """The statistics of [1, 2, 3] are mean 2, unbiased variance 1, and
    standard error 1 / sqrt(3), under keys prefixed by the entity."""
    result = simple_stats("energy", np.array([1.0, 2.0, 3.0]))
    # float64 round-off of a few sums and a square root:
    assert result[key] == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "statistic",
    [
        pytest.param(partial(simple_stats, "energy"), id="simple-stats"),
        pytest.param(sem, id="sem"),
    ],
)
def test_statistics_reject_empty_data(statistic):
    """Statistics of an empty array are rejected rather than returned as
    NaN."""
    with pytest.raises(ValueError, match="must not be empty"):
        statistic(np.array([]))


def test_sem_is_sample_std_over_root_n():
    """The standard error of [1, 2, 3] is its unbiased standard deviation, 1,
    over sqrt(3)."""
    result = sem(np.array([1.0, 2.0, 3.0]))
    # float64 round-off of a few sums and a square root:
    assert result == pytest.approx(1 / np.sqrt(3), rel=1e-12)


@pytest.mark.parametrize(
    "r,a,b,expected",
    [
        pytest.param(3, 1, 2, 20 * np.pi / 3, id="slab"),
        pytest.param(3, -3, 3, 36 * np.pi, id="full-sphere"),
        pytest.param(3, 4, 5, 0.0, id="both-planes-outside"),
    ],
)
def test_spherical_segment_volume(r, a, b, expected):
    """The volume between two planes cutting a sphere is
    pi (r^2 (b - a) - (b^3 - a^3) / 3), with the planes clipped to the
    sphere, and zero when both lie outside it on one side."""
    result = spherical_segment(r, a, b)
    # float64 round-off of a few products; atol for the zero case:
    assert result == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_spherical_segment_rejects_nonpositive_radius():
    """A sphere of negative radius is rejected."""
    with pytest.raises(ValueError, match="must be positive"):
        spherical_segment(-3, 1, 2)


@pytest.mark.parametrize(
    "r1,r2,d,expected",
    [
        pytest.param(3, 4, 2, 725 * np.pi / 24, id="partial-overlap"),
        pytest.param(3, 4, 10, 0.0, id="no-overlap"),
        pytest.param(3, 4, 1, 36 * np.pi, id="smaller-inside-larger"),
        pytest.param(3, 0, 1, 0.0, id="second-radius-zero"),
        pytest.param(0, 4, 1, 0.0, id="first-radius-zero"),
    ],
)
def test_sphere_sphere_intersection_volume(r1, r2, d, expected):
    """The lens volume follows MathWorld's closed form; it is zero for
    disjoint or zero-radius spheres and the smaller sphere's volume when one
    contains the other."""
    result = sphere_sphere_intersection(r1, r2, d)
    # float64 round-off of a few products; atol for the zero cases:
    assert result == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_create_bin_edge_and_hist_unit_edges():
    """Unit bins over [0, 5] have edges at the integers 0 to 5."""
    edges, _ = create_bin_edge_and_hist(1.0, 0.0, 5.0)
    # float64 round-off of np.linspace over 5 unit steps:
    assert_allclose(edges, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], rtol=0, atol=1e-12)


def test_create_bin_edge_and_hist_histogram_is_zero_per_bin():
    """The histogram has one zero count per bin."""
    _, hist = create_bin_edge_and_hist(1.0, 0.0, 5.0)
    assert_array_equal(hist, np.zeros(5, dtype=np.int64))


def test_create_bin_edge_and_hist_rejects_inverted_range():
    """A lower bound above the upper bound is rejected."""
    with pytest.raises(ValueError, match="must be less than"):
        create_bin_edge_and_hist(1.0, 5.0, 0.0)


@pytest.mark.parametrize(
    "bin_type,expected_n_bins",
    [
        pytest.param("ordinary", 5, id="ordinary"),
        pytest.param("nonnegative", 5, id="nonnegative"),
        pytest.param("periodic", 5, id="periodic"),
    ],
)
def test_fixedsize_bins_unit_bins_over_five(bin_type, expected_n_bins):
    """Unit bins over [0, 5] give five bins and six edges for every bin
    type."""
    result = fixedsize_bins(1.0, 0.0, 5.0, bin_type=bin_type)
    assert result["n_bins"] == expected_n_bins
    assert len(result["bin_edges"]) == expected_n_bins + 1


def test_fixedsize_bins_rejects_inverted_range():
    """A lower bound above the upper bound is rejected."""
    with pytest.raises(ValueError, match="must be less than"):
        fixedsize_bins(1.0, 5.0, 0.0)


def test_fixedsize_bins_rejects_unknown_bin_type():
    """A bin type other than ordinary, nonnegative, or periodic is
    rejected."""
    with pytest.raises(ValueError, match="is an invalid option"):
        fixedsize_bins(1.0, 0.0, 5.0, bin_type="logarithmic")


@pytest.mark.parametrize(
    "make_edges",
    [
        pytest.param(
            lambda path: create_bin_edge_and_hist(
                1.0, 0.0, 5.0, output=path
            )[0],
            id="create-bin-edge-and-hist",
        ),
        pytest.param(
            lambda path: fixedsize_bins(
                1.0, 0.0, 5.0, save_bin_edges=path
            )["bin_edges"],
            id="fixedsize-bins",
        ),
    ],
)
def test_bin_edges_saved_to_npy_equal_returned_edges(tmp_path, make_edges):
    """The edges saved to the named .npy file are the edges returned (a
    round-trip property)."""
    path = tmp_path / "edges"
    edges = make_edges(str(path))
    # .npy stores float64 bit for bit, so the round trip is exact:
    assert_array_equal(np.load(tmp_path / "edges.npy"), edges)


@pytest.mark.parametrize(
    "make_array,expected_dtype",
    [
        pytest.param(
            lambda: fixedsize_bins(1.0, 0.0, 2.0)["collector"],
            np.int64,
            id="fixedsize-collector",
        ),
        pytest.param(
            lambda: fixedsize_bins(1.0, 0.0, 2.0)["collector_std"],
            np.float64,
            id="fixedsize-collector-std",
        ),
        pytest.param(
            lambda: create_bin_edge_and_hist(1.0, 0.0, 2.0)[1],
            np.int64,
            id="create-bin-edge-and-hist",
        ),
    ],
)
def test_histogram_collectors_are_64_bit(make_array, expected_dtype):
    """The empty collectors are 64-bit, so per-frame counts accumulated in
    place cannot overflow a 16-bit type."""
    assert make_array().dtype == expected_dtype


@pytest.mark.parametrize(
    "make_array,power,expected",
    [
        pytest.param(
            lambda: fixedsize_bins(1.0, 0.0, 2.0)["collector"],
            1,
            [60000, 0],
            id="fixedsize-collector",
        ),
        pytest.param(
            lambda: fixedsize_bins(1.0, 0.0, 2.0)["collector_std"],
            2,
            [2.0 * 30000**2, 0.0],
            id="fixedsize-collector-std",
        ),
        pytest.param(
            lambda: create_bin_edge_and_hist(1.0, 0.0, 2.0)[1],
            1,
            [60000, 0],
            id="create-bin-edge-and-hist",
        ),
    ],
)
def test_histogram_collectors_accumulate_past_int16_limit(
    make_array, power, expected
):
    """Two frames of 30000 counts in the first bin add to 60000, past the
    int16 limit of 32767, without wrapping."""
    collector = make_array()
    frame_counts = np.array([30000, 0])
    collector += frame_counts**power
    collector += frame_counts**power
    # integer counts below 2**53 are exact in float64, so no tolerance:
    assert_allclose(collector, expected, rtol=0, atol=0)


@pytest.mark.parametrize(
    "bin_size,lmin,lmax,expected_edges",
    [
        pytest.param(
            0.3, -6.0, 6.0, -6.0 + 0.3 * np.arange(41), id="multiple"
        ),
        pytest.param(
            0.7, 0.0, 2.0, np.array([0.0, 0.7, 1.4, 2.1]), id="non-multiple"
        ),
        pytest.param(
            0.3, 0.0, 2.1, 0.3 * np.arange(8), id="multiple-ratio-rounds-up"
        ),
    ],
)
def test_create_bin_edge_and_hist_edges_have_bin_size_width(
    bin_size, lmin, lmax, expected_edges
):
    """Edges step by bin_size from lmin to the first multiple covering lmax,
    and a range within round-off of a whole number of bins gets no extra
    bin."""
    edges, hist = create_bin_edge_and_hist(bin_size, lmin, lmax)
    # float64 round-off on at most 41 steps of order 1:
    np.testing.assert_allclose(edges, expected_edges, rtol=0, atol=1e-12)
    assert len(hist) == len(expected_edges) - 1


@pytest.mark.parametrize(
    "bin_size,lmin,lmax",
    [
        pytest.param(0.3, 0.0, 2.1, id="ratio-rounds-up"),
        pytest.param(0.03, 0.0, 15.9, id="ratio-rounds-down"),
        pytest.param(0.3, -6.0, 6.0, id="ratio-exact"),
    ],
)
def test_create_bin_edge_and_hist_exact_multiple_ends_at_lmax(
    bin_size, lmin, lmax
):
    """When the range is a whole number of bins up to round-off, the first and
    last edges are exactly lmin and lmax."""
    edges, _ = create_bin_edge_and_hist(bin_size, lmin, lmax)
    assert edges[0] == lmin
    assert edges[-1] == lmax


@pytest.mark.parametrize(
    "bin_type,bin_size,lmin,lmax,expected_n_bins",
    [
        pytest.param("ordinary", 0.3, 0.0, 2.1, 7, id="ordinary-rounds-up"),
        pytest.param(
            "ordinary", 0.03, 0.0, 15.9, 530, id="ordinary-rounds-down"
        ),
        pytest.param(
            "nonnegative", 0.3, 0.0, 2.1, 7, id="nonnegative-rounds-up"
        ),
        pytest.param("periodic", 0.3, 0.0, 2.1, 7, id="periodic-rounds-up"),
    ],
)
def test_fixedsize_bins_exact_multiple_tiles_range(
    bin_type, bin_size, lmin, lmax, expected_n_bins
):
    """A range that is a whole number of bins up to float round-off gets
    exactly that many bins, with the first and last edges at lmin and lmax."""
    result = fixedsize_bins(bin_size, lmin, lmax, bin_type=bin_type)
    assert result["n_bins"] == expected_n_bins
    assert result["bin_edges"][0] == lmin
    assert result["bin_edges"][-1] == lmax


@pytest.mark.parametrize(
    "bin_size,lmin,lmax,expected_n_bins",
    [
        pytest.param(0.3, -6.0, 6.0, 40, id="step-0.3"),
        pytest.param(np.pi / 36, -np.pi, np.pi, 72, id="azimuth-5-degrees"),
    ],
)
def test_fixedsize_bins_periodic_spans_one_period(
    bin_size, lmin, lmax, expected_n_bins
):
    """A period that is a multiple of bin_size gives exactly one period of
    bins, ending at lmax, with no warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        result = fixedsize_bins(bin_size, lmin, lmax, bin_type="periodic")
    assert result["n_bins"] == expected_n_bins
    assert len(result["bin_edges"]) == expected_n_bins + 1
    assert len(result["collector"]) == expected_n_bins
    assert result["bin_edges"][0] == lmin
    assert result["bin_edges"][-1] == lmax


def test_fixedsize_bins_nonnegative_negative_lmin_raises():
    """A negative lower bound is rejected for nonnegative bins, since clipping
    it to zero would leave the range below zero uncovered."""
    with pytest.raises(ValueError, match="must be nonnegative"):
        fixedsize_bins(1.0, -1.0, 2.0, bin_type="nonnegative")


def test_fixedsize_bins_periodic_non_multiple_raises():
    """A period that is not a multiple of bin_size is rejected."""
    with pytest.raises(ValueError, match="multiple of 'bin_size'"):
        fixedsize_bins(0.7, 0.0, 2.0, bin_type="periodic")


@pytest.mark.parametrize(
    "lmin,lmax,expected_edges",
    [
        pytest.param(0.1, 2.5, [0.0, 1.0, 2.0, 3.0], id="lmin-below-dl"),
        pytest.param(0.25, 2.75, [0.0, 1.0, 2.0, 3.0], id="lmin-equals-dl"),
        pytest.param(0.0, 2.5, [0.0, 1.0, 2.0, 3.0], id="lmin-zero"),
        pytest.param(0.5, 3.0, [0.25, 1.25, 2.25, 3.25], id="no-clipping"),
    ],
)
def test_fixedsize_bins_nonnegative_edges(lmin, lmax, expected_edges):
    """Nonnegative bins keep width bin_size when the lower end is clipped
    to zero."""
    result = fixedsize_bins(1.0, lmin, lmax, bin_type="nonnegative")
    # float64 round-off on 3 unit steps:
    np.testing.assert_allclose(
        result["bin_edges"], expected_edges, rtol=0, atol=1e-12
    )
    assert result["range"] == pytest.approx(
        (expected_edges[0], expected_edges[-1]), abs=1e-12
    )


@pytest.mark.parametrize(
    "bin_type,bin_size,lmin,lmax",
    [
        pytest.param("ordinary", 0.3, -6.0, 6.0, id="ordinary-multiple"),
        pytest.param("ordinary", 0.7, 0.0, 2.0, id="ordinary-non-multiple"),
        pytest.param(
            "ordinary", 0.3, 0.0, 2.1, id="ordinary-ratio-rounds-up"
        ),
        pytest.param(
            "ordinary", 0.03, 0.0, 15.9, id="ordinary-ratio-rounds-down"
        ),
        pytest.param("nonnegative", 1.0, 0.1, 2.5, id="nonnegative-clipped"),
        pytest.param("nonnegative", 0.3, 0.0, 6.0, id="nonnegative-zero"),
        pytest.param("periodic", 0.3, -6.0, 6.0, id="periodic-step-0.3"),
        pytest.param(
            "periodic", np.pi / 36, -np.pi, np.pi, id="periodic-azimuth"
        ),
    ],
)
def test_fixedsize_bins_edges_cover_range_with_bin_size_width(
    bin_type, bin_size, lmin, lmax
):
    """Every bin type returns n_bins + 1 edges of width bin_size whose first
    and last edges enclose lmin and lmax exactly, and a collector of n_bins
    slots (a property test)."""
    result = fixedsize_bins(bin_size, lmin, lmax, bin_type=bin_type)
    edges = result["bin_edges"]
    n_bins = result["n_bins"]
    assert len(edges) == n_bins + 1
    assert len(result["collector"]) == n_bins
    # float64 round-off on at most 530 steps of order 1:
    np.testing.assert_allclose(
        np.diff(edges), np.full(n_bins, bin_size), rtol=0, atol=1e-12
    )
    # no slack: a position at exactly lmin or lmax must fall inside the edges
    assert edges[0] <= lmin
    assert edges[-1] >= lmax


def test_radial_histogram_bins_distances_from_origin():
    """Distances from the origin are binned; the last edge is inclusive."""
    positions = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    edges = np.array([0, 1, 2, 3])
    # distances 1, 1, 2: two in [1, 2), one in [2, 3].
    result = radial_histogram(positions, edges)
    assert_array_equal(result, np.array([0, 2, 1]))


def test_radial_cyl_histogram_bins_distances_from_axis():
    """Distances from the y axis are binned in the x-z plane."""
    positions = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3], [0, 3, 4]])
    edges = np.array([0, 1, 2, 3, 4])
    # distances from the y axis: 1, 0, 3, 4; the last edge is inclusive.
    hist = radial_cyl_histogram(positions, edges, dim=1)
    assert_array_equal(hist, np.array([1, 1, 0, 2]))


def test_axial_histogram_bins_coordinates_along_axis():
    """Coordinates along the x axis are binned."""
    positions = np.array([
        [1, 2, 3],
        [2, 2, 3],
        [3, 2, 3],
        [4, 2, 3]
    ])
    edges = np.array([0, 1, 2, 3, 4, 5])
    hist = axial_histogram(positions, edges, dim=0)
    assert_array_equal(hist, np.array([0, 1, 1, 1, 1]))


def test_azimuth_cyl_histogram_bins_angles_about_axis():
    """Azimuth angles about the z axis are binned over [-pi, pi]."""
    pos = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0]])  # x-y plane
    edges = np.linspace(-np.pi, np.pi, 5)
    # angles 0, pi, pi/2; bins [-pi, -pi/2), [-pi/2, 0), [0, pi/2),
    # [pi/2, pi], the last edge inclusive.
    result = azimuth_cyl_histogram(pos, edges, dim=2)
    assert_array_equal(result, np.array([0, 0, 1, 2]))


def test_planar_cartesian_histogram_bins_in_plane_normal_to_axis():
    """Positions are binned in the x-y plane when the axis is z."""
    positions = np.array([
        [1, 1, 0],  # falls into bin (0,0)
        [2, 2, 0],  # falls into bin (1,1)
        [3, 3, 0]   # falls into bin (1,1)
    ])
    edges = [np.array([0, 2, 4]), np.array([0, 2, 4])]
    hist = planar_cartesian_histogram(positions, edges, dim=2)
    # counts are whole numbers, exact in float64:
    assert_array_equal(hist, np.array([[1.0, 0.0], [0.0, 2.0]]))


@pytest.mark.parametrize(
    "histogram,kwargs",
    [
        pytest.param(radial_histogram, {}, id="radial"),
        pytest.param(radial_cyl_histogram, {"dim": 1}, id="radial-cyl"),
        pytest.param(axial_histogram, {"dim": 0}, id="axial"),
        pytest.param(azimuth_cyl_histogram, {"dim": 2}, id="azimuth-cyl"),
        pytest.param(
            planar_cartesian_histogram, {"dim": 2}, id="planar-cartesian"
        ),
    ],
)
def test_histograms_reject_one_dimensional_positions(histogram, kwargs):
    """A one-dimensional positions array is rejected."""
    with pytest.raises(ValueError, match="must be a 2D array"):
        histogram(np.array([1.0, 2.0, 3.0]), np.array([0.0, 1.0]), **kwargs)


@pytest.mark.parametrize(
    "histogram",
    [
        pytest.param(radial_cyl_histogram, id="radial-cyl"),
        pytest.param(axial_histogram, id="axial"),
        pytest.param(azimuth_cyl_histogram, id="azimuth-cyl"),
        pytest.param(planar_cartesian_histogram, id="planar-cartesian"),
    ],
)
def test_histograms_reject_axis_outside_three_dimensions(
    positions_3d, histogram
):
    """An axis outside {0, 1, 2} is rejected."""
    with pytest.raises(ValueError, match="must be one of"):
        histogram(positions_3d, np.array([0.0, 1.0]), dim=5)


def test_planar_cartesian_histogram_rejects_other_than_two_edge_arrays(
    positions_3d,
):
    """The planar histogram needs exactly one edge array per planar axis."""
    edges = [np.array([0.0, 1.0])] * 3
    with pytest.raises(ValueError, match="must contain two arrays"):
        planar_cartesian_histogram(positions_3d, edges, dim=2)


@pytest.mark.parametrize(
    "histogram,kwargs",
    [
        pytest.param(
            radial_histogram, {"bin_range": (0, 3)}, id="radial"
        ),
        pytest.param(
            radial_cyl_histogram,
            {"bin_range": (0, 3), "dim": 2},
            id="radial-cyl",
        ),
        pytest.param(
            axial_histogram, {"bin_range": (0, 3), "dim": 0}, id="axial"
        ),
        pytest.param(
            azimuth_cyl_histogram,
            {"bin_range": (0, 3), "dim": 2},
            id="azimuth-cyl",
        ),
        pytest.param(
            planar_cartesian_histogram,
            {"bin_ranges": [(0, 3), (0, 3)], "dim": 2},
            id="planar-cartesian",
        ),
    ],
)
def test_histograms_reject_range_argument(positions_3d, histogram, kwargs):
    """The histogram functions no longer accept a range argument, which
    np.histogram ignored when given edges."""
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        histogram(positions_3d, edges, **kwargs)


@pytest.mark.parametrize(
    "histogram,kwargs",
    [
        pytest.param(radial_histogram, {}, id="radial"),
        pytest.param(radial_cyl_histogram, {"dim": 2}, id="radial-cyl"),
        pytest.param(axial_histogram, {"dim": 0}, id="axial"),
        pytest.param(azimuth_cyl_histogram, {"dim": 2}, id="azimuth-cyl"),
    ],
)
def test_histograms_reject_positions_beyond_edges(histogram, kwargs):
    """A position outside the outermost edges raises instead of being
    silently dropped; the message counts the positions that fell outside."""
    positions = np.array(
        [[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [2.5, 0.0, 0.0], [-9.0, 0.0, 0.0]]
    )
    # (-9, 0, 0) lies beyond the last edge for every histogram: its radial
    # and cylindrical-radial distances are 9, its x coordinate is -9, and its
    # azimuth about z is pi > 3.
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="1 of 4 positions"):
        histogram(positions, edges, **kwargs)


def test_planar_cartesian_histogram_rejects_positions_beyond_edges():
    """A position outside the planar edges raises instead of being silently
    dropped."""
    positions = np.array([[1.0, 1.0, 0.0], [3.0, 3.0, 0.0], [9.0, 1.0, 0.0]])
    edges = [np.array([0.0, 2.0, 4.0]), np.array([0.0, 2.0, 4.0])]
    with pytest.raises(ValueError, match="1 of 3 positions"):
        planar_cartesian_histogram(positions, edges, dim=2)


def test_radial_histogram_counts_positions_on_both_outer_edges():
    """Positions at exactly the first and the last edge are both counted, so
    nothing on the boundary of the range is lost."""
    positions = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    hist = radial_histogram(positions, edges)
    assert_array_equal(hist, np.array([1, 0, 1]))


# Bounds keep the bin count at most 1000 and every range far wider than the
# float64 spacing near |lmin| <= 100 (about 1e-14).
@given(
    bin_size=st.floats(min_value=0.01, max_value=10.0),
    lmin=st.floats(min_value=-100.0, max_value=100.0),
    n_bins_real=st.floats(min_value=0.01, max_value=1000.0),
    dim=st.sampled_from([0, 1, 2]),
    data=st.data(),
)
def test_axial_histogram_counts_every_position_in_range_once(
    bin_size: float,
    lmin: float,
    n_bins_real: float,
    dim: AxisT,
    data: st.DataObject,
) -> None:
    """Every position in [lmin, lmax] is counted exactly once over the edges
    create_bin_edge_and_hist builds for that range (a conservation
    property)."""
    lmax = lmin + n_bins_real * bin_size
    positions = data.draw(
        arrays(
            np.float64,
            st.tuples(st.integers(min_value=1, max_value=50), st.just(3)),
            elements=st.floats(min_value=lmin, max_value=lmax),
        )
    )
    edges, _ = create_bin_edge_and_hist(bin_size, lmin, lmax)
    hist = axial_histogram(positions, edges, dim=dim)
    assert hist.sum() == positions.shape[0]
