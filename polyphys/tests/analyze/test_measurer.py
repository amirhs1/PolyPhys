import warnings

import pytest
import numpy as np
from numpy.testing import assert_array_equal

from polyphys.analyze.measurer import (
    apply_pbc_orthogonal,
    pair_distance,
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


@pytest.mark.parametrize(
    "pbc,expected_lengths,expected_inv",
    [
        ({0: 10.0, 2: 20.0}, [10.0, 0.0, 20.0], [0.1, 0.0, 0.05]),
    ],
)
def test_apply_pbc_orthogonal(pbc, expected_lengths, expected_inv):
    lengths = np.zeros(3)
    inv = np.zeros(3)
    updated, updated_inv = apply_pbc_orthogonal(lengths, inv, pbc)
    assert_array_equal(updated, expected_lengths)
    assert_array_equal(updated_inv, expected_inv)


def test_apply_pbc_zero_length_raises():
    with pytest.raises(ZeroDivisionError):
        apply_pbc_orthogonal(np.zeros(3), np.zeros(3), {1: 0.0})


def test_apply_pbc_negative_length_raises():
    with pytest.raises(ValueError):
        apply_pbc_orthogonal(np.zeros(3), np.zeros(3), {1: -5.0})


def test_pair_distance_without_pbc():
    pos = np.array([[1.0, 2.0, 3.0], [4.0, 2.0, 3.0]])
    result = pair_distance(pos)
    assert_array_equal(result, np.array([3.0, 0.0, 0.0]))


def test_pair_distance_invalid_shape_raises():
    with pytest.raises(ValueError):
        pair_distance(np.array([[1.0, 2.0, 3.0]]))  # only one position


def test_transverse_size():
    pos = np.array([[1, 0, 0], [1, 2, 2]])
    result = transverse_size(pos, axis=0)
    assert np.isclose(result, 2.82842712)


def test_transverse_size_invalid_axis():
    pos = np.array([[1, 0, 0], [1, 2, 2]])
    with pytest.raises(IndexError):
        transverse_size(pos, axis=3)


def test_max_distance():
    pos = np.array([[0, 0, 0], [4, 2, 6]])
    result = max_distance(pos)
    assert_array_equal(result, np.array([4.0, 2.0, 6.0]))


def test_fsd():
    pos = np.array([[1, 2, 3], [4, 8, 6]])
    result = fsd(pos, axis=1)
    assert result == 6.0


def test_fsd_invalid_axis():
    pos = np.array([[1, 2, 3]])
    with pytest.raises(IndexError):
        fsd(pos, axis=5)


def test_simple_stats_valid():
    data = np.array([1.0, 2.0, 3.0])
    result = simple_stats("energy", data)
    assert result["energy_mean"] == 2.0
    assert result["energy_var"] == 1.0
    assert result["energy_sem"] == pytest.approx(0.57735, rel=1e-5)


def test_simple_stats_empty():
    with pytest.raises(ValueError, match="must not be empty"):
        simple_stats("test", np.array([]))


def test_sem_valid():
    result = sem(np.array([1.0, 2.0, 3.0]))
    assert result == pytest.approx(0.57735, rel=1e-5)


def test_sem_empty():
    with pytest.raises(ValueError, match="must not be empty"):
        sem(np.array([]))


def test_spherical_segment_basic():
    assert spherical_segment(3, 1, 2) == \
        pytest.approx(20.94395102393195, rel=1e-5)


def test_spherical_segment_full_sphere():
    assert spherical_segment(3, -3, 3) == \
        pytest.approx(113.09733552923255, rel=1e-5)


def test_spherical_segment_negative_radius():
    with pytest.raises(ValueError, match="must be positive"):
        spherical_segment(-3, 1, 2)


def test_sphere_sphere_intersection_partial_overlap():
    assert sphere_sphere_intersection(3, 4, 2) == \
        pytest.approx(94.90227807719167, rel=1e-5)


def test_sphere_sphere_intersection_no_overlap():
    assert sphere_sphere_intersection(3, 4, 10) == 0.0


def test_sphere_sphere_intersection_fully_contained():
    assert sphere_sphere_intersection(3, 4, 1) == \
        pytest.approx(113.0973, rel=1e-5)


def test_sphere_sphere_intersection_zero_radius():
    assert sphere_sphere_intersection(3, 0, 1) == 0.0
    assert sphere_sphere_intersection(0, 4, 1) == 0.0


def test_create_bin_edge_and_hist():
    edges, hist = create_bin_edge_and_hist(1.0, 0.0, 5.0)
    assert_array_equal(edges, np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]))
    assert_array_equal(hist, np.zeros(5, dtype=np.int16))


def test_create_bin_edge_and_hist_invalid_range():
    with pytest.raises(ValueError):
        create_bin_edge_and_hist(1.0, 5.0, 0.0)


def test_radial_histogram():
    """Distances from the origin are binned; the last edge is inclusive."""
    pos = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    edges = np.array([0, 1, 2, 3])
    result = radial_histogram(pos, edges)
    assert_array_equal(result, np.array([0, 2, 1]))


def test_radial_histogram_invalid_input():
    """A one-dimensional positions array is rejected."""
    with pytest.raises(ValueError, match="must be a 2D array"):
        radial_histogram(np.array([1, 2, 3]), np.array([0, 1]))


@pytest.mark.parametrize(
    "bin_type,expected_n_bins",
    [("ordinary", 5), ("nonnegative", 5), ("periodic", 5)],
)
def test_fixedsize_bins(bin_type, expected_n_bins):
    result = fixedsize_bins(1.0, 0.0, 5.0, bin_type=bin_type)
    assert result["n_bins"] == expected_n_bins
    assert len(result["bin_edges"]) == expected_n_bins + 1


def test_histogram_collectors_do_not_overflow_int16():
    """Collectors accumulate counts past the int16 limit without wrapping."""
    n_per_frame = 30000
    frame = np.full(n_per_frame, 0.5)  # every sample falls in bin [0, 1)
    info = fixedsize_bins(1.0, 0.0, 2.0)
    _, hist = create_bin_edge_and_hist(1.0, 0.0, 2.0)
    assert info["collector"].dtype == np.int64
    assert info["collector_std"].dtype == np.float64
    assert hist.dtype == np.int64
    for _ in range(2):
        counts, _ = np.histogram(frame, bins=info["bin_edges"])
        info["collector"] += counts
        info["collector_std"] += counts**2
        hist += np.histogram(frame, bins=np.array([0.0, 1.0, 2.0]))[0]
    expected = np.array([2 * n_per_frame, 0])
    assert_array_equal(info["collector"], expected)
    assert_array_equal(hist, expected)
    assert_array_equal(
        info["collector_std"], np.array([2.0 * n_per_frame**2, 0.0])
    )


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


@pytest.fixture
def positions_3d():
    """Three atom positions in 3D, for tests where their values are
    incidental."""
    return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]])


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


def test_radial_cyl_histogram_valid():
    """Distances from the y axis are binned in the x-z plane."""
    positions = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3], [0, 3, 4]])
    edges = np.array([0, 1, 2, 3, 4])
    # distances from the y axis: 1, 0, 3, 4; the last edge is inclusive.
    hist = radial_cyl_histogram(positions, edges, dim=1)
    assert_array_equal(hist, np.array([1, 1, 0, 2]))


def test_radial_cyl_histogram_invalid_dim(positions_3d):
    """An axis outside {0, 1, 2} is rejected."""
    with pytest.raises(ValueError, match="must be one of"):
        radial_cyl_histogram(positions_3d, np.array([0, 1]), dim=5)


def test_radial_cyl_histogram_invalid_ndim():
    """A one-dimensional positions array is rejected."""
    pos = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="must be a 2D array"):
        radial_cyl_histogram(pos, np.array([0, 1]), dim=1)


def test_axial_histogram_valid():
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


def test_axial_histogram_invalid_dim(positions_3d):
    """An axis outside {0, 1, 2} is rejected."""
    with pytest.raises(ValueError, match="must be one of"):
        axial_histogram(positions_3d, np.array([0, 1]), dim=4)


def test_axial_histogram_invalid_ndim():
    """A one-dimensional positions array is rejected."""
    pos = np.array([1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="must be a 2D array"):
        axial_histogram(pos, np.array([0, 1]), dim=0)


def test_azimuth_cyl_histogram():
    """Azimuth angles about the z axis are binned over [-pi, pi]."""
    pos = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0]])  # x-y plane
    edges = np.linspace(-np.pi, np.pi, 5)
    # angles 0, pi, pi/2; bins [-pi, -pi/2), [-pi/2, 0), [0, pi/2),
    # [pi/2, pi], the last edge inclusive.
    result = azimuth_cyl_histogram(pos, edges, dim=2)
    assert_array_equal(result, np.array([0, 0, 1, 2]))


def test_planar_cartesian_histogram_normal_case():
    """Positions are binned in the x-y plane when the axis is z."""
    positions = np.array([
        [1, 1, 0],  # falls into bin (0,0)
        [2, 2, 0],  # falls into bin (1,1)
        [3, 3, 0]   # falls into bin (1,1)
    ])
    edges = [np.array([0, 2, 4]), np.array([0, 2, 4])]

    # Project onto plane perpendicular to z-axis (dim=2), i.e. use x and y
    hist = planar_cartesian_histogram(positions, edges, dim=2)

    expected = np.array([
        [1.0, 0.0],
        [0.0, 2.0]
    ])
    np.testing.assert_array_equal(hist, expected)
