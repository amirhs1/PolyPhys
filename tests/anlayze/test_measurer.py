import pytest
import numpy as np
from numpy.testing import assert_array_equal

from polyphys.analyze.measurer import (
    apply_pbc_orthogonal,
    pair_distance,
    transverse_size,
    max_distance,
    fsd,
    create_bin_edge_and_hist,
    fixedsize_bins,
    radial_histogram,
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


def test_apply_pbc_orthogonal():
    lengths = np.zeros(3)
    inv = np.zeros(3)
    updated, updated_inv = apply_pbc_orthogonal(
        lengths, inv, {0: 10.0, 2: 20.0}
    )
    assert_array_equal(updated, np.array([10.0, 0.0, 20.0]))
    assert_array_equal(updated_inv, np.array([0.1, 0.0, 0.05]))


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


def test_create_bin_edge_and_hist():
    edges, hist = create_bin_edge_and_hist(1.0, 0.0, 5.0)
    assert_array_equal(edges, np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]))
    assert_array_equal(hist, np.zeros(5, dtype=np.int16))


def test_create_bin_edge_and_hist_invalid_range():
    with pytest.raises(ValueError):
        create_bin_edge_and_hist(1.0, 5.0, 0.0)


def test_radial_histogram():
    pos = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    edges = np.array([0, 1, 2, 3])
    result = radial_histogram(pos, edges, (0, 3))
    assert_array_equal(result, np.array([0, 2, 1]))


def test_radial_histogram_invalid_input():
    with pytest.raises(ValueError):
        radial_histogram(np.array([1, 2, 3]), np.array([0, 1]), (0, 1))


@pytest.mark.parametrize(
    "bin_type,expected_n_bins",
    [("ordinary", 5), ("nonnegative", 5), ("periodic", 5)],
)
def test_fixedsize_bins(bin_type, expected_n_bins):
    result = fixedsize_bins(1.0, 0.0, 5.0, bin_type=bin_type)
    assert result["n_bins"] == expected_n_bins
    assert len(result["bin_edges"]) == expected_n_bins + 1


def test_azimuth_cyl_histogram():
    from polyphys.analyze.measurer import azimuth_cyl_histogram

    pos = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0]])  # x-y plane
    edges = np.linspace(-np.pi, np.pi, 5)
    result = azimuth_cyl_histogram(pos, edges, (-np.pi, np.pi), dim=2)
    assert result.sum() == 3
