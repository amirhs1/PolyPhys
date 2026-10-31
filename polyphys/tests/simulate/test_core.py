"""Tests for :mod:`polyphys.simulate.core`.

Every expected value follows by hand from the box formulas of the LAMMPS
``Howto_triclinic`` page and ``Domain::lamda2x`` and ``Domain::unmap``: edge
vectors A = (lx, 0, 0), B = (xy, ly, 0), C = (xz, yz, lz); scaled
coordinates s map to o + s_a A + s_b B + s_c C; image flags i unwrap r to
r + i_x A + i_y B + i_z C.
"""

import numpy as np
import pytest

from polyphys.simulate.core import Box, LammpsFormatError
from polyphys.simulate.data import BrokenDataError
from polyphys.simulate.dump import BrokenDumpError
from polyphys.simulate.log import BrokenLogError


@pytest.fixture
def make_box():
    """Return a factory of boxes; keyword arguments replace the defaults.

    The default is the triclinic box with lo (1, 2, 3), hi (5, 7, 9), and
    tilt factors (0.5, -1.0, 2.0), so lengths (4, 5, 6).
    """
    def _make(**overrides):
        fields = {"lo": np.array([1.0, 2.0, 3.0]),
                  "hi": np.array([5.0, 7.0, 9.0]),
                  "tilt": np.array([0.5, -1.0, 2.0]), "triclinic": True}
        fields.update(overrides)
        return Box(**fields)
    return _make


def test_box_edges_follow_the_tilt_factors(make_box):
    """The edges are A = (lx, 0, 0), B = (xy, ly, 0), C = (xz, yz, lz)."""
    box = make_box()

    edges = box.edges()

    np.testing.assert_array_equal(
        edges, [[4.0, 0.0, 0.0], [0.5, 5.0, 0.0], [-1.0, 2.0, 6.0]])


def test_box_to_cartesian_adds_the_scaled_edges_to_the_origin(make_box):
    """Scaled (0.5, 0.2, 0.5) maps to o + 0.5 A + 0.2 B + 0.5 C."""
    box = make_box()

    cartesian = box.to_cartesian([[0.5, 0.2, 0.5]])

    # x = 1 + 2 + 0.1 - 0.5, y = 2 + 1 + 1, z = 3 + 3, by hand; the products
    # of these binary fractions are exact in float64.
    np.testing.assert_allclose(cartesian, [[2.6, 4.0, 6.0]], rtol=0,
                               atol=1e-15)


def test_box_unwrap_adds_whole_edges_per_image_flag(make_box):
    """Image flags (1, -1, 2) move r by A - B + 2 C."""
    box = make_box()

    unwrapped = box.unwrap([[0.0, 0.0, 0.0]], [[1, -1, 2]])

    # A - B + 2C = (4 - 0.5 - 2, -5 + 4, 12), by hand.
    np.testing.assert_array_equal(unwrapped, [[1.5, -1.0, 12.0]])


def test_box_tilt_defaults_to_zero():
    """An orthogonal box made from its bounds alone has zero tilt."""
    box = Box(lo=np.zeros(3), hi=np.ones(3))

    np.testing.assert_array_equal(box.tilt, [0.0, 0.0, 0.0])
    assert not box.triclinic
    assert box.boundary is None


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        pytest.param({"hi": np.ones(2)}, r"'hi' must have shape \(3,\)",
                     id="hi-shape"),
        pytest.param({"hi": np.array([5.0, 2.0, 9.0])}, "must be below",
                     id="hi-not-above-lo"),
        pytest.param({"triclinic": False}, "zero tilt factors",
                     id="orthogonal-with-tilt"),
    ],
)
def test_box_rejects_an_invalid_box(make_box, overrides, match):
    """Misshapen or misordered bounds, or tilt without triclinic, raise."""
    with pytest.raises(ValueError, match=match):
        make_box(**overrides)


def test_box_arrays_are_read_only(make_box):
    """A box's bounds cannot be changed in place."""
    box = make_box()

    with pytest.raises(ValueError, match="read-only"):
        box.lo[0] = 0.0


@pytest.mark.parametrize(
    "error", [BrokenLogError, BrokenDumpError, BrokenDataError])
def test_broken_file_errors_share_one_base(error):
    """Each broken-file error is a LammpsFormatError and a ValueError."""
    assert issubclass(error, LammpsFormatError)
    assert issubclass(error, ValueError)
