"""Tests for :mod:`polyphys.viz.colors`.

The expected colors of the three Tol palettes are read from Figures 1, 3,
and 4 of Paul Tol, *Colour Schemes*, SRON/EPS/TN/09-002, issue 3.2
(18 August 2021), in the fixed-sequence order its page 6 gives; those of
``DUTCH_FIELD`` from the wpDataTables article the module docstring cites,
read on 2026-11-06. The form test states a property of every palette. The
expected endpoints of a truncated map are the parent map's own values at
the cut points, a reproducible computation on the parent; the count and
name follow from the arguments.
"""
import re

import pytest

from polyphys.viz.colors import (
    DUTCH_FIELD,
    STUDY_COLORS,
    TOL_BRIGHT,
    TOL_MUTED,
    TOL_VIBRANT,
    truncate_colormap,
)


@pytest.fixture
def viridis():
    """matplotlib's registered viridis map, the parent of every cut."""
    matplotlib = pytest.importorskip("matplotlib")
    return matplotlib.colormaps["viridis"]


@pytest.mark.parametrize("palette, expected", [
    pytest.param(
        TOL_BRIGHT,
        ("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377",
         "#BBBBBB"),
        id="tol-bright-fig1"),
    pytest.param(
        TOL_VIBRANT,
        ("#EE7733", "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988",
         "#BBBBBB"),
        id="tol-vibrant-fig3"),
    pytest.param(
        TOL_MUTED,
        ("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255",
         "#44AA99", "#999933", "#AA4499", "#DDDDDD"),
        id="tol-muted-fig4"),
    pytest.param(
        DUTCH_FIELD,
        ("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300",
         "#dc0ab4", "#b3d4ff", "#00bfa0"),
        id="dutch-field-wpdatatables"),
])
def test_palette_matches_its_source(palette, expected):
    """Each cited palette holds its source's colors in the source's order."""
    assert palette == expected


@pytest.mark.parametrize("palette", [
    pytest.param(TOL_BRIGHT, id="tol-bright"),
    pytest.param(TOL_VIBRANT, id="tol-vibrant"),
    pytest.param(TOL_MUTED, id="tol-muted"),
    pytest.param(DUTCH_FIELD, id="dutch-field"),
    pytest.param(STUDY_COLORS, id="study-colors"),
])
def test_palette_entries_are_distinct_hex_colors(palette):
    """Every palette entry is a #RRGGBB string, and no color repeats."""
    hex_color = re.compile(r"^#[0-9A-Fa-f]{6}$")

    assert all(hex_color.match(color) for color in palette)
    assert len(set(color.lower() for color in palette)) == len(palette)


def test_truncate_colormap_endpoints_are_the_parent_values(viridis):
    """The cut map starts at cmap(lower) and ends at cmap(upper)."""
    np = pytest.importorskip("numpy")

    cut = truncate_colormap(viridis, 0.2, 0.8, n=5)

    # Exact: the endpoint colors are copied from the parent's lookup table.
    np.testing.assert_array_equal(cut(0.0), viridis(0.2))
    np.testing.assert_array_equal(cut(1.0), viridis(0.8))


def test_truncate_colormap_has_n_colors_and_a_named_range(viridis):
    """The cut map samples n colors and is named trunc(<name>,lower,upper)."""
    cut = truncate_colormap(viridis, 0.2, 0.8, n=5)

    assert cut.N == 5
    assert cut.name == "trunc(viridis,0.2,0.8)"


def test_truncate_colormap_accepts_a_registered_name(viridis):
    """A color map name is looked up in matplotlib's registry."""
    np = pytest.importorskip("numpy")

    cut = truncate_colormap("viridis", 0.2, 0.8, n=5)

    np.testing.assert_array_equal(cut(0.0), viridis(0.2))


@pytest.mark.parametrize("lower, upper, n", [
    pytest.param(0.8, 0.2, 5, id="lower-above-upper"),
    pytest.param(-0.1, 1.0, 5, id="lower-below-zero"),
    pytest.param(0.0, 1.1, 5, id="upper-above-one"),
    pytest.param(0.0, 1.0, 1, id="one-color"),
])
def test_truncate_colormap_rejects_a_bad_range_or_count(viridis, lower,
                                                        upper, n):
    """A range outside 0 <= lower < upper <= 1, or n below 2, is refused."""
    with pytest.raises(ValueError, match="lower|n must"):
        truncate_colormap(viridis, lower, upper, n=n)
