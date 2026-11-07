"""Tests for :mod:`polyphys.viz.helpers`.

Every expected tick value is hand arithmetic on the module's rule: ticks
every `step` from `lower` to `upper`, limits widened by ``step / 10``
rounded down to its first non-zero digit (``0.025 -> 0.02``,
``0.05 -> 0.05``), labels rounded to `decimals`. Every expected legend value
is the input the handle or legend was built from: the labels, colors,
styles, and titles given, read back from the artists.
"""
import pytest

pytest.importorskip("matplotlib")

import numpy as np  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from polyphys.viz.helpers import (  # noqa: E402
    band_under_lines,
    color_handles,
    line_handles,
    move_legend_to_axes,
    set_xticks,
    set_yticks,
)


@pytest.fixture
def two_panels(figure):
    """Two side-by-side axes; the left one holds two labeled lines and the
    figure a legend of them titled 'curves'."""
    left, right = figure.subplots(ncols=2)
    left.plot([0, 1], [0, 1], label="rise")
    left.plot([0, 1], [1, 0], label="fall")
    legend = figure.legend(title="curves")
    return left, right, legend


def test_set_xticks_places_ticks_every_step_to_upper(ax):
    """Ticks run from lower to upper in steps of step, upper included."""
    ticks = set_xticks(ax, 0.0, 1.0, 0.25)

    # Exact: the positions are lower + k * step for a step with a short
    # binary expansion.
    np.testing.assert_array_equal(ticks, [0.0, 0.25, 0.5, 0.75, 1.0])
    np.testing.assert_array_equal(ax.get_xticks(), ticks)


def test_set_xticks_widens_the_limits_by_the_tolerance(ax):
    """The limits sit one rounded-down tenth of a step beyond the ticks."""
    set_xticks(ax, 0.0, 1.0, 0.25)

    # 0.25 / 10 = 0.025 rounds down to 0.02; float64 round-off tolerance.
    assert ax.get_xlim() == pytest.approx((-0.02, 1.02), abs=1e-12)


def test_set_xticks_writes_a_near_zero_tick_as_exactly_zero(ax):
    """A tick that float arithmetic puts near zero is exactly 0.0."""
    # -0.3 + 3 * 0.1 is 5.55e-17 in float64, not 0.0.
    ticks = set_xticks(ax, -0.3, 0.3, 0.1)

    assert ticks[3] == 0.0
    assert ax.get_xticklabels()[3].get_text() == "0.0"


def test_set_xticks_hides_labels_on_request(ax):
    """labels=False leaves every tick label empty."""
    set_xticks(ax, 0.0, 1.0, 0.25, labels=False)

    assert [label.get_text() for label in ax.get_xticklabels()] == [""] * 5


def test_set_xticks_rounds_labels_to_decimals(ax):
    """Labels are the ticks rounded to `decimals` places."""
    set_xticks(ax, 0.0, 0.9, 0.3, decimals=1)

    assert ([label.get_text() for label in ax.get_xticklabels()]
            == ["0.0", "0.3", "0.6", "0.9"])


def test_set_xticks_forwards_text_kwargs_to_the_labels(ax):
    """Keyword arguments such as fontsize reach every tick label."""
    set_xticks(ax, 0.0, 1.0, 0.5, fontsize=7)

    assert all(label.get_fontsize() == 7 for label in ax.get_xticklabels())


@pytest.mark.parametrize("lower, upper, step", [
    pytest.param(0.0, 1.0, 0.0, id="zero-step"),
    pytest.param(0.0, 1.0, -0.1, id="negative-step"),
    pytest.param(0.0, 1.0, 1.0, id="step-equals-span"),
    pytest.param(1.0, 0.0, 0.1, id="lower-above-upper"),
])
def test_set_xticks_rejects_a_bad_step_or_order(ax, lower, upper, step):
    """A step that is not positive and below the span, or lower >= upper,
    is refused."""
    with pytest.raises(ValueError, match="step must|lower must"):
        set_xticks(ax, lower, upper, step)


def test_set_yticks_applies_the_same_rule_on_the_y_axis(ax):
    """set_yticks places the ticks and widens the limits of the y axis."""
    ticks = set_yticks(ax, 0.0, 1.0, 0.5)

    np.testing.assert_array_equal(ticks, [0.0, 0.5, 1.0])
    np.testing.assert_array_equal(ax.get_yticks(), ticks)
    # 0.5 / 10 = 0.05 is already one significant digit.
    assert ax.get_ylim() == pytest.approx((-0.05, 1.05), abs=1e-12)


def test_color_handles_make_one_patch_per_label_in_order():
    """Each label gets a Patch of its color, in the order given."""
    handles = color_handles(["a", "b"], ["#ff0000", "#0000ff"])

    assert all(isinstance(handle, Patch) for handle in handles)
    assert [handle.get_label() for handle in handles] == ["a", "b"]
    assert [handle.get_facecolor() for handle in handles] == [
        to_rgba("#ff0000"), to_rgba("#0000ff")]


@pytest.mark.parametrize("labels, color", [
    pytest.param(["a", "b", "c"], (0.1, 0.2, 0.3), id="rgb-tuple"),
    pytest.param(["a", "b"], (0.1, 0.2, 0.3, 0.5), id="rgba-tuple"),
    pytest.param(["a", "b"], "red", id="name"),
])
def test_color_handles_repeat_one_color_for_every_label(labels, color):
    """One color, a name or an RGB(A) tuple, is given to every label."""
    handles = color_handles(labels, color)

    assert [handle.get_facecolor() for handle in handles] == \
        [to_rgba(color)] * len(labels)


def test_color_handles_reject_a_list_of_the_wrong_length():
    """A list of colors must have one entry per label."""
    with pytest.raises(ValueError, match="colors has 3 entries"):
        color_handles(["a", "b"], ["red", "green", "blue"])


def test_line_handles_repeat_a_numpy_integer_width():
    """A NumPy integer width is one value, given to every label."""
    handles = line_handles(["a", "b"], linewidths=np.int64(2))

    assert [h.get_linewidth() for h in handles] == [2, 2]


def test_line_handles_repeat_a_dash_tuple_for_every_label():
    """An (offset, on-off sequence) dash tuple is one line style."""
    handles = line_handles(["a", "b"], linestyles=(0, (5, 2)))

    assert [h.get_linestyle() for h in handles] == ["--", "--"]


def test_line_handles_carry_the_style_marker_and_color_given():
    """Each label gets a Line2D with its own style and marker and the shared
    color and width."""
    handles = line_handles(["a", "b"], linestyles=["-", "--"],
                           markers=["o", "^"], linewidths=2, color="red")

    assert all(isinstance(handle, Line2D) for handle in handles)
    assert [(h.get_label(), h.get_linestyle(), h.get_marker())
            for h in handles] == [("a", "-", "o"), ("b", "--", "^")]
    assert [h.get_linewidth() for h in handles] == [2, 2]
    assert [h.get_color() for h in handles] == ["red", "red"]


def test_line_handles_reject_a_sequence_of_the_wrong_length():
    """A per-label sequence must have one entry per label."""
    with pytest.raises(ValueError, match="markers has 3 entries"):
        line_handles(["a", "b"], markers=["o", "^", "s"])


def test_move_legend_to_axes_moves_entries_and_title_into_the_axes(
        two_panels):
    """The figure legend disappears and each target axes gets a legend with
    the same labels and title."""
    left, right, legend = two_panels

    new = move_legend_to_axes(legend, [left, right],
                              ["upper left", "lower right"])

    assert left.figure.legends == []
    assert [leg.axes for leg in new] == [left, right]
    for panel in (left, right):
        texts = [t.get_text() for t in panel.get_legend().get_texts()]
        assert texts == ["rise", "fall"]
        assert panel.get_legend().get_title().get_text() == "curves"


def test_move_legend_to_axes_takes_a_new_title(two_panels):
    """A title argument replaces the moved legend's title."""
    left, right, legend = two_panels

    move_legend_to_axes(legend, [right], ["upper left"], title="new")

    assert right.get_legend().get_title().get_text() == "new"


def test_move_legend_to_axes_rejects_mismatched_locations(two_panels):
    """One location per target axes is required."""
    left, right, legend = two_panels

    with pytest.raises(ValueError, match="locs has 1 entries for 2 axes"):
        move_legend_to_axes(legend, [left, right], ["upper left"])


def test_band_under_lines_shades_each_band_in_its_line_color(ax):
    """One band per line is filled between the given edges, over the line's
    x data, in the line's color."""
    ax.plot([0, 1, 2], [0.5, 0.5, 0.5], color="#ff0000")
    ax.plot([0, 1, 2], [1.5, 1.5, 1.5], color="#0000ff")

    bands = band_under_lines(
        ax, [([0, 0, 0], [1, 1, 1]), ([1, 1, 1], [2, 2, 2])], alpha=0.5)

    assert len(bands) == 2
    assert [tuple(b.get_facecolor()[0]) for b in bands] == [
        to_rgba("#ff0000", 0.5), to_rgba("#0000ff", 0.5)]
    for band, (lower, upper) in zip(bands, [(0, 1), (1, 2)]):
        vertices = band.get_paths()[0].vertices
        # Exact: the polygon vertices are the integer edges given.
        assert vertices[:, 1].min() == lower
        assert vertices[:, 1].max() == upper
        np.testing.assert_array_equal(np.unique(vertices[:, 0]), [0, 1, 2])


@pytest.mark.parametrize("bands, message", [
    pytest.param([([0], [1]), ([0], [1])], "2 bands were given for 1 lines",
                 id="more-bands-than-lines"),
    pytest.param([([0, 0], [1, 1])], "a band has edges of shape",
                 id="edge-length-differs-from-line"),
])
def test_band_under_lines_rejects_a_band_that_does_not_fit(ax, bands,
                                                           message):
    """More bands than lines, or an edge of the wrong length, is refused."""
    ax.plot([0, 1, 2], [0, 0, 0])

    with pytest.raises(ValueError, match=message):
        band_under_lines(ax, bands)
