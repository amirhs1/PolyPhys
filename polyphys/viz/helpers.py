"""
==========================================================
:mod:`polyphys.viz.helpers`
==========================================================

The :mod:`polyphys.viz.helpers` module holds the matplotlib helpers the
studies' multi-panel figures use: a rule for placing an axis's ticks, proxy
legend handles, a legend moved into chosen panels, and a band shaded under
each line. It needs matplotlib, which the ``viz`` extra installs, and raises
:exc:`ImportError` with that hint when it is missing.

Ticks
=====

:func:`set_xticks` and :func:`set_yticks` place the major ticks of one axis
by one rule. Given `lower`, `upper`, and `step`:

- the ticks are :math:`t_k = \\mathrm{lower} + k \\, \\mathrm{step}` for
  :math:`k = 0, 1, \\dots, \\lfloor (\\mathrm{upper} - \\mathrm{lower}) /
  \\mathrm{step} \\rfloor`, so the last tick is `upper` when `step` divides
  the span;
- the *tolerance* is ``step / 10`` rounded down to its first non-zero
  digit, by :func:`polyphys.manage.utils.round_down_first_non_zero`; a tick
  closer to zero than the tolerance is written as exactly ``0.0``, so
  ``-5.6e-17`` never appears as a label;
- the axis limits are ``(lower - tolerance, upper + tolerance)``, a frame
  that keeps the first and last ticks off the axes' edges;
- the labels are the ticks rounded to `decimals` places, or hidden.

The ticks are not rounded to multiples of a power of ten; the caller
chooses `lower`, `upper`, and `step` with that in mind.

Legends
=======

- :func:`color_handles` and :func:`line_handles` make *proxy* handles,
  artists that draw nothing on the axes and exist only to sit in a legend,
  so a figure can carry one legend for its colors and another for its line
  styles or markers, whatever the plotting call produced.
- :func:`move_legend_to_axes` takes an existing legend, such as the figure
  legend a facet grid draws outside its panels, and re-creates it inside
  chosen panels.
- :func:`band_under_lines` shades a band, such as a confidence interval,
  under each line already drawn on an axes, in that line's color.

Examples
========

>>> import matplotlib
>>> matplotlib.use("Agg")
>>> import matplotlib.pyplot as plt
>>> from polyphys.viz.helpers import color_handles, set_xticks
>>> fig, ax = plt.subplots()
>>> set_xticks(ax, 0.0, 1.0, 0.25)
array([0.  , 0.25, 0.5 , 0.75, 1.  ])
>>> [float(limit) for limit in ax.get_xlim()]
[-0.02, 1.02]
>>> handles = color_handles(["bug", "crowder"], ["#4477AA", "#EE6677"])
>>> [text.get_text() for text in ax.legend(handles=handles).get_texts()]
['bug', 'crowder']
>>> plt.close(fig)
"""
import numbers
from typing import Any, Literal, Optional, Sequence, Union

import numpy as np

from ..manage.utils import round_down_first_non_zero

try:
    import matplotlib  # noqa: F401
except ImportError as error:  # pragma: no cover - depends on the install
    raise ImportError(
        "polyphys.viz.helpers needs matplotlib; install it with "
        "python -m pip install 'polyphys[viz]'.") from error

from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

__all__ = ["band_under_lines", "color_handles", "line_handles",
           "move_legend_to_axes", "set_xticks", "set_yticks"]


def _ticks(lower: float, upper: float, step: float) -> tuple[Any, float]:
    """Return the tick positions and the tolerance of the ticks rule."""
    if not lower < upper:
        raise ValueError(
            f"lower must be below upper, got lower={lower} and "
            f"upper={upper}.")
    span = upper - lower
    if step <= 0 or step >= span:
        raise ValueError(
            f"step must be positive and smaller than the span "
            f"upper - lower = {span}, got step={step}.")
    tolerance = round_down_first_non_zero(step / 10)
    # The 1e-9 absorbs the float error of a span that step divides exactly.
    count = int(np.floor(span / step + 1e-9)) + 1
    ticks = lower + step * np.arange(count)
    ticks = np.where(np.abs(ticks) < tolerance, 0.0, ticks)
    return ticks, tolerance


def _set_ticks(
    ax: Axes,
    which: Literal["x", "y"],
    lower: float,
    upper: float,
    step: float,
    labels: bool,
    decimals: int,
    text_kwargs: dict[str, Any],
) -> Any:
    """Apply the ticks rule to the x or y axis of `ax`."""
    ticks, tolerance = _ticks(lower, upper, step)
    texts = ([str(v) for v in np.around(ticks, decimals)] if labels
             else [""] * len(ticks))
    if which == "x":
        ax.set_xlim(lower - tolerance, upper + tolerance)
        ax.set_xticks(ticks)
        ax.set_xticklabels(texts, **text_kwargs)
    else:
        ax.set_ylim(lower - tolerance, upper + tolerance)
        ax.set_yticks(ticks)
        ax.set_yticklabels(texts, **text_kwargs)
    return ticks


def set_xticks(
    ax: Axes,
    lower: float,
    upper: float,
    step: float,
    *,
    labels: bool = True,
    decimals: int = 3,
    **text_kwargs: Any,
) -> Any:
    """
    Place the major x ticks of `ax` every `step` from `lower` to `upper`.

    The ticks, the tolerance that frames them, and the labels follow the
    rule in the module docstring. Time and space are :math:`O(n)` in the
    number of ticks.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes whose x axis is set.
    lower, upper : float
        The first tick and the largest tick value, in data units;
        ``lower < upper``.
    step : float
        The distance between ticks, positive and smaller than
        ``upper - lower``.
    labels : bool, default True
        Whether the ticks are labeled; ``False`` hides every label.
    decimals : int, default 3
        The number of decimals the labels are rounded to.
    **text_kwargs
        Passed to :meth:`matplotlib.axes.Axes.set_xticklabels`, such as
        ``fontsize``.

    Returns
    -------
    numpy.ndarray
        The tick positions, shape ``(n_ticks,)``.

    Raises
    ------
    ValueError
        If ``lower >= upper``, or `step` is not positive and smaller than
        ``upper - lower``.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> set_xticks(ax, -0.2, 0.2, 0.1, decimals=1)
    array([-0.2, -0.1,  0. ,  0.1,  0.2])
    >>> [label.get_text() for label in ax.get_xticklabels()]
    ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    >>> plt.close(fig)
    """
    return _set_ticks(ax, "x", lower, upper, step, labels, decimals,
                      text_kwargs)


def set_yticks(
    ax: Axes,
    lower: float,
    upper: float,
    step: float,
    *,
    labels: bool = True,
    decimals: int = 3,
    **text_kwargs: Any,
) -> Any:
    """
    Place the major y ticks of `ax` every `step` from `lower` to `upper`.

    The y-axis twin of :func:`set_xticks`, with the same rule, parameters,
    return value, and errors; `text_kwargs` go to
    :meth:`matplotlib.axes.Axes.set_yticklabels`.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes whose y axis is set.
    lower, upper : float
        The first tick and the largest tick value; ``lower < upper``.
    step : float
        The distance between ticks, positive and smaller than
        ``upper - lower``.
    labels : bool, default True
        Whether the ticks are labeled.
    decimals : int, default 3
        The number of decimals the labels are rounded to.
    **text_kwargs
        Passed to :meth:`matplotlib.axes.Axes.set_yticklabels`.

    Returns
    -------
    numpy.ndarray
        The tick positions, shape ``(n_ticks,)``.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> set_yticks(ax, 0.0, 1.0, 0.5, labels=False)
    array([0. , 0.5, 1. ])
    >>> [float(limit) for limit in ax.get_ylim()]
    [-0.05, 1.05]
    >>> plt.close(fig)
    """
    return _set_ticks(ax, "y", lower, upper, step, labels, decimals,
                      text_kwargs)


def _per_label(
    name: str,
    values: Union[None, str, float, Sequence[Any]],
    n_labels: int,
) -> list[Any]:
    """Return one value per label: one value repeated, or a matching list.

    One value is None, a string, a real number (NumPy scalars included), a
    tuple of real numbers such as an RGB or RGBA color, or a dash tuple
    ``(offset, (on, off, ...))``. Anything else is a sequence with one entry
    per label.
    """
    if _is_one_value(values):
        return [values] * n_labels
    values = list(values)
    if len(values) != n_labels:
        raise ValueError(
            f"{name} has {len(values)} entries for {n_labels} labels.")
    return values


def _is_one_value(value: Any) -> bool:
    """Return True if `value` is one style, not one per label."""
    if value is None or isinstance(value, (str, numbers.Real)):
        return True
    if not isinstance(value, tuple):
        return False
    if all(isinstance(item, numbers.Real) for item in value):
        return True
    return (len(value) == 2 and isinstance(value[0], numbers.Real)
            and isinstance(value[1], (tuple, list))
            and all(isinstance(item, numbers.Real) for item in value[1]))


def color_handles(
    labels: Sequence[str],
    colors: Any,
    **patch_kwargs: Any,
) -> list[Patch]:
    """
    Make one filled :class:`~matplotlib.patches.Patch` per label.

    The patches are proxy handles for a legend keyed by color. Time and
    space are :math:`O(n)` in the number of labels.

    Parameters
    ----------
    labels : sequence of str
        The legend entries.
    colors : matplotlib color or sequence of matplotlib colors
        One color for every label, a name or a tuple of numbers such as the
        RGB triple ``(0.1, 0.2, 0.3)``; or a list, or another sequence, with
        one color per label, in the same order.
    **patch_kwargs
        Passed to every :class:`~matplotlib.patches.Patch`, such as
        ``alpha``.

    Returns
    -------
    list of matplotlib.patches.Patch
        The handles, one per label, for the ``handles`` argument of
        :meth:`matplotlib.axes.Axes.legend`.

    Raises
    ------
    ValueError
        If `colors` is a sequence of colors whose length differs from that
        of `labels`.

    Examples
    --------
    >>> patches = color_handles(["low", "high"], ["#4477AA", "#EE6677"],
    ...                         alpha=0.5)
    >>> [patch.get_label() for patch in patches]
    ['low', 'high']
    >>> patches[0].get_alpha()
    0.5

    One RGB tuple colors every label:

    >>> patches = color_handles(["a", "b"], (0.0, 0.5, 1.0))
    >>> [patch.get_facecolor()[:3] for patch in patches]
    [(0.0, 0.5, 1.0), (0.0, 0.5, 1.0)]
    """
    colors = _per_label("colors", colors, len(labels))
    return [Patch(label=label, color=color, **patch_kwargs)
            for label, color in zip(labels, colors)]


def line_handles(
    labels: Sequence[str],
    *,
    linestyles: Union[None, str, Sequence[Any]] = None,
    linewidths: Union[None, float, Sequence[float]] = None,
    markers: Union[None, str, Sequence[Any]] = None,
    color: Any = "black",
    **line_kwargs: Any,
) -> list[Line2D]:
    """
    Make one empty :class:`~matplotlib.lines.Line2D` per label.

    The lines are proxy handles for a legend keyed by line style, width,
    or marker. Each of `linestyles`, `linewidths`, and `markers` is either
    one value, used for every label, a sequence with one entry per label,
    or ``None``, which leaves that attribute at matplotlib's default. One
    value is a string, a number, or a dash tuple such as ``(0, (5, 2))``;
    a list, or a tuple of other items, is one entry per label. Time
    and space are :math:`O(n)` in the number of labels.

    Parameters
    ----------
    labels : sequence of str
        The legend entries.
    linestyles : str or sequence, optional
        The line style of each handle, such as ``'-'`` or ``'--'``.
    linewidths : float or sequence of float, optional
        The line width of each handle, in points.
    markers : str or sequence, optional
        The marker of each handle, such as ``'o'`` or ``'^'``.
    color : matplotlib color, default 'black'
        The color of every handle.
    **line_kwargs
        Passed to every :class:`~matplotlib.lines.Line2D`, such as
        ``markersize``.

    Returns
    -------
    list of matplotlib.lines.Line2D
        The handles, one per label.

    Raises
    ------
    ValueError
        If a sequence among `linestyles`, `linewidths`, and `markers` does
        not have one entry per label.

    Examples
    --------
    >>> lines = line_handles(["ring", "linear"], linestyles=["-", "--"],
    ...                      markers="o", linewidths=2)
    >>> [(line.get_label(), line.get_linestyle(), line.get_marker())
    ...  for line in lines]
    [('ring', '-', 'o'), ('linear', '--', 'o')]
    """
    n_labels = len(labels)
    styles = _per_label("linestyles", linestyles, n_labels)
    widths = _per_label("linewidths", linewidths, n_labels)
    marks = _per_label("markers", markers, n_labels)
    handles = []
    for label, style, width, mark in zip(labels, styles, widths, marks):
        kwargs = {"linestyle": style, "linewidth": width, "marker": mark}
        kwargs = {k: v for k, v in kwargs.items() if v is not None}
        handles.append(
            Line2D([], [], label=label, color=color, **kwargs,
                   **line_kwargs))
    return handles


def move_legend_to_axes(
    legend: Legend,
    axes: Sequence[Axes],
    locs: Sequence[Any],
    *,
    title: Optional[str] = None,
    **legend_kwargs: Any,
) -> list[Legend]:
    """
    Remove `legend` and re-create it inside each of `axes`.

    The new legends carry the handles, labels, and title of `legend`; each
    replaces the legend its axes had. This moves, for instance, the legend
    a seaborn facet grid draws beside its panels into one or more panels.
    Time and space are :math:`O(n_\\mathrm{axes} \\, n_\\mathrm{entries})`.

    Parameters
    ----------
    legend : matplotlib.legend.Legend
        The legend to move; it is removed from its figure or axes.
    axes : sequence of matplotlib.axes.Axes
        The axes that receive a copy each.
    locs : sequence
        One legend location per axes, as the ``loc`` argument of
        :meth:`matplotlib.axes.Axes.legend` takes it.
    title : str, optional
        A new title; by default the title of `legend` is kept.
    **legend_kwargs
        Passed to every :meth:`matplotlib.axes.Axes.legend` call, such as
        ``frameon``.

    Returns
    -------
    list of matplotlib.legend.Legend
        The new legends, in the order of `axes`.

    Raises
    ------
    ValueError
        If `axes` and `locs` differ in length.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt
    >>> fig, (left, right) = plt.subplots(ncols=2)
    >>> _ = left.plot([0, 1], [0, 1], label="rise")
    >>> outer = fig.legend(title="curves")
    >>> inner = move_legend_to_axes(outer, [right], ["upper left"])
    >>> fig.legends
    []
    >>> [text.get_text() for text in right.get_legend().get_texts()]
    ['rise']
    >>> right.get_legend().get_title().get_text()
    'curves'
    >>> plt.close(fig)
    """
    if len(axes) != len(locs):
        raise ValueError(
            f"locs has {len(locs)} entries for {len(axes)} axes.")
    entries = [(handle, text.get_text())
               for handle, text in zip(legend.legend_handles,
                                       legend.get_texts())
               if handle is not None]
    handles = [handle for handle, _ in entries]
    labels = [label for _, label in entries]
    new_title = legend.get_title().get_text() if title is None else title
    legend.remove()
    return [ax.legend(handles, labels, loc=loc, title=new_title,
                      **legend_kwargs)
            for ax, loc in zip(axes, locs)]


def band_under_lines(
    ax: Axes,
    bands: Sequence[tuple[Any, Any]],
    *,
    alpha: float = 0.25,
    **fill_kwargs: Any,
) -> list[PolyCollection]:
    """
    Shade one band under each line of `ax`, in that line's color.

    ``bands[i]`` is the pair ``(lower, upper)`` of the band under
    ``ax.lines[i]``, each an array over that line's x data; the band is
    drawn with :meth:`matplotlib.axes.Axes.fill_between`. Lines beyond
    ``len(bands)`` get no band. Time and space are :math:`O(n_\\mathrm{lines}
    \\, n_\\mathrm{points})`.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes whose lines are shaded, in drawing order.
    bands : sequence of pairs of array_like
        The lower and upper edge of each band, shape ``(n_points,)`` each,
        matching the x data of the line.
    alpha : float, default 0.25
        The opacity of every band.
    **fill_kwargs
        Passed to every :meth:`matplotlib.axes.Axes.fill_between` call,
        such as ``linewidth``.

    Returns
    -------
    list of matplotlib.collections.PolyCollection
        The bands, in the order of the lines.

    Raises
    ------
    ValueError
        If there are more bands than lines, or a band edge's length differs
        from its line's number of points.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> _ = ax.plot([0, 1, 2], [1, 1, 1], color="#4477AA")
    >>> bands = band_under_lines(ax, [([0.5, 0.5, 0.5], [1.5, 1.5, 1.5])])
    >>> len(bands)
    1
    >>> bands[0].get_facecolor()[0].round(2).tolist()
    [0.27, 0.47, 0.67, 0.25]
    >>> plt.close(fig)
    """
    lines = list(ax.lines)
    if len(bands) > len(lines):
        raise ValueError(
            f"{len(bands)} bands were given for {len(lines)} lines.")
    collections: list[PolyCollection] = []
    for line, (lower, upper) in zip(lines, bands):
        x = np.asarray(line.get_xdata(), dtype=float)
        lower = np.asarray(lower, dtype=float)
        upper = np.asarray(upper, dtype=float)
        if lower.shape != x.shape or upper.shape != x.shape:
            raise ValueError(
                f"a band has edges of shape {lower.shape} and "
                f"{upper.shape} for a line of {x.shape[0]} points.")
        collections.append(
            ax.fill_between(x, lower, upper, color=line.get_color(),
                            alpha=alpha, **fill_kwargs))
    return collections
