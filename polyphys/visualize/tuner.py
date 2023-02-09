"""
a collection of classes and functions used by other modules and submodules.
"""
from typing import List
import matplotlib.colors as mplc
from ..manage.utilizer import round_down_first_non_zero
import numpy as np
import matplotlib as mpl
import seaborn as sns
import inspect


def add_legend_sub_axis(
    grid: sns.FacetGrid,
    axes_idx: List[int],
    locs: List[str],
    **kwargs
) -> sns.FacetGrid:
    """
    Remove the current legend in a seaborn facte grid and add it to that
    grid's `axes` at given locations `loc`.
    """
    old_legend = grid.legend
    # Extract the components of the legend we need to reuse
    handles = old_legend.legendHandles
    labels = [t.get_text() for t in old_legend.get_texts()]
    # Extract legend properties that can be passed to the recreation method
    # (Vexingly, these don't all round-trip)
    legend_kws = inspect.signature(mpl.legend.Legend).parameters
    props = {
        k: v for k, v in old_legend.properties().items() if k in legend_kws
    }
    # Delegate default bbox_to_anchor rules to matplotlib
    props.pop("bbox_to_anchor")
    # Try to propagate the existing title and font properties;
    # respect new ones too
    title = props.pop("title")
    if "title" in kwargs:
        title.set_text(kwargs.pop("title"))
    title_kwargs = {k: v for k, v in kwargs.items() if k.startswith("title_")}
    for key, val in title_kwargs.items():
        title.set(**{key[6:]: val})
        kwargs.pop(key)
    # Try to respect the frame visibility
    kwargs.setdefault("frameon", old_legend.legendPatch.get_visible())
    # Remove the old legend and create the new one
    props.update(kwargs)
    old_legend.remove()
    for (row, col), loc in zip(axes_idx, locs):
        ax_legend = grid.axes[row, col].legend(handles, labels, loc=loc, **props)
        grid.axes[row, col].add_artist(ax_legend)
    # Let the Grid object continue to track the correct legend object
    return grid


def set_facet_grid_legend(
    facet_grid: sns.FacetGrid,
    new_title: str,
    new_labels: List[str]
) -> None:
    """Sets new labels in a facet gird plot.

    Parameters
    ----------
    facet_grid: Seaborn FacetGrid
        A FacetGrid object
    new_title: str
        A new title for the legend.
    new_labels: list of str
        A list of new labels.
    """
    # check axes and find which is have legend
    for ax in facet_grid.axes.flat:
        legend = facet_grid.axes.flat[0].get_legend()
        if legend is not None:
            break
    # or legend may be on a figure
    if legend is None:
        legend = facet_grid._legend

    # change legend texts
    legend.set_title(new_title)
    for text, label in zip(legend.texts, new_labels):
        text.set_text(label)


def color_patcher(colors, alpha=1.0):
    """return color patches for legend.

    Parameters:
    color: matplotlib color object

    Return:
    color_patches; matplotlib color patch object
    """
    color_patches = []
    for kwargs in mpl.cycler(color=colors):
        color_patches.append(mpl.patches.Patch(**kwargs, alpha=alpha))
    return color_patches


def yticks(axis, limits, code=False, decimals=3, **kwargs):
    """set yticks of a given axis based on the lower and upper limits
    as well as the step size between ticks and a frame-like limit below
    the lower limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in
    which ticks are multiples of 10^m where m is an integer number.

    Parameters
    ----------
    axis: matplotlib.axes.Axes
        The axis for which the ticks are changed.
    limits: tuple of 3 floats
        The first number is lower limit, the second one is the upper limit,
        the 3rd is the step size  between ticks.
    code: bool
        Whether shows the tick lables or not.
    decimals: int
        Number of decimal in tick labels.
    **kwargs:
        The same as the kwargs for matplotlib.pyplot.Axis.set_yticklabels
    """
    lower, upper, step = limits
    tolerance = round_down_first_non_zero(step / 10)
    if (step >= abs(upper - lower)) or (step <= 0):
        raise ValueError(
            f"'{step}'"
            " should be positive and smaller than the absolute difference of"
            " lower and upper limits."
        )
    if (tolerance >= step) or (tolerance < 0):
        raise ValueError(
            f"'{tolerance}'"
            " should be positive and smaller than the absolute values of"
            " lower and upper limits."
        )
    axis.set_ylim(lower - tolerance, upper + tolerance)
    axis.set_yticks([
        0 if (i < tolerance) and (i > -1 * tolerance) else i for i in
        np.arange(lower, upper+step, step)])
    if code is False:
        axis.set_yticklabels(labels=[])
    else:
        axis.set_yticklabels(
            labels=np.around(axis.get_yticks(), decimals=decimals), **kwargs)


def xticks(axis, limits, code=False, decimals=3, **kwargs):
    """set xticks of a given axis based on the lower and upper limits
    as well as the step size between ticks and a frame-like limit below
    the lower limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in
    which ticks are multiples of 10^m where m is an integer number.

    Parameters
    ----------
    axis: matplotlib.axes.Axes
        The axis for which the ticks are changed.
    limits: tuple of 3 floats
        The first number is lower limit, the second one is the upper limit,
        the 3rd is the step size  between ticks.
    code: bool
        Whether shows the tick lables or not.
    decimals: int
        Number of decimal in tick labels.
    **kwargs:
        The same as the kwargs for matplotlib.pyplot.Axis.set_xticklabels
    """
    lower, upper, step = limits
    tolerance = round_down_first_non_zero(step / 10)
    if (step >= abs(upper - lower)) or (step <= 0):
        raise ValueError(
            f"'{step}'"
            ",should be positive and smaller than the absolute difference of"
            "lower and upper limits."
        )
    axis.set_xlim(lower - tolerance, upper + tolerance)
    axis.set_xticks([
        0 if (i < tolerance) and (i > -1 * tolerance) else i for i in
        np.arange(lower, upper+step, step)])
    if code is False:
        axis.set_xticklabels(labels=[])
    else:
        axis.set_xticklabels(
            labels=np.around(axis.get_xticks(), decimals=decimals), **kwargs)


def yticks_nstep(axis, limits, code=False, decimals=3, **kwargs):
    """sets yticks of of a given axis based on the lower and upper limits
    as well as the number of ticks and a frame-like limit below the lower
    limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float,
    int, float)): the first number is lower limit, the secone one is the
    upper limit, the 3rd is the number of ticks, and the last one is
    a tolerance.
    code  (bool): whether shows the tick lables or not.
    decimals (int): number of decimal in tick labels.
    **kwargs: the same as the kwargs for set_yticklabels
    """
    lower, upper, nstep = limits
    if not isinstance(nstep, int):
        raise TypeError(
            f"'{nstep}' "
            "is not an integer.")
    step = round_down_first_non_zero(abs(upper - lower) / nstep)
    tolerance = round_down_first_non_zero(step / 10)
    axis.set_ylim(lower - tolerance, upper + tolerance)
    axis.set_yticks([
        0 if (i < tolerance) and (i > -1 * tolerance) else i for i in
        np.linspace(lower, upper, nstep)])
    if code is False:
        axis.set_yticklabels(labels=[])
    else:
        axis.set_yticklabels(
            labels=np.around(axis.get_yticks(), decimals=decimals), **kwargs)


def xticks_nstep(axis, limits, code=False, decimals=3, **kwargs):
    """set xticks of of a given axis based on the lower and upper limits
    as well as the number of ticks and a frame-like limit below the lower
    limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float,
    int, float)): the first number is lower limit, the secone one is the
    upper limit, the 3rd is the number of ticks, and the last one is a
    tolerance.
    code  (bool): whether shows the tick lables or not.
    decimals (int): number of decimal in tick labels.
    **kwargs: the same as the kwargs for set_xticklabels

    Return: -
    """
    lower, upper, nstep = limits
    if not isinstance(nstep, int):
        raise TypeError(f"'{nstep}' is not an integer.")
    step = round_down_first_non_zero(abs(upper - lower) / nstep)
    tolerance = round_down_first_non_zero(step / 10)
    axis.set_xlim(lower - tolerance, upper + tolerance)
    axis.set_xticks([
        0 if (i < tolerance) and (i > -1*tolerance) else i for i in
        np.linspace(lower, upper, nstep)])
    if code is False:
        axis.set_xticklabels(labels=[])
    else:
        axis.set_xticklabels(
            labels=np.around(axis.get_xticks(), decimals=decimals), **kwargs)


def change_legend_name(line, legend_names, **kwargs):
    """changes the legend lables in a seaborn born axis.

    Caution:
    The legend names are combined in one list of strings. The order of names
    in this list is the order by which the hue, style, and size keywords are
    defined in saeborn.lineplot. For each keyword, the sub-list starts with
    the title name of that keyword in the legend table and then the names
    of the various values of that keyword.

    Parameters:
    line (matplotlib.axes.Axes): The line that contains the old legend names.
    legend_names (list of str): A list in which each item is either name of
    a seaborn lineplot keyword (hue, style, size) or an item for the keyword
    -- See the caution above.
    **kwargs : The standard kwargs for matplotlib.legend method.
    """
    legned_old_names = line.legend(
        fontsize=16, bbox_to_anchor=(1.005, 1), loc=2, edgecolor='black',
        title_fontsize=16, markerscale=1.5, **kwargs).texts
    for new_name, line in zip(legned_old_names, legend_names):
        new_name.set_text(line)


def color_handler(attributes, colors, **kwargs):
    """creates handles for matplotlib legend functions based on the
    line styles. It recieves a list  physical attributes and a list of
    linestyles (the two list have the same size) and returns a list of
    handles.
    """
    handles = []
    for attr, color in zip(attributes, colors):
        mline = mpl.patches.Patch(label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles


def ls_handler(attributes, linstyles, color='black', lw=2, **kwargs):
    """creates handles for matplotlib legend functions based on the
    line styles. It recieves a list  physical attributes and a list of
    linestyles (the two list have the same size) and returns a list of
    handles.
    """
    handles = []
    for attr, linestyle in zip(attributes, linstyles):
        mline = mpl.lines.Line2D(
            [], [], ls=linestyle, label=attr, color=color, lw=lw, **kwargs)
        handles.append(mline)
    return handles


def lw_handler(attributes, linewidths, color='block', **kwargs):
    """creates handles for matplotlib legend functions based on the
    line widths. It recieves a list  physical attributes and a list of
    linewidths (the two list have the same size) and returns a list of
    handles.
    """
    handles = []
    for attr, linewidth in zip(attributes, linewidths):
        mline = mpl.lines.Line2D(
            [], [], lw=linewidth, label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles


def marker_handler(attributes, markers, color='black', **kwargs):
    """creates handles for matplotlib legend functions based on the
    line widths. It recieves a list  physical attributes and a list of
    linewidths (the two list have the same size) and returns a list of
    handles.
    """
    handles = []
    for attr, marker in zip(attributes, markers):
        mline = mpl.lines.Line2D(
            [], [],  marker=marker, label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles


def marker_ls_handler(
    attributes,
    linestyles,
    markers,
    color='black',
    **kwargs
) -> list:
    """creates handles for matplotlib legend functions based on the markers and
    linestyles. It recieves a list  physical attributes and a list of
    linewidths (the two list have the same size) and returns a list of
    handles.
    """
    handles = []
    for attr, ls, marker in zip(attributes, linestyles, markers):
        mline = mpl.lines.Line2D(
            [], [],  ls=ls, marker=marker, label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles


def truncated_colormap(
    cmap: mplc.Colormap,
    min_value: float = 0.0,
    max_value: float = 1.0,
    ncolors: float = 200
) -> mplc.Colormap:
    """create a linear segmented color map from a given color map between a
    pair min and max values into n points. This function allows the user to
    created discrete color maps from contiunous ones.

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap
        an equidistantly (matplotlib or seaborn) color map in the range [0,1]
        where i.e. 0 maps to colors[0] and 1 maps to colors[-1] of the colors
        in the cmap.
    min_value: float, default 0.0
        The minimum value of the truncated cmap in the range [0,1]
    max_value: float, default 1.0
        The maximum value of the truncated cmap in the range [0,1]
    ncolors: int, default 200
        Number of the colors in the truncated cmap.

    Return
    ------
    matplotlib.colors.ColorMap:
        A truncated color map.
    """
    return mplc.LinearSegmentedColormap.from_list(
        f'trunc({cmap.name},{min_value},{max_value})',
        cmap(
            np.linspace(
                min_value, max_value, ncolors
                )
            )
        )
