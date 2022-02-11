"""
a collection of classes and functions used by other modules and submodules.
"""

import numpy as np
import matplotlib as mpl


def round_down_first_non_zero(x: float) -> float:
    """rounds down a number to its first non-zero digit.

    Parameters:
    x (float or int): number which is rounded.

    Returns:
    a float or integer number to which to which x is rounded down to its \
        first non-zero digit.
    """
    if x == 0:
        return x
    else:
        exponent = np.floor(np.log10(abs(x)))
        non_zero = 10 ** exponent
        return round(np.floor(x/non_zero)*non_zero, abs(exponent))


def round_up_nearest(dividend: float, diviser: float) -> float:
    """rounds up the floating-point number dividend by the diviser.

    Parameters:
    dividend (float): The number should be rounded up.
    diviser (float): The number used as the diviser.

    Return
    A floating-point number that is a ceil for dividend and is divisible \
        by diviser.
    """
    return np.ceil(dividend / diviser) * diviser


def color_patcher(colors):
    """return color patches for legend.

    Parameters:
    color: matplotlib color object

    Return:
    color_patches; matplotlib color patch object
    """
    color_patches = []
    for kwargs in mpl.cycler(color=colors):
        color_patches.append(mpl.patches.Patch(**kwargs))
    return color_patches


def yticks(axis, limits, code=False, decimals=3, **kwargs):
    """set yticks of of a given axis based on the lower and upper limits \
    as well as the step size between ticks and a frame-like limit below \
    the lower limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in \
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float, float, float)):
    the first number is lower limit, the second one is the upper limit, \
    the 3rd is the step size  between ticks, and the last one is a tolerance
    code  (bool): whether shows the tick lables or not.
    decimals (int): number of decimal in tick labels.
    **kwargs: the same as the kwargs for set_yticklabels

    Return: --
    """
    lower, upper, step = limits
    tolerance = round_down_first_non_zero(step / 10)
    if (step >= abs(upper - lower)) or (step <= 0):
        raise ValueError(
            f"'{step}' "
            "should be positive and smaller than the absolute difference \
            of lower and upper limits.")
    if (tolerance >= step) or (tolerance < 0):
        raise ValueError(
            f"'{tolerance}' "
            "should be positive and smaller than the absolute values of \
            lower and upper limits.")
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
    """set xticks of of a given axis based on the lower and upper limits \
    as well as the step size between ticks and a frame-like limit below \
    the lower limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in \
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float, float, \
    float)): the first number is lower limit, the second one is the upper \
    limit, the 3rd is the step size  between ticks, and the last one is a \
    tolerance.
    code  (bool): whether shows the tick lables or not.
    decimals (int): number of decimal in tick labels.
    **kwargs: the same as the kwargs for set_xticklabels

    Return: -
    """
    lower, upper, step = limits
    tolerance = round_down_first_non_zero(step / 10)
    if (step >= abs(upper - lower)) or (step <= 0):
        raise ValueError(
            f"'{step}' "
            "should be positive and smaller than the absolute \
            difference of lower and upper limits.")
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
    """sets yticks of of a given axis based on the lower and upper limits \
    as well as the number of ticks and a frame-like limit below the lower \
    limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in \
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float, \
    int, float)): the first number is lower limit, the secone one is the \
    upper limit, the 3rd is the number of ticks, and the last one is \
    a tolerance.
    code  (bool): whether shows the tick lables or not.
    decimals (int): number of decimal in tick labels.
    **kwargs: the same as the kwargs for set_yticklabels

    Return: -

    Requirements:
    matplotlib
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
    """set xticks of of a given axis based on the lower and upper limits \
    as well as the number of ticks and a frame-like limit below the lower \
    limit and above the upper limit.

    Caution: The numbers are not rounded based on standard convention in \
    which ticks are multiples of 10^m where m is an integer number.

    Parameters:
    axis (matplotlib.axes.Axes): The axis for which the ticks are changed.
    limits (a tuple for 4 numbers in this order: (float, float, \
    int, float)): the first number is lower limit, the secone one is the \
    upper limit, the 3rd is the number of ticks, and the last one is a \
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
    The legend names are combined in one list of strings. The order of names \
    in this list is the order by which the hue, style, and size keywords are \
    defined in saeborn.lineplot. For each keyword, the sub-list starts with \
    the title name of that keyword in the legend table and then the names \
    of the various values of that keyword.

    Parameters:
    line (matplotlib.axes.Axes): The line that contains the old legend names.
    legend_names (list of str): A list in which each item is either name of \
    a seaborn lineplot keyword (hue, style, size) or an item for the keyword \
    -- See the caution above.
    **kwargs : The standard kwargs for matplotlib.legend method.

    Return: -

    Requirements:
    Matplotlib, Seaborn.
    """
    legned_old_names = line.legend(
        fontsize=16, bbox_to_anchor=(1.005, 1), loc=2, edgecolor='black',
        title_fontsize=16, markerscale=1.5, **kwargs).texts
    for new_name, line in zip(legned_old_names, legend_names):
        new_name.set_text(line)


def ls_handler(attributes, linstyles, color='black', lw=2, **kwargs):
    """creates handles for matplotlib legend functions based on the \
    line styles. It recieves a list  physical attributes and a list of \
    linestyles (the two list have the same size) and returns a list of \
    handles.
    """
    handles = []
    for attr, linestyle in zip(attributes, linstyles):
        mline = mpl.lines.Line2D(
            [], [], ls=linestyle, label=attr, color=color, lw=lw, **kwargs)
        handles.append(mline)
    return handles


def lw_handler(attributes, linewidths, color='block', **kwargs):
    """creates handles for matplotlib legend functions based on the \
    line widths. It recieves a list  physical attributes and a list of \
    linewidths (the two list have the same size) and returns a list of \
    handles.
    """
    handles = []
    for attr, linewidth in zip(attributes, linewidths):
        mline = mpl.lines.Line2D(
            [], [], lw=linewidth, label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles


def marker_handler(attributes, markers, color='black', **kwargs):
    """creates handles for matplotlib legend functions based on the \
    line widths. It recieves a list  physical attributes and a list of \
    linewidths (the two list have the same size) and returns a list of \
    handles.
    """
    handles = []
    for attr, marker in zip(attributes, markers):
        mline = mpl.lines.Line2D(
            [], [],  marker=marker, label=attr, color=color, **kwargs)
        handles.append(mline)
    return handles
