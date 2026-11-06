"""
==========================================================
:mod:`polyphys.viz.colors`
==========================================================

The :mod:`polyphys.viz.colors` module holds the qualitative color lists the
studies' figures use, as tuples of ``#RRGGBB`` strings, and
:func:`truncate_colormap`, which cuts a matplotlib color map down to a
sub-range of its domain so that, for instance, the pale end of a sequential
map does not vanish against a white background.

The module imports no plotting library at import time, so the color lists
can be read without the ``viz`` extra and passed to any plotting call that
accepts matplotlib color strings; :func:`truncate_colormap` imports
matplotlib when called and raises :exc:`ImportError` with the install hint
when it is missing.

Sources
=======

``TOL_BRIGHT``, ``TOL_VIBRANT``, ``TOL_MUTED``
    Paul Tol, *Colour Schemes*, SRON Technical Note SRON/EPS/TN/09-002,
    issue 3.2, 18 August 2021, Figures 1, 3, and 4
    (https://sronpersonalpages.nl/~pault/data/colourschemes.pdf). The
    schemes are designed to be distinct for color-blind readers, distinct
    from black and white, and distinct on screen and paper. Each tuple is
    in the order the note gives for picking colors in a fixed sequence
    (page 6), and ends with the scheme's gray, which the note reserves for
    bad or missing data.
``DUTCH_FIELD``
    A nine-color palette named "Dutch Field" in a HEAVY.AI (formerly
    OmniSci) blog post, "12 color palettes for telling better stories with
    your data", which could not be opened when this module was written. The
    values here are the ones listed, under the same name, by
    M. Jovanovic, "12 Great Data Visualization Color Palettes to Use",
    wpDataTables, 17 December 2023
    (https://wpdatatables.com/data-visualization-color-palette/), read on
    2026-11-06. Neither page states terms of reuse.
``STUDY_COLORS``
    The maintainer's own twelve-color list, assembled from the lists above
    and CSS named colors for the published figures; its fifth entry is
    matplotlib's ``gray``, ``#808080``.

Examples
========

>>> from polyphys.viz.colors import TOL_BRIGHT, truncate_colormap
>>> TOL_BRIGHT[0]
'#4477AA'
>>> len(TOL_BRIGHT)
7
>>> truncate_colormap("viridis", 0.2, 0.8, n=5).name
'trunc(viridis,0.2,0.8)'
"""
from typing import TYPE_CHECKING, Union

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - imported for the annotations only
    from matplotlib.colors import Colormap, LinearSegmentedColormap

__all__ = [
    "DUTCH_FIELD",
    "STUDY_COLORS",
    "TOL_BRIGHT",
    "TOL_MUTED",
    "TOL_VIBRANT",
    "truncate_colormap",
]

#: Tol's *bright* scheme, his default for lines and their labels: blue, red,
#: green, yellow, cyan, purple, grey.
TOL_BRIGHT: tuple[str, ...] = (
    "#4477AA",
    "#EE6677",
    "#228833",
    "#CCBB44",
    "#66CCEE",
    "#AA3377",
    "#BBBBBB",
)

#: Tol's *vibrant* scheme, designed for TensorBoard: orange, blue, cyan,
#: magenta, red, teal, grey.
TOL_VIBRANT: tuple[str, ...] = (
    "#EE7733",
    "#0077BB",
    "#33BBEE",
    "#EE3377",
    "#CC3311",
    "#009988",
    "#BBBBBB",
)

#: Tol's *muted* scheme, with more colors but no clear red or medium blue:
#: rose, indigo, sand, green, cyan, wine, teal, olive, purple, pale grey.
TOL_MUTED: tuple[str, ...] = (
    "#CC6677",
    "#332288",
    "#DDCC77",
    "#117733",
    "#88CCEE",
    "#882255",
    "#44AA99",
    "#999933",
    "#AA4499",
    "#DDDDDD",
)

#: The nine-color "Dutch Field" palette; see the module docstring for its
#: provenance.
DUTCH_FIELD: tuple[str, ...] = (
    "#e60049",
    "#0bb4ff",
    "#50e991",
    "#e6d800",
    "#9b19f5",
    "#ffa300",
    "#dc0ab4",
    "#b3d4ff",
    "#00bfa0",
)

#: The maintainer's twelve-color list for the published figures.
STUDY_COLORS: tuple[str, ...] = (
    "#00008b",
    "#e60049",
    "#9b19f5",
    "#228b22",
    "#808080",
    "#ffe119",
    "#ff00ff",
    "#00ffff",
    "#ffa300",
    "#00ff00",
    "#8b4513",
    "#000000",
)


def truncate_colormap(
    cmap: "Union[str, Colormap]",
    lower: float = 0.0,
    upper: float = 1.0,
    n: int = 256,
) -> "LinearSegmentedColormap":
    """
    Cut `cmap` down to the sub-range ``[lower, upper]`` of its domain.

    The result is a linear segmented color map through the `n` colors
    that `cmap` takes at equally spaced points of ``[lower, upper]``, so its
    value at 0 is ``cmap(lower)`` and its value at 1 is ``cmap(upper)``.
    Time and space are :math:`O(n)`.

    Parameters
    ----------
    cmap : str or matplotlib.colors.Colormap
        The color map to cut, or the name of a registered one, such as
        ``'viridis'``.
    lower, upper : float, default 0.0 and 1.0
        The sub-range of ``[0, 1]`` to keep; ``0 <= lower < upper <= 1``.
    n : int, default 256
        The number of colors sampled from `cmap`, at least 2.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        The truncated map, named ``'trunc(<cmap.name>,<lower>,<upper>)'``.

    Raises
    ------
    ImportError
        If matplotlib is not installed; the ``viz`` extra installs it.
    ValueError
        If the range is not ``0 <= lower < upper <= 1`` or `n` is below 2.

    Examples
    --------
    >>> import matplotlib
    >>> dark = truncate_colormap("viridis", 0.5, 1.0, n=4)
    >>> dark.name
    'trunc(viridis,0.5,1.0)'
    >>> dark(0.0) == matplotlib.colormaps["viridis"](0.5)
    True
    """
    try:
        import matplotlib
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError as error:  # pragma: no cover - depends on the install
        raise ImportError(
            "polyphys.viz.colors.truncate_colormap needs matplotlib; "
            "install it with python -m pip install 'polyphys[viz]'."
        ) from error
    if isinstance(cmap, str):
        cmap = matplotlib.colormaps[cmap]
    if not 0.0 <= lower < upper <= 1.0:
        raise ValueError(
            f"the range must satisfy 0 <= lower < upper <= 1, got "
            f"lower={lower} and upper={upper}.")
    if n < 2:
        raise ValueError(f"n must be at least 2, got {n}.")
    colors = cmap(np.linspace(lower, upper, n))
    return LinearSegmentedColormap.from_list(
        f"trunc({cmap.name},{lower},{upper})", colors, N=n)
