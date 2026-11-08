"""
==========================================================
:mod:`polyphys.analyze.density`
==========================================================


The :mod:`polyphys.analyze.density` module turns a histogram of bead centers
along one coordinate into two profiles: the local number density, which
counts centers per bin volume, and the local volume fraction, which counts
the volume the beads occupy. The second is the harder one, because a bead of
finite size reaches into the bins next to the one that holds its center.
Every routine is a numpy or scipy kernel on arrays: reading frames, binning
positions, and averaging over frames are the caller's work, and the only
geometric inputs are the bin edges, the bead diameter, and the box lengths.

Notation
========

:math:`x`
    The binned coordinate: the axial position :math:`z` of a cylinder, or
    the distance :math:`r` from the axis of a cylinder or from a point in a
    cube. Lengths are in the caller's unit, the monomer diameter in this
    package.
:math:`e_0 < e_1 < \\dots < e_B`
    The bin edges, ``edges``; bin :math:`i` is :math:`[e_i, e_{i+1})`, so
    there are :math:`B` bins, ``n_bins``.
:math:`V_i`
    The volume of bin :math:`i`, ``volumes``: a slab
    :math:`\\pi D^2 (e_{i+1} - e_i) / 4`, a cylindrical shell
    :math:`\\pi L (e_{i+1}^2 - e_i^2)`, or a spherical shell
    :math:`4 \\pi (e_{i+1}^3 - e_i^3) / 3`.
:math:`H_i`, :math:`F`
    The number of bead centers counted in bin :math:`i` over :math:`F`
    frames, ``counts`` and ``n_frames``.
:math:`\\rho_i = H_i / (F V_i)`
    The local number density of bin :math:`i`, ``rho``, in beads per cubic
    length unit.
:math:`a`, :math:`R_b = a / 2`, :math:`v = \\pi a^3 / 6`
    The bead diameter ``bead_diameter``, its radius, and its volume.
:math:`s(x \\to j)`
    The *share kernel*: the volume of one bead whose center is at coordinate
    :math:`x` that lies inside bin :math:`j`. It sums to :math:`v` over the
    bins when the bead is inside the domain.
:math:`S_{ij}`
    The *share matrix*, ``share``: the volume of a test bead centered in bin
    :math:`i` that lies in bin :math:`j`, with the bead at the bin center
    (``shares='center'``) or averaged over the bin (``shares='average'``).
:math:`\\phi_i = \\sum_j S_{ij} \\rho_j`
    The local volume fraction of bin :math:`i`, ``phi``, dimensionless.

Contents
========

- :func:`bin_volumes`: the closed-form bin volumes.
- :func:`share_matrix`: the share matrix, built once per ``(edges,
  bead_diameter)``.
- :func:`sphere_cylinder_intersection`: the volume of a ball inside a
  coaxial cylinder, the kernel of the radial direction in a cylinder.
- :func:`number_density` and :func:`volume_fraction`: the two profiles from
  counts, volumes, and the share matrix.
- :func:`normalize_profile`: an explicit normalization.
- :func:`density_profiles` and :class:`Profiles`: the one-call form.

The volume fraction
===================

The volume fraction field is the indicator of "inside some bead", and in the
continuum its average over a region is the density convolved with a ball of
radius :math:`R_b`. With the density piecewise constant on the bins, the
value at a point :math:`x` is :math:`\\sum_j \\rho_j\\, s(x \\to j)`. The
research code evaluated this at the bin center, :math:`\\phi_i = \\sum_j
\\rho_j\\, s(c_i \\to j)` (thesis, appendix "Discrete local distributions
in the cylindrical coordinate system", eq. 7); ``shares='center'`` keeps
that form. The default, ``shares='average'``, averages the kernel over the
bead's position in its bin with the bin's own measure,

.. math::

    \\bar{s}(i \\to j) = \\frac{1}{V_i} \\int_{e_i}^{e_{i+1}}
    s(x \\to j)\\, w(x)\\, dx, \\qquad
    w(x)\\,dx = A\\,dx,\\; 2 \\pi L x\\,dx,\\; 4 \\pi x^2\\,dx,

which is the exact bin average under the one assumption already made for
:math:`\\rho`. Because the indicator is symmetric in the two points,
:math:`V_i\\, \\bar{s}(i \\to j) = V_j\\, \\bar{s}(j \\to i)`, and from
that one identity the averaged form conserves the total bead volume,
:math:`\\sum_i V_i \\phi_i = \\sum_j H_j v / F`, and returns
:math:`\\phi = \\rho v` for a uniform fluid, both to round-off. The
centered form gets the second right and misses the first by
:math:`O(h^2)` along :math:`r`: up to 2.65 % of a bead volume for
cylindrical shells and 7.27 % for spherical shells at a bin width
:math:`h = a / 4`, and 0.71 % and 1.88 % at :math:`a / 8`, largest in the
innermost bins. Both forms return :math:`\\rho v` only in the bins whose
test bead stays inside the last edge; in the outermost radial bins the
bead reaches past it and the row sums of the share matrix fall short of
:math:`v`.

The kernels
===========

A slab's kernel is a spherical segment
(:func:`~polyphys.analyze.measurer.spherical_segment`), and the averaged slab
kernel is its second antiderivative in closed form. A spherical shell's kernel
is a difference of sphere-sphere intersection volumes
(:func:`~polyphys.analyze.measurer.sphere_sphere_intersection`). A cylindrical
shell's kernel is a difference of ball-in-cylinder volumes,
:func:`sphere_cylinder_intersection`, which has no closed form and is
integrated numerically (``intersection='exact'``); the research code used the
sphere-sphere volume in its place (``intersection='sphere'``), whose error is
about :math:`0.047\\, a / R` of a bead volume at an edge of radius :math:`R`,
with :math:`a = 2 R_b` the bead diameter, first order and independent of the
bin width. The shells' averaged kernels are integrated by Gauss-Legendre rules
on the pieces between the kernels' breakpoints, exactly where the integrand is
polynomial.

The axial direction of a cylinder is periodic: a bead that reaches past one
end of the domain reaches into the bins at the other end. The radial
directions are not.

Conventions the kernels do not fix
==================================

:func:`bin_volumes` takes the slab diameter the caller passes.
:func:`density_profiles` passes :math:`D_\\mathrm{eff} = D - a`, the
diameter of the region a bead center can reach, so that a uniform crowder
fluid returns the bulk volume fraction as the parsers define it,
:math:`\\phi_c = N_c v_c / (\\pi D_\\mathrm{eff}^2 L / 4)`. Along
:math:`r` no such choice arises. A spherical shell beyond half the box
length of a periodic cube overstates the volume inside the box, so
:func:`bin_volumes` warns there when it is told the box length.

Dependencies
============
- `numpy` and `scipy.sparse`.
- :mod:`polyphys.analyze.measurer`:
  :func:`~polyphys.analyze.measurer.spherical_segment` and
  :func:`~polyphys.analyze.measurer.sphere_sphere_intersection`.
- :mod:`polyphys.manage.types`: :data:`~polyphys.manage.types.GeometryT`
  and :data:`~polyphys.manage.types.DirectionT`.

References
==========
- Weisstein, Eric W. "Spherical Segment." From MathWorld--A Wolfram Web
  Resource, eq. 14. https://mathworld.wolfram.com/SphericalSegment.html
- Weisstein, Eric W. "Sphere-Sphere Intersection." From MathWorld--A
  Wolfram Web Resource, eq. 16.
  https://mathworld.wolfram.com/Sphere-SphereIntersection.html
- Weisstein, Eric W. "Circle-Circle Intersection." From MathWorld--A
  Wolfram Web Resource, eq. 15.
  https://mathworld.wolfram.com/Circle-CircleIntersection.html
- A. H. Sadeghi Isfahani, PhD thesis, University of Waterloo (2024),
  appendix "Discrete local distributions in the cylindrical coordinate
  system", for the centered form.
- C. Jeon, C. Hyeon, Y. Jung, and B.-Y. Ha, Soft Matter 12, 9786-9796
  (2016), https://doi.org/10.1039/c6sm01924b, for the use of these profiles
  in the density sum rule (:mod:`polyphys.analyze.sumrule`).

.. versionadded:: 0.6.0
"""

import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from scipy.sparse import coo_array, csr_array

from ..manage.types import DirectionT, GeometryT
from .measurer import sphere_sphere_intersection, spherical_segment

__all__ = [
    'IntersectionMode',
    'NormalizeMode',
    'Profiles',
    'ShareMode',
    'bin_volumes',
    'density_profiles',
    'normalize_profile',
    'number_density',
    'share_matrix',
    'sphere_cylinder_intersection',
    'volume_fraction',
]

FloatArray = npt.NDArray[np.float64]

ShareMode = Literal['average', 'center']
"""How the test bead sits in its bin, averaged over the bin or at its
center as in the research code.

.. versionadded:: 0.6.0
"""

IntersectionMode = Literal['exact', 'sphere']
"""The kernel of a cylindrical shell, the ball-in-cylinder volume or the
sphere-sphere volume the research code used in its place.

.. versionadded:: 0.6.0
"""

NormalizeMode = Literal['none', 'sum', 'integral']
"""No normalization, division by the sum over bins, or division by the
integral against the bin volumes.

.. versionadded:: 0.6.0
"""

# The bin scheme behind each supported (geometry, direction): the power p of
# the coordinate in the bin measure w(x) dx = x**p dx (up to a constant).
_SchemeT = Literal['slab', 'cylinder', 'sphere']
_SCHEMES: dict[tuple[str, str], _SchemeT] = {
    ('cylindrical', 'z'): 'slab',
    ('cylindrical', 'r'): 'cylinder',
    ('cubic', 'r'): 'sphere',
}
_MEASURE_POWER: dict[_SchemeT, int] = {'slab': 0, 'cylinder': 1, 'sphere': 2}

# Gauss-Legendre nodes per polynomial piece. The position average of a
# shell's kernel is a polynomial of degree at most 7 on each piece, so 16
# nodes are exact there; the ball-in-cylinder integrand is smooth on each
# piece with a half-power end point, where 32 nodes give about 1e-8.
_N_GAUSS_POSITION = 16
_N_GAUSS_AXIAL = 32


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _check_edges(edges: npt.ArrayLike, scheme: _SchemeT) -> FloatArray:
    """Return `edges` as float64 after checking the bin-edge contract."""
    out = np.asarray(edges, dtype=np.float64)
    if out.ndim != 1 or out.size < 2:
        raise ValueError(
            'edges must be a one-dimensional array of at least two values;'
            f' got shape {out.shape}.'
        )
    if not np.all(np.isfinite(out)):
        raise ValueError('edges must be finite.')
    if not np.all(np.diff(out) > 0):
        raise ValueError('edges must be strictly increasing.')
    if scheme != 'slab' and out[0] < 0:
        raise ValueError(
            f'a radial edge cannot be negative; got edges[0] = {out[0]}.'
        )
    return out


def _check_scheme(geometry: str, direction: str) -> _SchemeT:
    """Return the bin scheme of a supported (geometry, direction) pair."""
    try:
        return _SCHEMES[(geometry, direction)]
    except KeyError:
        supported = ', '.join(f'{g!r}, {d!r}' for g, d in _SCHEMES)
        raise ValueError(
            f"({geometry!r}, {direction!r}) is not supported; the supported"
            f' (geometry, direction) pairs are ({supported}).'
        ) from None


def _check_positive(name: str, value: float | None) -> float:
    """Return `value` as a float after checking that it is given and > 0."""
    if value is None:
        raise ValueError(f'{name} is needed here.')
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be positive; got {value}.')
    return float(value)


def _check_choice(name: str, value: str, choices: Sequence[str]) -> None:
    """Raise if `value` is not one of `choices`."""
    if value not in choices:
        raise ValueError(
            f'{name} must be one of {", ".join(map(repr, choices))};'
            f' got {value!r}.'
        )


# ---------------------------------------------------------------------------
# Bin volumes
# ---------------------------------------------------------------------------


def bin_volumes(
    edges: npt.ArrayLike,
    geometry: GeometryT,
    direction: DirectionT,
    *,
    lcyl: float | None = None,
    dcyl: float | None = None,
    lcube: float | None = None,
) -> FloatArray:
    """
    Compute the volume of each bin along one coordinate, in closed form.

    Parameters
    ----------
    edges : array_like
        The bin edges, shape ``(n_bins + 1,)``, finite and strictly
        increasing; non-negative for a radial direction.
    geometry : :py:obj:`~.GeometryT`
        ``'cylindrical'`` or ``'cubic'``.
    direction : :py:obj:`~.DirectionT`
        ``'r'`` for shells, or ``'z'`` for the slabs of a cylinder.
    lcyl : float, optional
        The cylinder length, needed for ``('cylindrical', 'r')``.
    dcyl : float, optional
        The diameter of the slab, needed for ``('cylindrical', 'z')``. The
        kernel takes it as given; :func:`density_profiles` passes the
        effective diameter :math:`D - a`.
    lcube : float, optional
        The side of the periodic cube, for ``('cubic', 'r')``. When given,
        an outermost edge beyond half the side by more than one bin width
        raises a warning, since such a shell overstates the volume inside
        the box; the volumes are still returned.

    Returns
    -------
    numpy.ndarray
        The bin volumes, shape ``(n_bins,)``, in the cube of the length
        unit: :math:`4\\pi (e_{i+1}^3 - e_i^3) / 3` for a spherical shell,
        :math:`\\pi L (e_{i+1}^2 - e_i^2)` for a cylindrical shell, and
        :math:`\\pi D^2 (e_{i+1} - e_i) / 4` for a slab.

    Raises
    ------
    ValueError
        If `edges` breaks its contract, if the pair is not supported, or if
        the box length the pair needs is missing or not positive.

    Warns
    -----
    UserWarning
        If `lcube` is given and ``edges[-1] > lcube / 2 + max bin width``.

    Notes
    -----
    Time and space are :math:`O(n_\\mathrm{bins})`. The half-side warning
    follows MAICoS's ``DensitySphere``, which warns when its ``rmax``
    exceeds half the smallest box vector and computes the shell volumes
    regardless.

    Examples
    --------
    >>> import numpy as np
    >>> from polyphys.analyze.density import bin_volumes
    >>> np.round(bin_volumes([0, 1, 2], 'cubic', 'r'), 5)
    array([ 4.18879, 29.32153])
    >>> np.round(bin_volumes([0, 1, 2], 'cylindrical', 'r', lcyl=3), 5)
    array([ 9.42478, 28.27433])
    >>> np.round(bin_volumes([-1, 0, 1], 'cylindrical', 'z', dcyl=2), 5)
    array([3.14159, 3.14159])

    .. versionadded:: 0.6.0
    """
    scheme = _check_scheme(geometry, direction)
    e = _check_edges(edges, scheme)
    if scheme == 'sphere':
        volumes = 4.0 * np.pi * np.diff(e**3) / 3.0
        if lcube is not None:
            side = _check_positive('lcube', lcube)
            if e[-1] > side / 2.0 + np.max(np.diff(e)):
                warnings.warn(
                    f'the outermost edge {e[-1]} lies beyond half the box'
                    f' side {side / 2.0}; the spherical shells there'
                    ' overstate the volume inside the box.',
                    UserWarning,
                    stacklevel=2,
                )
    elif scheme == 'cylinder':
        length = _check_positive('lcyl', lcyl)
        volumes = np.pi * length * np.diff(e**2)
    else:
        diameter = _check_positive('dcyl', dcyl)
        volumes = np.pi * diameter**2 * np.diff(e) / 4.0
    return np.asarray(volumes, dtype=np.float64)


# ---------------------------------------------------------------------------
# Kernels
# ---------------------------------------------------------------------------


def _gauss_pieces(
    lo: float, hi: float, breakpoints: npt.ArrayLike, n_points: int
) -> tuple[FloatArray, FloatArray]:
    """
    Return Gauss-Legendre nodes and weights for ``[lo, hi]`` split at the
    breakpoints that fall strictly inside it.
    """
    cuts = np.asarray(breakpoints, dtype=np.float64)
    cuts = np.unique(cuts[(cuts > lo) & (cuts < hi)])
    bounds = np.concatenate(([lo], cuts, [hi]))
    t, w = np.polynomial.legendre.leggauss(n_points)
    half = 0.5 * np.diff(bounds)
    mid = 0.5 * (bounds[:-1] + bounds[1:])
    nodes = mid[:, None] + half[:, None] * t[None, :]
    weights = half[:, None] * w[None, :]
    return nodes.ravel(), weights.ravel()


def _lens_area(r: FloatArray, r_cyl: float, d: FloatArray) -> FloatArray:
    """
    Return the area of the intersection of a disc of radius `r` with a disc
    of radius `r_cyl` whose centers are `d` apart (MathWorld, Circle-Circle
    Intersection, eq. 15), broadcast over `r` and `d`.
    """
    r, d = np.broadcast_arrays(r, d)
    small = np.minimum(r, r_cyl)
    area = np.where(d <= np.abs(r_cyl - r), np.pi * small**2, 0.0)
    lens = (d > np.abs(r_cyl - r)) & (d < r + r_cyl)
    rl, dl = r[lens], d[lens]
    cos_r = np.clip((dl**2 + rl**2 - r_cyl**2) / (2.0 * dl * rl), -1.0, 1.0)
    cos_c = np.clip((dl**2 + r_cyl**2 - rl**2) / (2.0 * dl * r_cyl), -1, 1)
    root = np.sqrt(
        np.maximum(
            (-dl + rl + r_cyl) * (dl + rl - r_cyl) * (dl - rl + r_cyl)
            * (dl + rl + r_cyl),
            0.0,
        )
    )
    area[lens] = (
        rl**2 * np.arccos(cos_r) + r_cyl**2 * np.arccos(cos_c) - 0.5 * root
    )
    return np.asarray(area, dtype=np.float64)


def _ball_in_cylinder(
    r_ball: float, r_cyl: float, d: FloatArray
) -> FloatArray:
    """
    Return the volume of a ball of radius `r_ball` inside the coaxial
    cylinder of radius `r_cyl`, for an array of axial distances `d` of the
    ball's center from the axis.

    The ball is sliced by planes normal to the axis; the slice at height
    ``z`` is a disc of radius ``sqrt(r_ball**2 - z**2)`` and its share is a
    lens area. The integral over ``z`` is a Gauss-Legendre rule on the
    pieces between the heights where the lens changes regime.
    """
    d = np.asarray(d, dtype=np.float64)
    full = 4.0 * np.pi * r_ball**3 / 3.0
    out = np.where(r_cyl >= d + r_ball, full, 0.0)
    cut = (r_cyl > 0.0) & (r_cyl < d + r_ball) & (r_cyl > d - r_ball)
    if not np.any(cut):
        return np.asarray(out, dtype=np.float64)
    dc = d[cut]
    # Heights where the slice disc stops containing the cylinder's disc
    # (z_in) and where it stops being cut by it (z_out).
    z_in = np.sqrt(np.maximum(r_ball**2 - (r_cyl + dc) ** 2, 0.0))
    z_out = np.sqrt(np.maximum(r_ball**2 - (r_cyl - dc) ** 2, 0.0))
    bounds = np.stack(
        [np.zeros_like(dc), z_in, z_out, np.full_like(dc, r_ball)], axis=1
    )
    t, w = np.polynomial.legendre.leggauss(_N_GAUSS_AXIAL)
    half = 0.5 * np.diff(bounds, axis=1)  # (n, 3)
    mid = 0.5 * (bounds[:, :-1] + bounds[:, 1:])
    z = mid[:, :, None] + half[:, :, None] * t  # (n, 3, m)
    radius = np.sqrt(np.maximum(r_ball**2 - z**2, 0.0))
    area = _lens_area(radius, r_cyl, dc[:, None, None])
    out[cut] = 2.0 * np.sum(area * w * half[:, :, None], axis=(1, 2))
    return np.asarray(out, dtype=np.float64)


def sphere_cylinder_intersection(
    r_sphere: float, r_cylinder: float, d: float
) -> float:
    """
    Compute the volume of a ball that lies inside a cylinder whose axis
    passes at distance `d` from the ball's center.

    Parameters
    ----------
    r_sphere : float
        Radius of the ball. Must be positive.
    r_cylinder : float
        Radius of the cylinder; zero or negative gives an empty cylinder.
    d : float
        Distance of the ball's center from the cylinder's axis, ``>= 0``.

    Returns
    -------
    float
        The volume of the ball inside the cylinder.

    Raises
    ------
    ValueError
        If `r_sphere` is not positive or `d` is negative.

    Notes
    -----
    The ball is sliced by planes normal to the axis. At height :math:`z`
    the slice is a disc of radius :math:`\\sqrt{R_b^2 - z^2}` whose center
    is :math:`d` from the axis, and the cylinder cuts the plane in a disc of
    radius :math:`R` about the axis; the slice's share is the lens area of
    the two discs (MathWorld, Circle-Circle Intersection, eq. 15), and the
    volume is its integral over :math:`z`. The integral has no closed form,
    so it is a Gauss-Legendre rule with 32 nodes on each of the three
    pieces between the heights where the lens changes regime; the result
    agrees with an adaptive quadrature to about :math:`10^{-8}` of the
    ball's volume. Time is :math:`O(1)`.

    Two limits are exact: the ball lies inside the cylinder when
    :math:`R \\ge d + R_b`, and outside it when :math:`R \\le d - R_b`. On
    the axis, :math:`d = 0`, the volume inside a cylinder of radius
    :math:`R < R_b` is :math:`\\tfrac{4\\pi}{3}\\,[R_b^3 - (R_b^2 -
    R^2)^{3/2}]`.

    Examples
    --------
    >>> from polyphys.analyze.density import sphere_cylinder_intersection
    >>> round(sphere_cylinder_intersection(1.0, 2.0, 0.5), 10)
    4.1887902048
    >>> sphere_cylinder_intersection(1.0, 2.0, 3.5)
    0.0

    On the axis with a narrower cylinder, the closed form above gives
    :math:`\\tfrac{4\\pi}{3}(1 - 0.8^3)`:

    >>> round(sphere_cylinder_intersection(1.0, 0.6, 0.0), 6)
    2.04413
    >>> round(4 * 3.141592653589793 / 3 * (1 - 0.8**3), 6)
    2.04413

    .. versionadded:: 0.6.0
    """
    radius = _check_positive('r_sphere', r_sphere)
    if d < 0:
        raise ValueError(f'd must be non-negative; got {d}.')
    return float(
        _ball_in_cylinder(radius, float(r_cylinder), np.array([float(d)]))[0]
    )


def _segment_second_antiderivative(u: FloatArray, r_ball: float) -> FloatArray:
    """
    Return :math:`G_2(u) = \\int_{-\\infty}^{u} G`, where :math:`G(u)` is
    the volume of a ball of radius `r_ball` below the plane at signed height
    `u` from its center (the spherical segment from :math:`-R_b` to
    :math:`u`); so :math:`G_2` is the second antiderivative of the ball's
    cross-section area, used to average a slab's share over the bead's
    position in closed form.
    """
    r = r_ball
    full = 4.0 * np.pi * r**3 / 3.0
    uc = np.clip(u, -r, r)
    inside = np.pi * (
        r**2 * (uc + r) ** 2 / 2.0 - (uc**4 - r**4) / 12.0
        - r**3 * (uc + r) / 3.0
    )
    return np.asarray(
        np.where(u >= r, full * u, np.where(u <= -r, 0.0, inside)),
        dtype=np.float64,
    )


_segment = np.vectorize(spherical_segment, otypes=[np.float64])
_sphere_sphere = np.vectorize(sphere_sphere_intersection, otypes=[np.float64])


def _slab_entries(
    edges: FloatArray, r_ball: float, shares: ShareMode
) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp], FloatArray]:
    """
    Return the ``(row, col, value)`` triples of the share matrix of periodic
    slabs: the volume of a test bead in slab `row` that lies in slab `col`,
    with the slabs that wrap around the ends counted through their images.
    """
    n_bins = edges.size - 1
    period = edges[-1] - edges[0]
    if 2.0 * r_ball >= period:
        raise ValueError(
            f'a bead of diameter {2 * r_ball} does not fit the periodic'
            f' length {period}.'
        )
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    reach = int(np.ceil((r_ball + widths.max()) / widths.min())) + 1
    rows_i = np.arange(n_bins)
    rows, cols, vals = [], [], []
    for offset in range(-reach, reach + 1):
        unwrapped = rows_i + offset
        image = np.floor_divide(unwrapped, n_bins)
        col = unwrapped - image * n_bins
        lo = edges[col] + image * period
        hi = edges[col + 1] + image * period
        if shares == 'center':
            value = _segment(r_ball, lo - centers, hi - centers)
        else:
            g2 = _segment_second_antiderivative
            value = (
                g2(hi - edges[:-1], r_ball) - g2(hi - edges[1:], r_ball)
                - g2(lo - edges[:-1], r_ball) + g2(lo - edges[1:], r_ball)
            ) / widths
        keep = value > 0.0
        rows.append(rows_i[keep])
        cols.append(col[keep])
        vals.append(value[keep])
    return (
        np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    )


def _shell_entries(
    edges: FloatArray,
    r_ball: float,
    scheme: _SchemeT,
    shares: ShareMode,
    inside: Callable[[float, FloatArray], FloatArray],
) -> tuple[npt.NDArray[np.intp], npt.NDArray[np.intp], FloatArray]:
    """
    Return the ``(row, col, value)`` triples of the share matrix of
    concentric shells. `inside(edge, x)` is the volume of a bead centered at
    the radial positions `x` that lies within radius `edge`; a shell's share
    is the difference between its two edges.
    """
    n_bins = edges.size - 1
    power = _MEASURE_POWER[scheme]
    centers = 0.5 * (edges[:-1] + edges[1:])
    rows, cols, vals = [], [], []
    for i in range(n_bins):
        lo, hi = edges[i], edges[i + 1]
        j_lo = max(int(np.searchsorted(edges, lo - r_ball, 'right')) - 1, 0)
        j_hi = min(
            int(np.searchsorted(edges, hi + r_ball, 'left')) - 1, n_bins - 1
        )
        shell_edges = edges[j_lo: j_hi + 2]
        if shares == 'center':
            x = np.array([centers[i]])
            weights = np.array([1.0])
            norm = 1.0
        else:
            cuts = np.concatenate(
                (np.abs(shell_edges - r_ball), shell_edges + r_ball)
            )
            x, weights = _gauss_pieces(lo, hi, cuts, _N_GAUSS_POSITION)
            weights = weights * x**power
            norm = (hi ** (power + 1) - lo ** (power + 1)) / (power + 1)
        within = np.stack([inside(edge, x) for edge in shell_edges])
        value = (np.diff(within, axis=0) @ weights) / norm
        keep = value > 0.0
        rows.append(np.full(int(keep.sum()), i))
        cols.append(np.arange(j_lo, j_hi + 1)[keep])
        vals.append(value[keep])
    return (
        np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    )


def share_matrix(
    edges: npt.ArrayLike,
    bead_diameter: float,
    geometry: GeometryT,
    direction: DirectionT,
    *,
    shares: ShareMode = 'average',
    intersection: IntersectionMode = 'exact',
) -> csr_array:
    """
    Build the share matrix: the volume of a test bead centered in each bin
    that lies in each other bin.

    Parameters
    ----------
    edges : array_like
        The bin edges, shape ``(n_bins + 1,)``, finite and strictly
        increasing; non-negative for a radial direction.
    bead_diameter : float
        The bead diameter :math:`a`, positive; along the periodic ``'z'`` it
        must be smaller than the domain length.
    geometry : :py:obj:`~.GeometryT`
        ``'cylindrical'`` or ``'cubic'``.
    direction : :py:obj:`~.DirectionT`
        ``'r'`` for shells, or ``'z'`` for the periodic slabs of a cylinder.
    shares : :py:data:`ShareMode`, default ``'average'``
        ``'average'`` averages the test bead over its bin with the bin's
        measure; ``'center'`` places it at the bin center, as the research
        code did.
    intersection : :py:data:`IntersectionMode`, default ``'exact'``
        The kernel of a cylindrical shell: :func:`sphere_cylinder_intersection`
        or the sphere-sphere volume the research code used. Ignored for the
        other pairs, whose kernels are exact closed forms.

    Returns
    -------
    scipy.sparse.csr_array
        ``S``, shape ``(n_bins, n_bins)``, with ``S[i, j]`` the volume of a
        test bead of bin ``i`` inside bin ``j``; ``S @ rho`` is the volume
        fraction. Each row has about ``a / h + 2`` nonzero entries. Along a
        radial direction the row sum falls short of the bead volume in the
        outermost bins, where the bead reaches past the last edge.

    Raises
    ------
    ValueError
        If `edges` breaks its contract, if the pair is not supported, if
        `bead_diameter` is not positive or does not fit the period, or if
        `shares` or `intersection` is not a known mode.

    Notes
    -----
    The nonzero band is found with :func:`numpy.searchsorted`; the entries
    are evaluated as the module docstring describes. Time is
    :math:`O(n_\\mathrm{bins} \\, k)` kernel evaluations with
    :math:`k \\approx a / h + 2` the band width (times the Gauss nodes for
    ``'average'``, and times the axial nodes for ``'exact'``), and space
    is :math:`O(n_\\mathrm{bins} \\, k)`. The matrix depends on the edges
    and the bead only, so build it once per space and reuse it for every
    whole and every frame.

    With ``shares='average'`` the matrix satisfies
    :math:`V_i S_{ij} = V_j S_{ji}`, which makes :func:`volume_fraction`
    conserve the total bead volume and return :math:`\\rho v` for a uniform
    fluid, both to round-off. With ``shares='center'`` the row sums are
    still the bead volume, but the symmetry, and with it the conservation,
    holds only to :math:`O(h^2)` along :math:`r`. The symmetry also needs
    the kernel to be the true indicator, so with ``intersection='sphere'``
    it holds only to :math:`O(a / R)`.

    Examples
    --------
    Four periodic slabs of width 1 and a bead of diameter 1: at the bin
    center the bead lies wholly in its own slab, so each row is the bead
    volume :math:`\\pi / 6` on the diagonal:

    >>> import numpy as np
    >>> from polyphys.analyze.density import share_matrix
    >>> s = share_matrix([0, 1, 2, 3, 4], 1.0, 'cylindrical', 'z',
    ...                  shares='center')
    >>> np.round(s.toarray(), 4)
    array([[0.5236, 0.    , 0.    , 0.    ],
           [0.    , 0.5236, 0.    , 0.    ],
           [0.    , 0.    , 0.5236, 0.    ],
           [0.    , 0.    , 0.    , 0.5236]])

    Averaged over the bin, the bead reaches into both neighbors, the last
    slab's into the first through the periodic boundary, and every row
    still sums to :math:`\\pi / 6`:

    >>> s = share_matrix([0, 1, 2, 3, 4], 1.0, 'cylindrical', 'z')
    >>> np.allclose(s.sum(axis=1), np.pi / 6)
    True
    >>> bool(s[3, 0] > 0 and s[3, 0] == s[3, 2])
    True

    .. versionadded:: 0.6.0
    """
    scheme = _check_scheme(geometry, direction)
    e = _check_edges(edges, scheme)
    r_ball = 0.5 * _check_positive('bead_diameter', bead_diameter)
    _check_choice('shares', shares, ('average', 'center'))
    _check_choice('intersection', intersection, ('exact', 'sphere'))
    if scheme == 'slab':
        rows, cols, vals = _slab_entries(e, r_ball, shares)
    else:
        if scheme == 'cylinder' and intersection == 'exact':
            def inside(edge: float, x: FloatArray) -> FloatArray:
                return _ball_in_cylinder(r_ball, edge, x)
        else:
            def inside(edge: float, x: FloatArray) -> FloatArray:
                return np.asarray(
                    _sphere_sphere(r_ball, edge, x), dtype=np.float64
                )
        rows, cols, vals = _shell_entries(e, r_ball, scheme, shares, inside)
    n_bins = e.size - 1
    return csr_array(
        coo_array((vals, (rows, cols)), shape=(n_bins, n_bins))
    )


# ---------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------


def number_density(
    counts: npt.ArrayLike,
    volumes: npt.ArrayLike,
    *,
    n_frames: int = 1,
    n_particles: int | None = None,
) -> FloatArray:
    """
    Compute the local number density from bin counts and bin volumes.

    Parameters
    ----------
    counts : array_like
        The bead centers counted per bin, shape ``(n_bins,)`` summed over
        the frames, or ``(n_frames, n_bins)`` with one row per frame (then
        leave `n_frames` at 1). Non-negative.
    volumes : array_like
        The bin volumes, shape ``(n_bins,)``, positive; see
        :func:`bin_volumes`.
    n_frames : int, default 1
        The number of frames the counts were summed over.
    n_particles : int, optional
        The number of beads of this species. When given, each row of
        `counts` must sum to ``n_frames * n_particles``, so that no bead
        was dropped when the positions were binned.

    Returns
    -------
    numpy.ndarray
        :math:`\\rho = H / (F V)`, the same shape as `counts`, in beads per
        cubic length unit.

    Raises
    ------
    ValueError
        If the shapes do not match, a count is negative, a volume is not
        positive, `n_frames` is not positive, or the counts do not sum to
        ``n_frames * n_particles``.

    Notes
    -----
    Time and space are :math:`O(n_\\mathrm{frames} \\, n_\\mathrm{bins})`.

    Examples
    --------
    >>> from polyphys.analyze.density import number_density
    >>> number_density([2, 4], [1.0, 2.0], n_frames=2)
    array([1., 1.])
    >>> number_density([2, 4], [1.0, 2.0], n_frames=2, n_particles=3)
    array([1., 1.])

    .. versionadded:: 0.6.0
    """
    h = np.asarray(counts, dtype=np.float64)
    v = np.asarray(volumes, dtype=np.float64)
    if h.ndim not in (1, 2):
        raise ValueError(
            f'counts must be one- or two-dimensional; got {h.ndim}'
            ' dimension(s).'
        )
    if v.ndim != 1 or v.shape[0] != h.shape[-1]:
        raise ValueError(
            f'volumes must have shape ({h.shape[-1]},) to match counts;'
            f' got {v.shape}.'
        )
    if np.any(h < 0):
        raise ValueError('counts must be non-negative.')
    if not np.all(v > 0):
        raise ValueError('volumes must be positive.')
    if n_frames < 1:
        raise ValueError(f'n_frames must be positive; got {n_frames}.')
    if n_particles is not None:
        expected = n_frames * n_particles
        totals = h.sum(axis=-1)
        if np.any(totals != expected):
            raise ValueError(
                f'counts sum to {totals.tolist()}, not to n_frames *'
                f' n_particles = {expected}; a bead was lost or double'
                ' counted when the positions were binned.'
            )
    return np.asarray(h / (n_frames * v), dtype=np.float64)


def volume_fraction(rho: npt.ArrayLike, share: csr_array) -> FloatArray:
    """
    Compute the local volume fraction from the number density and the share
    matrix.

    Parameters
    ----------
    rho : array_like
        The local number density, shape ``(n_bins,)`` or ``(n_frames,
        n_bins)``; see :func:`number_density`.
    share : scipy.sparse.csr_array
        The share matrix of the same bins and the bead's diameter; see
        :func:`share_matrix`.

    Returns
    -------
    numpy.ndarray
        :math:`\\phi_i = \\sum_j S_{ij} \\rho_j`, the same shape as `rho`,
        dimensionless.

    Raises
    ------
    ValueError
        If `rho` is not one- or two-dimensional, or its last axis does not
        match `share`.

    Notes
    -----
    One sparse matrix-vector product per row of `rho`: time
    :math:`O(n_\\mathrm{frames} \\, n_\\mathrm{bins} \\, k)` with :math:`k`
    the band width, space :math:`O(n_\\mathrm{frames} \\, n_\\mathrm{bins})`.

    Examples
    --------
    A uniform density through the averaged share matrix returns
    :math:`\\rho v` in every slab, here :math:`(3 / \\pi) (\\pi / 6) = 1/2`:

    >>> import numpy as np
    >>> from polyphys.analyze.density import share_matrix, volume_fraction
    >>> s = share_matrix([0, 1, 2, 3, 4], 1.0, 'cylindrical', 'z')
    >>> np.round(volume_fraction(np.full(4, 3 / np.pi), s), 12)
    array([0.5, 0.5, 0.5, 0.5])

    .. versionadded:: 0.6.0
    """
    density = np.asarray(rho, dtype=np.float64)
    if density.ndim not in (1, 2):
        raise ValueError(
            f'rho must be one- or two-dimensional; got {density.ndim}'
            ' dimension(s).'
        )
    if density.shape[-1] != share.shape[1]:
        raise ValueError(
            f'rho has {density.shape[-1]} bins but share has'
            f' {share.shape[1]}.'
        )
    if density.ndim == 1:
        return np.asarray(share @ density, dtype=np.float64)
    return np.asarray((share @ density.T).T, dtype=np.float64)


def normalize_profile(
    profile: npt.ArrayLike,
    mode: NormalizeMode,
    volumes: npt.ArrayLike | None = None,
) -> FloatArray:
    """
    Normalize a profile along its bins.

    Parameters
    ----------
    profile : array_like
        A profile, shape ``(n_bins,)`` or ``(n_frames, n_bins)``.
    mode : :py:data:`NormalizeMode`
        ``'none'`` returns a copy; ``'sum'`` divides each row by its sum
        over the bins, as the research code did; ``'integral'`` divides
        each row by :math:`\\sum_i p_i V_i`, so that the result integrates
        to 1 against the bin volumes.
    volumes : array_like, optional
        The bin volumes, shape ``(n_bins,)``; needed for ``'integral'``.

    Returns
    -------
    numpy.ndarray
        The normalized profile, the same shape as `profile`.

    Raises
    ------
    ValueError
        If `mode` is unknown, if ``'integral'`` is asked for without
        `volumes` or with volumes of the wrong shape, or if a row's
        normalizer is zero.

    Notes
    -----
    A normalized profile must not be passed to
    :mod:`polyphys.analyze.sumrule`, whose species sum needs the physical
    scale of every profile. Time and space are
    :math:`O(n_\\mathrm{frames} \\, n_\\mathrm{bins})`.

    Examples
    --------
    >>> from polyphys.analyze.density import normalize_profile
    >>> normalize_profile([1.0, 3.0], 'sum')
    array([0.25, 0.75])
    >>> normalize_profile([1.0, 3.0], 'integral', volumes=[2.0, 2.0])
    array([0.125, 0.375])

    .. versionadded:: 0.6.0
    """
    _check_choice('mode', mode, ('none', 'sum', 'integral'))
    p = np.array(profile, dtype=np.float64)
    if p.ndim not in (1, 2):
        raise ValueError(
            f'profile must be one- or two-dimensional; got {p.ndim}'
            ' dimension(s).'
        )
    if mode == 'none':
        return p
    if mode == 'sum':
        normalizer = p.sum(axis=-1, keepdims=True)
    else:
        if volumes is None:
            raise ValueError("volumes are needed for mode 'integral'.")
        v = np.asarray(volumes, dtype=np.float64)
        if v.shape != (p.shape[-1],):
            raise ValueError(
                f'volumes must have shape ({p.shape[-1]},); got {v.shape}.'
            )
        normalizer = (p * v).sum(axis=-1, keepdims=True)
    if np.any(normalizer == 0):
        raise ValueError('a profile with a zero normalizer cannot be'
                         ' normalized.')
    return np.asarray(p / normalizer, dtype=np.float64)


@dataclass(frozen=True)
class Profiles:
    """
    The local number density and volume fraction of one species along one
    coordinate.

    .. versionadded:: 0.6.0

    Attributes
    ----------
    rho : numpy.ndarray
        The local number density, shape ``(n_bins,)`` or ``(n_frames,
        n_bins)``, in beads per cubic length unit, normalized as asked.
    phi : numpy.ndarray
        The local volume fraction, the same shape, dimensionless,
        normalized as asked.
    volumes : numpy.ndarray
        The bin volumes, shape ``(n_bins,)``.
    """

    rho: FloatArray
    phi: FloatArray
    volumes: FloatArray


def density_profiles(
    counts: npt.ArrayLike,
    edges: npt.ArrayLike,
    bead_diameter: float,
    geometry: GeometryT,
    direction: DirectionT,
    *,
    n_frames: int = 1,
    n_particles: int | None = None,
    lcyl: float | None = None,
    dcyl: float | None = None,
    lcube: float | None = None,
    shares: ShareMode = 'average',
    intersection: IntersectionMode = 'exact',
    normalize: NormalizeMode = 'none',
) -> Profiles:
    """
    Compute the local number density and volume fraction of one species in
    one call.

    Parameters
    ----------
    counts : array_like
        The bead centers counted per bin, shape ``(n_bins,)`` summed over
        the frames or ``(n_frames, n_bins)``; see :func:`number_density`.
    edges : array_like
        The bin edges, shape ``(n_bins + 1,)``.
    bead_diameter : float
        The bead diameter :math:`a` of the species.
    geometry : :py:obj:`~.GeometryT`
        ``'cylindrical'`` or ``'cubic'``.
    direction : :py:obj:`~.DirectionT`
        ``'r'`` or ``'z'``.
    n_frames : int, default 1
        The number of frames summed into `counts`.
    n_particles : int, optional
        The number of beads of the species, to check that none was lost.
    lcyl, dcyl, lcube : float, optional
        The cylinder length and diameter and the cube side, as
        :func:`bin_volumes` needs them. Along ``'z'`` the slabs take the
        effective diameter ``dcyl - bead_diameter``, the diameter of the
        region a bead center can reach, so that a uniform fluid returns the
        bulk volume fraction as the parsers define it.
    shares : :py:data:`ShareMode`, default ``'average'``
        See :func:`share_matrix`.
    intersection : :py:data:`IntersectionMode`, default ``'exact'``
        See :func:`share_matrix`.
    normalize : :py:data:`NormalizeMode`, default ``'none'``
        Applied to both profiles; see :func:`normalize_profile`. Leave it
        at ``'none'`` for profiles that feed :mod:`polyphys.analyze.sumrule`.

    Returns
    -------
    Profiles
        The two profiles and the bin volumes.

    Raises
    ------
    ValueError
        As the four functions it calls raise, and if ``dcyl <=
        bead_diameter`` along ``'z'``.

    Notes
    -----
    This builds the share matrix for one call. For many wholes or frames on
    the same bins, call :func:`share_matrix` once and
    :func:`number_density` and :func:`volume_fraction` per profile. Time is
    that of :func:`share_matrix` plus :math:`O(n_\\mathrm{frames} \\,
    n_\\mathrm{bins} \\, k)`.

    Examples
    --------
    Three beads of diameter 1 per slab in a cylinder of effective diameter
    2 (so ``dcyl = 3``): the slab area is :math:`\\pi`, the density
    :math:`3 / \\pi`, and the volume fraction :math:`(3 / \\pi)(\\pi / 6) =
    1 / 2`:

    >>> import numpy as np
    >>> from polyphys.analyze.density import density_profiles
    >>> out = density_profiles([3, 3, 3, 3], [0, 1, 2, 3, 4], 1.0,
    ...                        'cylindrical', 'z', dcyl=3.0, n_particles=12)
    >>> np.round(out.rho, 6)
    array([0.95493, 0.95493, 0.95493, 0.95493])
    >>> np.round(out.phi, 12)
    array([0.5, 0.5, 0.5, 0.5])

    .. versionadded:: 0.6.0
    """
    scheme = _check_scheme(geometry, direction)
    a = _check_positive('bead_diameter', bead_diameter)
    slab_diameter = dcyl
    if scheme == 'slab':
        diameter = _check_positive('dcyl', dcyl)
        if diameter <= a:
            raise ValueError(
                f'dcyl = {diameter} must exceed the bead diameter {a}.'
            )
        slab_diameter = diameter - a
    volumes = bin_volumes(
        edges, geometry, direction, lcyl=lcyl, dcyl=slab_diameter,
        lcube=lcube,
    )
    share = share_matrix(
        edges, a, geometry, direction, shares=shares,
        intersection=intersection,
    )
    rho = number_density(
        counts, volumes, n_frames=n_frames, n_particles=n_particles
    )
    phi = volume_fraction(rho, share)
    return Profiles(
        rho=normalize_profile(rho, normalize, volumes),
        phi=normalize_profile(phi, normalize, volumes),
        volumes=volumes,
    )
