"""
==========================================================
:mod:`polyphys.analyze.measurer`
==========================================================


The :mod:`polyphys.analyze.measurer` module provides a suite of functions /
helper functions for computing geometric, structural, statistical, and
thermodynamic properties of a particle or a group of particles.

Contents
========

Core geometric and structural measurements:
:func:`apply_pbc_orthogonal`, :func:`pair_distance`, :func:`end_to_end`,
:func:`transverse_size`, :func:`max_distance`, and :func:`fsd`.

Summary statistics:
:func:`simple_stats` and :func:`sem`.

Advanced geometric calculations:
:func:`sphere_sphere_intersection` and :func:`spherical_segment`.

Binning and histogram processing:
:func:`create_bin_edge_and_hist`, :func:`fixedsize_bins`,
:func:`radial_histogram`, :func:`radial_cyl_histogram`,
:func:`axial_histogram`, :func:`azimuth_cyl_histogram`, and
:func:`planar_cartesian_histogram`.

Number-density and volume-fraction conversions live in
:mod:`polyphys.manage.utils`.

References
==========
For Feret's statistical diameter:

- Wang Y, Teraoka I, Hansen FY, Peters GH, Ole H. "A Theoretical Study of
  the Separation Principle in Size Exclusion Chromatography." Macromolecules
  2010, 43, 3, 1651-1659. https://doi.org/10.1021/ma902377g

For spherical segment calculations:

- Weisstein, Eric W. "Spherical Segment." From MathWorld--A Wolfram Web
  Resource. https://mathworld.wolfram.com/SphericalSegment.html

For sphere-sphere intersection:

- Weisstein, Eric W. "Sphere-Sphere Intersection." From MathWorld--A Wolfram
  Web Resource. https://mathworld.wolfram.com/Sphere-SphereIntersection.html
"""

from typing import Dict, Tuple, Optional, Union, Any, List
import numpy as np

from ..manage.types import AxisT, EntityT, PropertyT, BinT
from ..manage.utils import invalid_keyword

# Relative tolerance for deciding whether a range holds a whole number of
# bins, used by :func:`create_bin_edge_and_hist` and :func:`fixedsize_bins`.
# The ratio ``(lmax - lmin) / bin_size`` is computed in float64, so a range
# that is exactly a whole number of bins can come out a few ulps above or
# below that integer (``2.1 / 0.3`` is ``7.000000000000001``). A plain
# ``ceil`` then adds a bin, and a plain comparison can leave the edges a few
# 1e-15 inside the range. The tolerance has to sit between two scales:
#
# - the accumulated round-off in the ratio, a few times the float64 epsilon
#   (2.2e-16), so below about 1e-15 relative;
# - the smallest deviation a user can intend, which is set by the precision of
#   the inputs: box sizes and bead diameters are read from filenames with at
#   most seven significant digits, so about 1e-7 relative.
#
# Any value between those classifies every realistic input correctly. 1e-9
# keeps six orders of magnitude from the round-off side and two from the
# intent side; the margins are unequal on purpose, because the two errors are
# not equally bad. Calling a true multiple a non-multiple adds a bin and, for
# 'ordinary' bins, shifts the whole grid by half a bin; calling an intended
# non-multiple a multiple moves the edges by less than
# ``1e-9 * (lmax - lmin)``.
_BIN_COUNT_RTOL = 1e-9


def _count_bins(length: float, bin_size: float) -> Tuple[int, float]:
    """
    Return the number of bins of width `bin_size` that cover `length`, and
    the excess length they add.

    When ``length / bin_size`` is a whole number up to the relative tolerance
    :data:`_BIN_COUNT_RTOL`, the count is that whole number and the excess is
    exactly ``0.0``; otherwise the count is the ceiling of the ratio and the
    excess ``n_bins * bin_size - length`` is strictly positive, by a margin
    far above float64 round-off.

    Parameters
    ----------
    length : float
        Length of the range to cover, positive.
    bin_size : float
        Width of each bin, positive.

    Returns
    -------
    Tuple[int, float]
        ``(n_bins, excess)``.

    Examples
    --------
    >>> _count_bins(2.1, 0.3)
    (7, 0.0)
    >>> n_bins, excess = _count_bins(2.0, 0.7)
    >>> n_bins, round(excess, 12)
    (3, 0.1)
    """
    ratio = length / bin_size
    n_bins = int(np.rint(ratio))
    if n_bins > 0 and np.isclose(
        ratio, n_bins, rtol=_BIN_COUNT_RTOL, atol=0.0
    ):
        return n_bins, 0.0
    n_bins = int(np.ceil(ratio))
    return n_bins, n_bins * bin_size - length


def _histogram_or_raise(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """
    Histogram `values` over `edges`, raising if any value falls outside them.

    `np.histogram` silently drops values below the first edge or above the
    last one. For the spatial profiles this module serves, every position
    belongs to the system, so a dropped position is a defect in the caller's
    edges or an unwrapped coordinate, not data to ignore.

    Parameters
    ----------
    values : np.ndarray
        One-dimensional array of values to bin.
    edges : np.ndarray
        Monotonically increasing bin edges, including the rightmost edge.

    Returns
    -------
    np.ndarray
        Counts per bin, shape ``(len(edges) - 1,)``; their sum is
        ``values.size``.

    Raises
    ------
    ValueError
        If any value lies outside ``[edges[0], edges[-1]]``.

    Examples
    --------
    >>> edges = np.array([0., 1., 2., 3.])
    >>> _histogram_or_raise(np.array([0.5, 2.0, 3.0]), edges)
    array([1, 0, 2])
    >>> _histogram_or_raise(np.array([0.5, 9.0]), edges)
    Traceback (most recent call last):
        ...
    ValueError: 1 of 2 positions fall outside the edges [0.0, 3.0].
    """
    hist, _ = np.histogram(values, bins=edges)
    n_missed = values.size - int(hist.sum())
    if n_missed:
        raise ValueError(
            f"{n_missed} of {values.size} positions fall outside the edges "
            f"[{edges[0]}, {edges[-1]}]."
        )
    return hist


def apply_pbc_orthogonal(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: Dict[AxisT, float],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Updates the periodic boundary condition (PBC) lengths and their inverses
    based on specified axes in `pbc`.

    Parameters
    ----------
    pbc_lengths : numpy.ndarray
        Array of box lengths along each axis (typically shape (3,)).
    pbc_lengths_inverse : numpy.ndarray
        Array of inverse box lengths along each axis (shape (3,)).
    pbc : Dict[:py:data:`AxisT`, float]
        Dictionary with axis indices (0=x, 1=y, 2=z) as keys and box lengths
        as values. See :mod:`polyphys.manage.types`.

    Return
    ------
    pbc_lengths : numpy.ndarray
        The updated array of lengths in each direction.
    pbc_lengths_inverse : numpy.ndarray
        The updated array of the inverse lengths in each direction.

    Raises
    ------
    ZeroDivisionError
        If any specified length is zero, which would result in division by
        zero.
    ValueError
        If any specified length is negative, which is physically invalid for
        box dimensions.

    Examples
    --------
    >>> pbc_lengths = np.zeros(3)
    >>> pbc_lengths_inverse = np.zeros(3)
    >>> pbc = {0: 10.0, 2: 20.0}
    >>> apply_pbc_orthogonal(pbc_lengths, pbc_lengths_inverse, pbc)
    (array([10.,  0., 20.]), array([0.1 , 0.  , 0.05]))
    """
    for ax, length in pbc.items():
        if length == 0.0:
            raise ZeroDivisionError(f"Length for axis {ax} cannot be zero.")
        if length < 0.0:
            raise ValueError(
                f"Length for axis {ax} must be positive; got {length}."
            )
        pbc_lengths[ax] = length
        pbc_lengths_inverse[ax] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def pair_distance(
    positions: np.ndarray, pbc: Optional[Dict[AxisT, float]] = None
) -> np.ndarray:
    """
    Compute the pairwise distance vector between two particles in Cartesian
    coordinates, optionally applying the minimum image convention (PBC).

    Parameters
    ----------
    positions : numpy.ndarray
        Array of shape (2, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    pbc : Optional[Dict[:py:data:`AxisT`, float]], default None
        Dictionary with axis indices (0=x, 1=y, 2=z) as keys and box lengths as
        values. If provided, periodic boundary conditions (PBC) are applied;
        see :mod:`polyphys.manage.types`.

    Returns
    -------
    np.ndarray
        A 1D array with pair distances along each axis.

    Raises
    ------
    ValueError
        If `positions` does not contain exactly two sets of coordinates.

    Notes
    -----
    This function returns the component-wise distance vector, not the scalar
    Euclidean norm.

    Examples
    --------
    >>> positions = np.array([[1.0, 2.0, 3.0], [4.0, 2.0, 3.0]])
    >>> pair_distance(positions)
    array([3., 0., 0.])

    >>> pbc = {0: 8.0}
    >>> positions = np.array([[1.0, 2.0, 3.0], [9.5, 2.0, 3.0]])
    >>> pair_distance(positions, pbc)
    array([0.5, 0. , 0. ])
    """
    n_atoms, n_dims = positions.shape
    if n_atoms != 2:
        raise ValueError("'pair_distance' only works for two atoms.")
    # Calculation in the center of geometry (cog) of the atom group.
    delta = positions[1] - positions[0]

    if pbc is not None:
        n_dims = positions.shape[1]
        # Prepare length and inverse-length arrays
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
        )
        # Apply Minimum Image Convention
        delta -= pbc_lengths * np.round(delta * pbc_lengths_inverse)

    return delta


def end_to_end(positions: np.ndarray) -> Union[np.floating, np.ndarray]:
    """
    Computes the end-to-end distance of a linear polymer, in the frame of
    reference located at the polymer's center of geometry.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atomic positions, and sorted
        by atom number form 1 to N.

    Returns
    -------
    float
        End-to-end distance of the polymer.

    Raises
    ------
    ValueError
        If `positions` contains fewer than two atoms.

    Notes
    -----
    The distance is calculated in a reference frame centered at the center of
    geometry of atoms.

    Examples
    --------
    >>> positions = np.array([[0.0, 0.0, 0.0],
    ...                       [1.0, 0.0, 0.0],
    ...                       [2.0, 0.0, 0.0]])
    >>> end_to_end(positions)
    np.float64(2.0)

    >>> positions = np.array([[1.0, 1.0, 1.0],
    ...                       [2.0, 2.0, 2.0],
    ...                       [3.0, 3.0, 3.0]])
    >>> end_to_end(positions)
    np.float64(3.4641016151377544)
    """
    # calculation in the center of geometry of the atom group.
    if positions.shape[0] < 2:
        raise ValueError("`positions` must contain at least two atoms.")
    centered_positions = positions - np.mean(positions, axis=0)
    return np.linalg.norm(centered_positions[-1] - centered_positions[0])


def transverse_size(positions: np.ndarray, axis: AxisT) -> float:
    """
    Computes the mean transverse size of a group of atoms in the plane
    perpendicular to `axis`.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atomic positions, and sorted
        by atom number form 1 to N.
    axis: :py:data:`AxisT`
        The axis (0 for x, 1 for y, and 2 for z) in the plane perpendicular to
        which the transverse size is calculated. See
        :mod:`polyphys.manage.types`.

    Returns
    -------
    float
        Twice the mean transverse distance (a.k.a., diameter) in the plane
        perpendicular to `axis`.

    Raises
    ------
    IndexError
        If `axis` is out of bounds for `positions`.

    Notes
    -----
    The distance is calculated in a reference frame centered at the center of
    geometry of atoms in the plane perpendicular to a given axis.

    Examples
    --------
    >>> positions = np.array([[1, 0, 0], [1, 2, 2]])
    >>> transverse_size(positions, axis=0)
    np.float64(2.8284271247461903)

    >>> transverse_size(positions, axis=1)
    np.float64(2.0)
    """
    trans_axes = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
    if axis >= positions.shape[1]:
        raise IndexError(
            "Axis {axis} is out of bounds for positions array"
            f"with shape {positions.shape}."
        )
    trans_pos = positions[:, trans_axes[axis]] - np.mean(
        positions[:, trans_axes[axis]], axis=0
    )
    # Times by 2, so we have the diameter not radius:
    return 2 * np.mean(np.linalg.norm(trans_pos, axis=1))


def max_distance(positions: np.ndarray) -> np.ndarray:
    """
    Computes the maximum extent in each Cartesian axis.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.

    Returns
    -------
    np.ndarray
        Maximum extent in each dimension.

    Raises
    ------
    ValueError
        If `positions` is empty.

    Notes
    -----
    The distance is calculated in a reference frame centered at the center of
    geometry of atoms.

    Examples
    --------
    >>> positions = np.array([[0, 0, 0], [4, 2, 6], [2, 3, 0]])
    >>> max_distance(positions)
    array([4., 3., 6.])
    """
    # calculation in the center of geometry of the atom group.
    centered_positions = positions - np.mean(positions, axis=0)
    return np.ptp(centered_positions, axis=0)


def fsd(positions: np.ndarray, axis: AxisT) -> float:
    """
    Computes the mean Feret's statistical diameter (FSD) along a specified
    axis.

    fsd stands for Feret's statistical diameter: other names are the mean
    span dimension, the farthermost distance, or the mean caliper diameter.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.

    axis: :py:data:`AxisT`
        Axis index (0=x, 1=y, 2=z) along which to compute FSD; see
        :mod:`polyphys.manage.types`.

    Returns
    -------
    float
        Mean Feret's diameter along the specified axis.

    Raises
    ------
    IndexError
        If `axis` is out of bounds for `positions`.

    References
    ----------
    "A Theoretical Study of the Separation Principle in Size Exclusion
    Chromatography", Wang Y Teraoka
    I Hansen FY Peters GH Ole H. Macromolecules 2010, 43, 3, 1651-1659
    https://doi.org/10.1021/ma902377g

    Examples
    --------
    >>> positions = np.array([[1, 2, 3], [4, 8, 6]])
    >>> fsd(positions, axis=1)
    np.float64(6.0)
    """
    if axis >= positions.shape[1]:
        raise IndexError(
            "Axis {axis} is out of bounds for positions array"
            f"with shape {positions.shape}."
        )
    # calculation in the center of geometry of the atom group:
    positions = positions - np.mean(positions, axis=0)
    return np.ptp(positions[:, axis])


def simple_stats(
        entity: EntityT,
        data: np.ndarray
) -> Dict[PropertyT, float | np.floating]:
    """
    Compute basic statistics (mean, variance, SEM) for for a physical `entity`.

    Parameters
    ----------
    entity : EntityT
        Name of the physical entity.

    data: np.ndarray
        Array of entity values.

    Returns
    -------
    Dict[PropertyT, float]
        Dictionary with keys '<entity>_mean', '<entity>_var', and
        '<entity>_sem'.

    Raises
    ------
    ValueError
        If `data` is empty.

    Examples
    --------
    >>> import numpy as np
    >>> simple_stats('energy', np.array([1.0, 2.0, 3.0]))  # noqa: E501
    {'energy_mean': np.float64(2.0), 'energy_var': np.float64(1.0), 'energy_sem': np.float64(0.5773502691896258)}
    """
    if len(data) == 0:
        raise ValueError('Input array must not be empty.')

    # Unbiased std, var, and sem.
    return {
        entity + '_mean': np.mean(data),
        entity + '_var': np.var(data, ddof=1),
        entity + '_sem': np.std(data, ddof=1) / np.sqrt(len(data)),
    }


def sem(data: np.ndarray) -> float:
    """
    Calculates the standard error of the mean (SEM).

    Parameters
    ----------
    data: np.ndarray
        1D array of sample data.

    Returns
    -------
    float
        Standard error of the mean.

    Raises
    ------
    ValueError
        If `data` is empty.

    Notes
    -----
    SEM is calculated using sample standard deviation (ddof=1) and sample size.
    The estimate assumes independent samples, such as the means of independent
    replicates. For a time series of correlated frames it understates the
    error; use :func:`polyphys.analyze.stats.sem_blocked` instead.

    Examples
    --------
    >>> import numpy as np
    >>> sem(np.array([1.0, 2.0, 3.0]))
    np.float64(0.5773502691896258)
    """
    if len(data) == 0:
        raise ValueError('Input data array must not be empty.')

    return np.std(data, ddof=1) / len(data) ** 0.5


def spherical_segment(r: float, a: float, b: float) -> float:
    """
    Compute the volume of a spherical segment defined by two parallel planes.

    Parameters
    ----------
    r : float
        Radius of the sphere. Must be positive.
    a : float
        Distance of the first plane from the center of the sphere, along the
        axis of symmetry.
    b : float
        Distance of the second plane from the center of the sphere, along the
        axis of symmetry.

    Returns
    -------
    vol : float
        Volume of the spherical segment.

    Raises
    ------
    ValueError
        If `r` is not a positive number.

    Notes
    -----
    - `a` and `b` can be positive or negative values, as long as they fall
      within the range `[-r, r]`.
    - The function will adjust `a` and `b` to `-r` or `r` if they exceed these
      bounds.
    - If `a = r` or `b = r`, the spherical segment becomes a spherical cap.
    - If both `a` and `b` lie outside the range `[-r, r]` and share the same
      sign, the volume is zero.

    References
    ----------
    Weisstein, Eric W. "Spherical Segment." From MathWorld--A Wolfram Web
    Resource. https://mathworld.wolfram.com/SphericalSegment.html

    Examples
    --------
    >>> spherical_segment(3, 1, 2)
    20.94395102393195
    >>> spherical_segment(3, -3, 3)
    113.09733552923255
    """

    if r <= 0:
        raise ValueError(f"The radius 'r' must be positive. Got {r}.")

    # Ensure the bounds are within [-r, r]
    lower = max(min(a, b), -r)
    upper = min(max(a, b), r)

    # If both bounds are outside [-r, r] with the same sign, the volume is zero
    if lower * upper >= r**2:
        return 0.0
    # Calculate the volume of the spherical segment
    vol = np.pi * (r**2 * (upper - lower) - (upper**3 - lower**3) / 3)
    return vol


def sphere_sphere_intersection(r1: float, r2: float, d: float) -> float:
    """
    Compute the volume of intersection between two spheres.

    The sphere with radius `r1` is separated from the sphere with radius `r2`
    by a distance `d` along the x-axis. Thus, the vector form of the distance
    between their centers is `(d, 0, 0)` in Cartesian coordinates.

    Parameters
    ----------
    r1 : float
        Radius of the first sphere.
    r2 : float
        Radius of the second sphere.
    d : float
        Distance between the centers of the two spheres along the x-axis.

    Returns
    -------
    vol : float
        Volume of the intersection between the two spheres.

    References
    ----------
    Weisstein, Eric W. "Sphere-Sphere Intersection." From MathWorld--A Wolfram
    Web Resource. https://mathworld.wolfram.com/Sphere-SphereIntersection.html

    Examples
    --------
    >>> sphere_sphere_intersection(3, 4, 2)
    94.90227807719167
    >>> sphere_sphere_intersection(3, 4, 10)
    0.0
    """

    r_max = max(r1, r2)
    r_min = min(r1, r2)

    # Volume is zero if one sphere has a radius of zero or no intersection
    # occurs:
    if r1 == 0.0 or r2 == 0.0:
        return 0.0
    if d >= r_min + r_max:
        return 0.0
    if d <= r_max - r_min:
        # The smaller sphere is entirely contained within the larger one
        return 4 * np.pi * r_min**3 / 3

    # Calculate the intersection volume for overlapping spheres
    vol = (
        (np.pi / (12 * d))
        * (r_max + r_min - d) ** 2
        * (
            d**2
            + 2 * d * (r_max + r_min)
            - 3 * (r_min**2 + r_max**2)
            + 6 * r_min * r_max
        )
    )

    return vol


def create_bin_edge_and_hist(
    bin_size: float, lmin: float, lmax: float, output: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create bin edges and an empty histogram for data processing.

    The edges start at `lmin` and step by `bin_size` up to the first edge at
    or beyond `lmax`, so there are ``n_bins = ceil((lmax - lmin) / bin_size)``
    bins of width `bin_size`. When the range is a whole number of bins up to
    the relative tolerance :data:`_BIN_COUNT_RTOL` (1e-9), the last edge is
    exactly `lmax`; otherwise it exceeds `lmax` by less than one bin.

    Parameters
    ----------
    bin_size : float
        Size of each bin.
    lmin : float
        Lower bound of the system in the direction of interest.
    lmax : float
        Upper bound of the system in the direction of interest.
    output : Optional[str], default None
        Filename (including filepath) to which bin edges array is saved. A
        `.npy` extension will be appended to the filename if it does not
        already have one.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        - bin_edges : Array of bin edges, shape ``(n_bins + 1,)``, equal to
          ``lmin + bin_size * k`` for ``k = 0, ..., n_bins``.
        - hist : Zero-initialized histogram of dtype ``np.int64`` and shape
          ``(n_bins,)``. The 64-bit integer type keeps per-frame counts
          accumulated in place (``hist += np.histogram(...)[0]``) from
          overflowing.

    Raises
    ------
    ValueError
        If `lmin` is not less than `lmax`.

    Notes
    -----
    The `.npy` extension is handled internally by `np.save`.

    The edges come from `np.linspace`, so their number is exactly
    ``n_bins + 1``; a float-step `np.arange` can return one edge too many
    when `bin_size` is not exactly representable (for example 0.3). The bin
    count treats a ratio within 1e-9 of a whole number as that whole number;
    see :data:`_BIN_COUNT_RTOL` for why.

    Time and space are O(`n_bins`).

    Examples
    --------
    >>> edges, hist = \
        create_bin_edge_and_hist(bin_size=1.0, lmin=0.0, lmax=5.0)
    >>> edges
    array([0., 1., 2., 3., 4., 5.])
    >>> hist
    array([0, 0, 0, 0, 0])

    A step that is not exactly representable still gives one edge per bin
    plus one, ending exactly at `lmax`:

    >>> edges, hist = create_bin_edge_and_hist(0.3, -6.0, 6.0)
    >>> print(len(edges), len(hist), edges[0], edges[-1])
    41 40 -6.0 6.0

    A range that is not a multiple of `bin_size` is covered by one more bin:

    >>> create_bin_edge_and_hist(0.7, 0.0, 2.0)[0]
    array([0. , 0.7, 1.4, 2.1])
    """
    if lmin >= lmax:
        raise ValueError("'lmin' must be less than 'lmax'.")

    n_bins, excess = _count_bins(lmax - lmin, bin_size)
    upper = lmax if excess == 0.0 else lmin + n_bins * bin_size
    bin_edges = np.linspace(lmin, upper, n_bins + 1)
    hist = np.zeros(n_bins, dtype=np.int64)

    if output is not None:
        np.save(output, bin_edges)

    return bin_edges, hist


def fixedsize_bins(
    bin_size: float,
    lmin: float,
    lmax: float,
    bin_type: BinT = 'ordinary',
    save_bin_edges: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate bin edges and an empty collectors with consistent bin width.

    Every bin type returns ``n_bins + 1`` equally spaced edges of width
    `bin_size` that cover [`lmin`, `lmax`].

    Parameters
    ----------
    bin_size : float
        Width of each bin.
    lmin : float
        Start of the range
    lmax : float
        End of the range
    bin_type : :py:data:`BinT`, default 'ordinary'
        Type of binning. With ``n_bins = ceil((lmax - lmin) / bin_size)``
        bins, the excess length ``n_bins * bin_size - (lmax - lmin)`` is
        zero when the range is a whole number of bins and is otherwise
        shared out as follows:

        - 'ordinary': Extends `lmin` and `lmax` symmetrically, by half the
          excess each. Examples are Cartesian coordinates and spherical polar
          coordinate.
        - 'nonnegative': Extends `lmin` and `lmax` symmetrically if the
          extended `lmin` stays nonnegative; otherwise `lmin` is set to 0 and
          `lmax` to ``n_bins * bin_size``. Examples are radial directions in
          the polar, spherical, and cylindrical coordinate systems.
        - 'periodic': Periodic binning over one period, ``lmax - lmin``,
          which must be a whole number of bins; the edges run from `lmin` to
          `lmax` exactly. Examples are azimuthal directions in cylindrical
          and spherical coordinate systems.

        See the Notes for how "a whole number of bins" is decided.

        See :mod:`polyphys.manage.types`.
    save_bin_edges : Optional[str], default None
        Filename (including filepath) to which bin edges array is saved. A
        `.npy` extension will be appended to the filename if it does not
        already have one.

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:

        - 'n_bins' : int
          Number of bins.
        - 'bin_edges' : np.ndarray
          A monotonically increasing array of bin edges, including the
          rightmost edge.
        - 'collector' : np.ndarray
          Zero-initialized array of dtype ``np.int64`` and shape
          ``(n_bins,)`` for histogram counts. The 64-bit integer type keeps
          per-frame counts accumulated in place
          (``collector += np.histogram(...)[0]``) from overflowing.
        - 'collector_std' : np.ndarray
          Zero-initialized array of dtype ``np.float64`` and shape
          ``(n_bins,)`` for standard deviation values or the squared counts
          they are computed from.
        - 'range' : Tuple[float, float]
            Updated range of bins (`lmin`, `lmax`).

    Raises
    ------
    ValueError
        If `lmin` is not less than `lmax`.
    ValueError
        If `bin_type` is 'nonnegative' and `lmin` is negative: clipping the
        lower end to zero would leave the range below zero uncovered.
    ValueError
        If `bin_type` is 'periodic' and ``lmax - lmin`` is not a multiple of
        `bin_size`. A last bin running past the period would cover values
        that never occur, so a density normalized by its width would be
        biased low.

    Notes
    -----
    The `.npy` extension is handled internally by `np.save`.

    The edges come from `np.linspace`, so their number is exactly
    ``n_bins + 1``, and the first and last edges enclose `lmin` and `lmax`
    exactly: a position at either end of the range falls inside the bins,
    and `np.histogram` counts a value equal to the last edge in the last bin.

    **Deciding whether the range is a whole number of bins.** The ratio
    ``(lmax - lmin) / bin_size`` is computed in float64, so a range that is
    exactly a whole number of bins can come out a few ulps above or below
    that integer: ``2.1 / 0.3`` is ``7.000000000000001``, and
    ``530 * 0.03`` is a few 1e-15 below ``15.9``. Taking the ceiling of such
    a ratio adds a bin and, for 'ordinary' bins, shifts every edge by half
    a bin; comparing it directly can leave the end edges a few 1e-15 inside
    the range, where a position at exactly `lmin` or `lmax` is lost. The ratio
    is therefore compared with the nearest integer at the relative tolerance
    :data:`_BIN_COUNT_RTOL` = 1e-9. That tolerance sits between the two
    scales that matter: the accumulated round-off in the ratio, a few times
    the float64 epsilon (about 1e-15), and the smallest deviation a user can
    intend, which the seven significant digits of the box sizes and bead
    diameters in the simulation filenames put at about 1e-7. It is closer to
    the intent side on purpose, because the two misclassifications are not
    equally bad: calling a true multiple a non-multiple adds a bin and shifts
    the grid, while calling an intended non-multiple a multiple moves the
    edges by less than ``1e-9 * (lmax - lmin)``.

    The symmetric extension of the 'ordinary' bins follows MDAnalysis's
    ``fixedwidth_bins``; the nearest-integer bin count departs from it, for
    the reasons above.

    Time and space are O(`n_bins`).

    References
    ----------
    MDAnalysis ``lib.util.fixedwidth_bins``:
    https://github.com/MDAnalysis/mdanalysis/blob/ae9ff4b738606ba1e7c18a8a9b99a8a87eb69f74/package/MDAnalysis/lib/util.py#L1350

    Examples
    --------
    >>> result = fixedsize_bins(1.0, 0.0, 5.0)
    >>> result['n_bins']
    5
    >>> result['bin_edges']
    array([0., 1., 2., 3., 4., 5.])

    A range that is a whole number of bins up to round-off is tiled exactly:

    >>> result = fixedsize_bins(0.3, 0.0, 2.1)
    >>> edges = result['bin_edges']
    >>> print(result['n_bins'], edges[0], edges[-1])
    7 0.0 2.1

    A radial range whose lower end would fall below zero is clipped there,
    and the bins keep their width:

    >>> result = fixedsize_bins(1.0, 0.1, 2.5, bin_type='nonnegative')
    >>> result['bin_edges']
    array([0., 1., 2., 3.])

    Azimuthal bins of 5 degrees span exactly one period:

    >>> result = fixedsize_bins(np.pi / 36, -np.pi, np.pi, bin_type='periodic')
    >>> print(result['n_bins'], result['bin_edges'][-1] == np.pi)
    72 True
    """
    if lmin >= lmax:
        raise ValueError("'lmin' must be less than 'lmax'.")

    _length = lmax - lmin
    _delta = bin_size
    n_bins, excess = _count_bins(_length, _delta)
    dl = 0.5 * excess  # the excess length to share between the two ends

    if bin_type == 'ordinary':
        # add half of the excess to each end:
        lmin_adj = lmin - dl
        lmax_adj = lmax + dl
    elif bin_type == 'nonnegative':
        if lmin < 0.0:
            raise ValueError(
                "'lmin' must be nonnegative for the 'nonnegative' bin type."
            )
        if lmin - dl < 0.0:
            # clip the lower end at zero and span n_bins full widths from it
            lmin_adj = 0.0
            lmax_adj = n_bins * _delta
        else:
            lmin_adj = lmin - dl
            lmax_adj = lmax + dl
    elif bin_type == 'periodic':
        if excess != 0.0:
            raise ValueError(
                f"For the 'periodic' bin type, the period 'lmax - lmin' "
                f"({_length}) must be a multiple of 'bin_size' ({_delta})."
            )
        lmin_adj, lmax_adj = lmin, lmax
    else:
        invalid_keyword(bin_type, ['ordinary', 'nonnegative', 'periodic'])

    bin_edges = np.linspace(lmin_adj, lmax_adj, n_bins + 1, endpoint=True)
    collectors = np.zeros(n_bins, dtype=np.int64)
    collectors_std = np.zeros(n_bins, dtype=np.float64)

    if save_bin_edges is not None:
        np.save(save_bin_edges, bin_edges)

    return {
        'n_bins': n_bins,
        'bin_edges': bin_edges,
        'collector': collectors,
        'collector_std': collectors_std,
        'range': (lmin_adj, lmax_adj),
    }


def radial_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
) -> np.ndarray:
    """
    Compute the histogram of radial distances from the origin in the spherical
    coordinate system.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    edges : np.ndarray
        A monotonically increasing array of bin edges, including the rightmost
        edge. The last bin includes its right edge. Every position must lie
        within the outermost edges: a coordinate along a periodic axis must
        already be wrapped into the box, since this function does not wrap.

    Returns
    -------
    hist: np.ndarray
        Histogram data

    Raises
    ------
    ValueError
        If the positions array is not two-dimensional.
    ValueError
        If any position falls outside the outermost edges. `np.histogram`
        would drop it silently; here a dropped position means wrong edges or
        an unwrapped coordinate, so it is an error.

    Examples
    --------
    >>> positions= np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    >>> edges = np.array([0, 1, 2, 3])
    >>> radial_histogram(positions, edges)
    array([0, 2, 1])
    """
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim)."
        )

    rad_distances = np.linalg.norm(positions, axis=1)
    return _histogram_or_raise(rad_distances, edges)


def radial_cyl_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    dim: AxisT,
) -> np.ndarray:
    """
    Compute the histogram of radial distances from the longitudinal axis along
    `dim` in the cylindrical coordinate system.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    edges : np.ndarray
        A monotonically increasing array of bin edges, including the rightmost
        edge. The last bin includes its right edge. Every position must lie
        within the outermost edges: a coordinate along a periodic axis must
        already be wrapped into the box, since this function does not wrap.
    dim : :py:data:`AxisT`
        Axis direction (0=x, 1=y, 2=z); see :mod:`polyphys.manage.types`.

    Returns
    -------
    hist: np.ndarray
        Histogram data

    Raises
    ------
    ValueError
        If `dim` is not one of {0, 1, 2}.
    ValueError
        If the positions array is not two-dimensional.
    ValueError
        If any position falls outside the outermost edges. `np.histogram`
        would drop it silently; here a dropped position means wrong edges or
        an unwrapped coordinate, so it is an error.

    Examples
    --------
    >>> import numpy as np
    >>> positions = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3], [0, 3, 4]])
    >>> edges = np.array([0, 1, 2, 3, 4])
    >>> radial_cyl_histogram(positions, edges, dim=1)
    array([1, 1, 0, 2])
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim)."
        )

    trans_axes = np.roll(np.arange(3), -dim)[1:]  # selecting transverse axes
    trans_distances = np.linalg.norm(positions[:, trans_axes], axis=1)
    return _histogram_or_raise(trans_distances, edges)


def axial_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    dim: AxisT,
) -> np.ndarray:
    """
    Compute the histogram of distances from the origin an axis in direction
    `dim`.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    edges : np.ndarray
        A monotonically increasing array of bin edges, including the rightmost
        edge. The last bin includes its right edge. Every position must lie
        within the outermost edges: a coordinate along a periodic axis must
        already be wrapped into the box, since this function does not wrap.
    dim : :py:data:`AxisT`
        The longitudinal axis (0=x, 1=y, 2=z); see
        :mod:`polyphys.manage.types`.

    Returns
    -------
    hist: np.ndarray
        Histogram data

    Raises
    ------
    ValueError
        If `dim` is not one of {0, 1, 2}.
    ValueError
        If the positions array is not two-dimensional.
    ValueError
        If any position falls outside the outermost edges. `np.histogram`
        would drop it silently; here a dropped position means wrong edges or
        an unwrapped coordinate, so it is an error.

    Examples
    --------
    >>> positions= np.array([[1, 2, 3], [2, 3, 4], [1, 2, 2]])
    >>> edges = np.array([0, 1, 2, 3, 4, 5])
    >>> axial_histogram(positions, edges, dim=0)
    array([0, 2, 1, 0, 0])
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim)."
        )

    return _histogram_or_raise(positions[:, dim], edges)


def azimuth_cyl_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    dim: AxisT,
) -> np.ndarray:
    """
    Compute the histogram of azimuth angles in the cylindrical coordinate
    system with the longitudinal axis along `dim`.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    edges : np.ndarray
        A monotonically increasing array of bin edges, including the rightmost
        edge. The last bin includes its right edge. Every position must lie
        within the outermost edges: a coordinate along a periodic axis must
        already be wrapped into the box, since this function does not wrap.
    dim : :py:data:`AxisT`
        The longitudinal axis (0=x, 1=y, 2=z); see
        :mod:`polyphys.manage.types`.

    Returns
    -------
    hist: np.ndarray
        Histogram data

    Raises
    ------
    ValueError
        If `dim` is not one of {0, 1, 2}.
    ValueError
        If the positions array is not two-dimensional.
    ValueError
        If any position falls outside the outermost edges. `np.histogram`
        would drop it silently; here a dropped position means wrong edges or
        an unwrapped coordinate, so it is an error.

    Examples
    --------
    >>> positions= np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0]])
    >>> edges = np.linspace(-np.pi, np.pi, 5)
    >>> azimuth_cyl_histogram(positions, edges, dim=2)
    array([0, 0, 1, 2])
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim)."
        )

    transverse_axes = np.roll(np.arange(3), -dim)[1:]
    azimuthal_angles = np.arctan2(
        positions[:, transverse_axes[1]], positions[:, transverse_axes[0]]
    )
    return _histogram_or_raise(azimuthal_angles, edges)


def planar_cartesian_histogram(
    positions: np.ndarray,
    edges: List[np.ndarray],
    dim: AxisT,
) -> np.ndarray:
    """
    Compute the bi-dimensional histogram in the plan perpendicular to the axis
    along `dim` in the Cartesian coordinate system.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    edges : List[np.ndarray]
        Two monotonically increasing arrays of bin edges, one per planar axis,
        each including its rightmost edge. The last bin along each axis
        includes its right edge. Every position must lie within the outermost
        edges along both axes: a coordinate along a periodic axis must already
        be wrapped into the box, since this function does not wrap.
    dim : :py:data:`AxisT`
        Cartesian axis (0=x, 1=y, 2=z). The right-hand rule is used to pass the
        planar axes to the `np.histogram2d`: When `dim=0` (x), `dim=1` (y) and
        `dim=2` (z) values are passed respectively. When `dim=1` (y), `dim=2`
        (z) and `dim=0` (x) values are passed respectively. When `dim=2` (z),
        `dim=0` (x) and `dim=1` (y) values are passed respectively. See
        :mod:`polyphys.manage.types`.

    Returns
    -------
    hist: np.ndarray
        2D histogram data

    Raises
    ------
    ValueError
        If `dim` is not one of {0, 1, 2}.
    ValueError
        If the positions array is not two-dimensional.
    ValueError
        If any position falls outside the outermost edges. `np.histogram`
        would drop it silently; here a dropped position means wrong edges or
        an unwrapped coordinate, so it is an error.
    ValueError
        If the length of `edges` is not two.

    Examples
    --------
    >>> positions= np.array([[1, 1, 0], [2, 2, 0], [3, 3, 0]])
    >>> edges = [np.array([0, 2, 4]), np.array([0, 2, 4])]
    >>> planar_cartesian_histogram(positions, edges, dim=2)
    array([[1., 0.],
           [0., 2.]])
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim)."
        )
    if len(edges) != 2:
        raise ValueError("'edges' must contain two arrays, one for each axis.")

    t_dims = np.roll(np.arange(3), -dim)[1:]  # selecting transverse axes
    hist, _, _ = np.histogram2d(
        positions[:, t_dims[0]],
        positions[:, t_dims[1]],
        bins=edges,
    )
    n_missed = len(positions) - int(hist.sum())
    if n_missed:
        raise ValueError(
            f"{n_missed} of {len(positions)} positions fall outside the "
            f"edges [{edges[0][0]}, {edges[0][-1]}] x "
            f"[{edges[1][0]}, {edges[1][-1]}]."
        )
    return hist
