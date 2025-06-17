"""
==========================================================
:mod:`polyphys.analyze.measurer`
==========================================================


The :mod:`polyphys.analyze.measurer` module provides a suite of functions for
performing measurements and analyses on molecular simulation data.

Functions
=========

Core Geometric and Structural Measurements:
--------------------------------------------
.. autofunction:: apply_pbc_orthogonal
.. autofunction:: pair_distance
.. autofunction:: end_to_end
.. autofunction:: transverse_size
.. autofunction:: max_distance
.. autofunction:: fsd

Statistical Analysis:
---------------------
.. autofunction:: simple_stats
.. autofunction:: sem

Density and Volume Calculations:
--------------------------------
.. autofunction:: number_density_cube
.. autofunction:: volume_fraction_cube
.. autofunction:: number_density_cylinder
.. autofunction:: volume_fraction_cylinder

Advanced Geometric Calculations:
--------------------------------
.. autofunction:: spherical_segment
.. autofunction:: sphere_sphere_intersection

Binning for histogram processing
--------------------------------
.. autofunction:: create_bin_edge_and_hist
.. autofunction:: fixedsize_bins
.. autofunction:: radial_histogram
.. autofunction:: radial_cyl_histogram
.. autofunction:: axial_histogram
.. autofunction:: azimuth_histogram
.. autofunction:: planar_cartesian_histogram

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
import warnings
from typing import Dict, Tuple, Optional, Literal, Union, Any, List
import numpy as np
from ..manage.utilizer import invalid_keyword


def apply_pbc_orthogonal(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: Dict[int, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Updates the periodic boundary condition (PBC) lengths and their inverses
    based on specified dimensions and lengths in `pbc`.

    Parameters
    ----------
    pbc_lengths : numpy.ndarray
        An array of lengths in each direction.
    pbc_lengths_inverse : numpy.ndarray
        An array of the inverses of lengths in each direction.
    pbc : dict
        A dictionary where keys are dimension indices and values are the
        corresponding lengths for applying PBC.

    Return
    ------
    pbc_lengths : numpy.ndarray
        The updated array of lengths in each direction.
    pbc_lengths_inverse : numpy.ndarray
        The updated array of the inverse lengths in each direction.

    Raises
    ------
    ZeroDivisionError
        If a length value in `pbc` is zero, causing a division by zero in the
        inverse calculation.
    """
    for dim, length in pbc.items():
        if length == 0:
            raise ZeroDivisionError(
                f"Length for dimension {dim} cannot be zero.")
        pbc_lengths[dim] = length
        pbc_lengths_inverse[dim] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def pair_distance(
    positions: np.ndarray,
    pbc: Optional[Dict[int, float]] = None
) -> np.ndarray:
    """
    Computes the pair distance vector between two particles in Cartesian
    coordinates.

    Parameters
    ----------
    positions : numpy.ndarray
        Array of shape (2, n_dim) containing atom positions, and sorted
        by atom number form 1 to N.
    pbc : Optional[Dict[int, float]], default None
        A dictionary with dimension indices as keys (0 for x, 1 for y, and 2
        for z) and lengths of the simulation box in those dimensions as values.
        If provided, applies the minimum image convention.

    Returns
    -------
    dr_ij: np.ndarray
        1D array with pair distances along each axis.

    Raises
    ------
    ValueError
        If `positions` does not contain exactly two sets of coordinates.

    Notes
    -----
    The distance along each axis is the difference between the coordinates of
    the second and the first particle. This function does not compute the
    Euclidean (absolute) distance but instead returns the distance vector
    components.
    """
    n_atoms, n_dims = positions.shape
    if n_atoms != 2:
        raise ValueError("'pair_distance' only works for two atoms.")
    # Calculation in the center of geometry (cog) of the atom group.
    com_no_pbc = np.mean(positions, axis=0)
    shifts = -com_no_pbc
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
        )
        # Move the cog to the unit box
        com_pbc = pbc_lengths * np.around(pbc_lengths_inverse * com_no_pbc)
        # Find the shifts in the cog
        shifts = com_pbc - com_no_pbc
    # Apply the shifts
    positions += shifts
    return positions[1] - positions[0]


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
    geometry of atoms
    """
    # calculation in the center of geometry of the atom group.
    if positions.shape[0] < 2:
        raise ValueError("`positions` must contain at least two atoms.")
    centered_positions = positions - np.mean(positions, axis=0)
    return np.linalg.norm(centered_positions[-1] - centered_positions[0])


def transverse_size(
    positions: np.ndarray,
    axis: Literal[0, 1, 2]
) -> float:
    """
    Computes the mean transverse size of a group of atoms in the plane
    perpendicular to `axis`.

    Parameters
    ----------
    positions : np.ndarray
        Array of shape (n_atoms, n_dim) containing atomic positions, and sorted
        by atom number form 1 to N.
    axis: {0, 1, 2}
        The axis (0 for x, 1 for y, and 2 for z) in the plane perpendicular to
        which the transverse size is calculated.

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
    """
    trans_axes = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
    if axis >= positions.shape[1]:
        raise IndexError("Axis {axis} is out of bounds for positions array"
                         f"with shape {positions.shape}.")
    trans_pos = (positions[:, trans_axes[axis]]
                 - np.mean(positions[:, trans_axes[axis]], axis=0))
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
    """
    # calculation in the center of geometry of the atom group.
    centered_positions = positions - np.mean(positions, axis=0)
    return np.ptp(centered_positions, axis=0)


def fsd(positions: np.ndarray, axis: Literal[0, 1, 2]) -> float:
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

    axis: {0, 1, 2}
        Axis (0 for x, 1 for y, and 2 for z) along which to calculate FSD.

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
    """
    if axis >= positions.shape[1]:
        raise IndexError("Axis {axis} is out of bounds for positions array"
                         f"with shape {positions.shape}.")
    # calculation in the center of geometry of the atom group:
    positions = positions - np.mean(positions, axis=0)
    return np.ptp(positions[:, axis])


def simple_stats(prop: str, array: np.ndarray) -> Dict[str, float]:
    """
    Measures the mean, standard deviation, variance, and standard error
    of the mean (sem) for a property array.

    Parameters
    ----------
    prop : str
        Name of the physical property.
    array: np.ndarray
        Array of property values.

    Returns
    -------
    Dict[str, float]
        Dictionary with keys for mean, variance, and SEM of `prop`.

    Raises
    ------
    ValueError
        If `array` is empty.
    """
    if len(array) == 0:
        raise ValueError("Input array must not be empty.")
    # Unbiased std, var, and sem.
    return {
        prop + '_mean': np.mean(array),
        prop + '_var': np.var(array, ddof=1),
        prop + '_sem': np.std(array, ddof=1) / np.sqrt(len(array))
    }


def sem(data: np.ndarray) -> float:
    """
    Calculates the standard error of the mean (SEM) for a sample.

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
    """
    if len(data) == 0:
        raise ValueError("Input data array must not be empty.")

    return np.std(data, ddof=1) / len(data)**0.5


def number_density_cube(
    n_atom: float,
    d_atom: float,
    l_cube: float,
    pbc: bool = False
) -> float:
    """
    Compute the bulk number density of a species in a cubic box.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cube : float
        Length of one side of the cubic box.
    pbc : bool
        Periodic boundary conditions along all box sides. If `True`,
        :math:`v_{avail} = l_{cube}^3`; otherwise,
        :math:`v_{avail} = (l_{cube} - d_{atom})^3`. Defaults to `False`.

    Returns
    -------
    float
        Bulk number density of the species in the cubic box.

    Notes
    -----
    The bulk number density is calculated as :math:`n_{atom} / v_{avail}`,
    where `v_{avail}` is the available volume to the center of geometry of a
    particle based on the presence of periodic boundary conditions.
    """
    v_avail = (l_cube - int(pbc) * d_atom) ** 3
    return n_atom / v_avail


def volume_fraction_cube(
    n_atom: float,
    d_atom: float,
    l_cube: float,
    pbc: bool = False
) -> float:
    """
    Compute the volume fraction of a species in a cubic box.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cube : float
        Length of one side of the cubic box.
    pbc : bool
        Periodic boundary conditions along all box sides. If `True`,
        :math:`v_{avail} = l_{cube}^3`; otherwise,
        :math:`v_{avail} = (l_{cube} - d_{atom})^3`. Defaults to `False`.

    Returns
    -------
    float
        Volume fraction of the species in the cubic box.

    Notes
    -----
    The volume fraction is computed in the volume available to the center of
    geometry of each particle. For point-like particles, the notion of volume
    fraction is meaningless. For finite-size particles, the available volume
    depends on whether periodic boundary conditions (PBCs) are applied.
    """
    rho = number_density_cube(n_atom, d_atom, l_cube, pbc)
    return rho * np.pi * d_atom ** 3 / 6


def number_density_cylinder(
    n_atom: float,
    d_atom: float,
    l_cyl: float,
    d_cyl: float,
    pbc: bool = False
) -> float:
    """
    Compute the bulk number density of a species in a cylindrical confinement.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cyl : float
        Length of the cylindrical confinement.
    d_cyl : float
        Diameter of the cylindrical confinement.
    pbc : bool
        Periodic boundary conditions along the longitudinal axis. If `True`,
        :math:`v_{avail} = \\pi * l_{cyl} * (d_{cyl} - d_{atom})^2 / 4`;
        otherwise, :math:`v_{avail} = \\pi * (l_{cyl} - d_{atom}) * (d_{cyl}
        - d_{atom})^2 / 4`. Defaults to `False`.

    Returns
    -------
    float
        Bulk number density of the species in the cylindrical confinement.

    Notes
    -----
    The bulk number density is calculated as :math:`n_{atom} / v_{avail}`,
    where `v_{avail}` is the available volume to the center of geometry of a
    particle based on the presence of periodic boundary conditions.
    """
    v_avail = np.pi * (l_cyl - int(pbc) * d_atom) * (d_cyl - d_atom) ** 2 / 4
    return n_atom / v_avail


def volume_fraction_cylinder(
    n_atom: float,
    d_atom: float,
    l_cyl: float,
    d_cyl: float,
    pbc: bool = False
) -> float:
    """
    Compute the volume fraction of a species in a cylindrical confinement.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cyl : float
        Length of the cylindrical confinement.
    d_cyl : float
        Diameter of the cylindrical confinement.
    pbc : bool
        Periodic boundary conditions along the longitudinal axis. If `True`,
        :math:`v_{avail} = \\pi * l_{cyl} * (d_{cyl} - d_{atom})^2 / 4`;
        otherwise, :math:`v_{avail} = \\pi * (l_cyl - d_{atom}) * (d_{cyl}
        - d_{atom})^2 / 4`. Defaults to `False`.

    Returns
    -------
    float
        Volume fraction of the species in the cylindrical confinement.

    Notes
    -----
    The volume fraction is computed in the volume available to the center of
    geometry of each particle. For point-like particles, the notion of volume
    fraction is meaningless. For finite-size particles, the available volume
    depends on whether periodic boundary conditions (PBCs) are applied.
    """
    rho = number_density_cylinder(n_atom, d_atom, l_cyl, d_cyl, pbc)
    return rho * np.pi * d_atom ** 3 / 6


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
    .. [1] Weisstein, Eric W. "Spherical Segment." From MathWorld--A Wolfram
    Web Resource.
       https://mathworld.wolfram.com/SphericalSegment.html

    Examples
    --------
    >>> spherical_segment(3, 1, 2)
    37.69911184307752
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
    .. [1] Weisstein, Eric W. "Sphere-Sphere Intersection."
       From MathWorld--A Wolfram Web Resource.
       https://mathworld.wolfram.com/Sphere-SphereIntersection.html

    Examples
    --------
    >>> sphere_sphere_intersection(3, 4, 2)
    75.39822368615503
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
    vol = (np.pi / (12 * d)) * (r_max + r_min - d)**2 * (
        d**2 + 2 * d * (r_max + r_min)
        - 3 * (r_min**2 + r_max**2)
        + 6 * r_min * r_max
    )

    return vol


def create_bin_edge_and_hist(
    bin_size: float,
    lmin: float,
    lmax: float,
    output: Optional[str] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create bin edges and an empty histogram for data processing.

    Parameters
    ----------
    bin_size : float
        Size of each bin.
    lmin : float
        Lower bound of the system in the direction of interest.
    lmax : float
        Upper bound of the system in the direction of interest.
    output : str
        Filename (including filepath) to which bin edges array is saved. A
        `.npy` extension will be appended to the filename if it does not
        already have one.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        bin_edges : np.ndarray
            The edges to pass into a histogram.
        hist : np.ndarray
            An empty histogram array initialized with zeros.

    Raises
    ------
    ValueError
        If `lmin` is not less than `lmax`.

    Notes
    -----
    The `.npy` extension is handled internally by `np.save`.
    """
    if lmin >= lmax:
        raise ValueError("'lmin' must be less than 'lmax'.")
    bin_edges = np.arange(lmin, lmax + bin_size, bin_size)
    hist = np.zeros(len(bin_edges) - 1, dtype=np.int16)

    if output is not None:
        # Save bin edges to file if save_to path is provided
        # f"{save_to}{sim_name}-{edge_name}.npy"
        np.save(output, bin_edges)

    return bin_edges, hist


def fixedsize_bins(
    bin_size: float,
    lmin: float,
    lmax: float,
    bin_type: Literal['ordinary', 'nonnegative', 'periodic'] = 'ordinary',
    save_bin_edges: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generate bin edges and an empty histogram, adjusting `lmin` and `lmax` to
    ensure consistent `bin_size`.

    Parameters
    ----------
    bin_size : float
        Size of each bin.
    lmin : float
        Lower bound of the system in the direction of interest.
    lmax : float
        Upper bound of the system in the direction of interest.
    bin_type : {'ordinary', 'nonnegative', 'periodic'}, default 'ordinary'
        Type of bin:
            - 'ordinary': Extends `lmin` and `lmax` symmetrically. Examples are
              Cartesian coordinates and spherical polar coordinate.
            - 'nonnegative': Extends `lmin` are `lmax` are symmetrically  if
              adjusted `lmin` is nonnegative; otherwise, `lmin` is set to 0 and
              `lmax` is extended. Examples are radial directions in the polar,
              spherical, and cylindrical coordinate systems.
            - 'periodic': Periodic binning (e.g., azimuthal coordinate).
              Examples are azimuthal directions in cylindrical and spherical
            coordinate systems.
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
            Array initialized for histogram values.
        - 'collector_std' : np.ndarray
            Array initialized for standard deviation values.
        - 'range' : Tuple[float, float]
            Updated range of bins (`lmin`, `lmax`).

    Raises
    ------
    ValueError
        If `lmin` is not less than `lmax`.

    Notes
    -----
    The `.npy` extension is handled internally by `np.save`.

    References
    ----------
    https://docs.mdanalysis.org/1.1.1/documentation_pages/lib/util.html#MDAnalysis.analysis.density.fixedwidth_bins
    """
    if lmin >= lmax:
        raise ValueError("'lmin' must be less than 'lmax'.")

    _length = lmax - lmin
    _delta = bin_size

    if bin_type == 'ordinary':
        n_bins = int(np.ceil(_length / _delta))
        dl = 0.5 * (n_bins * _delta - _length)  # excess length
        # add half of the excess to each end:
        lmin_adj = lmin - dl
        lmax_adj = lmax + dl
        bin_edges = np.linspace(lmin_adj, lmax_adj, n_bins + 1, endpoint=True)
    elif bin_type == 'nonnegative':
        n_bins = int(np.ceil(_length / _delta))
        dl = 0.5 * (n_bins * _delta - _length)
        lmin_adj = max(0.0, lmin - dl)
        lmax_adj = lmax + 2 * dl if lmin_adj == 0 else lmax + dl
        bin_edges = np.linspace(lmin_adj, lmax_adj, n_bins + 1, endpoint=True)
    elif bin_type == 'periodic':
        n_bins = int(np.ceil(_length / _delta))
        lmin_adj, lmax_adj = lmin, lmax
        bin_edges = np.arange(lmin_adj, lmax_adj + _delta, _delta)
        if (len(bin_edges) - 1) != n_bins:
            # Number of bins (n_bins='{n_bins}') is different from the actual
            # number of bins (n_edges-1={len(bin_edges)-1}) for the 'periodic'
            # bin type because period, i.e., 'lmax-lmin={_length}', and
            # delta={_delta}, not 'n_bins', are used to created 'bin_edges'.
            warnings.warn("'n_bins' is set to 'len(bin_edges)-1'", UserWarning)
    else:
        invalid_keyword(bin_type, ['ordinary', 'nonnegative', 'periodic'])

    hist_collectors = np.zeros(n_bins, dtype=np.int16)
    hist_collectors_std = np.zeros(n_bins, dtype=np.int16)

    if save_bin_edges is not None:
        np.save(save_bin_edges, bin_edges)

    return {
        'n_bins': n_bins,
        'bin_edges': bin_edges,
        'collector': hist_collectors,
        'collector_std': hist_collectors_std,
        'range': (lmin_adj, lmax_adj)
    }


def radial_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    bin_range: Tuple[float, float],
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
        edge.
    bin_range : Tuple[float, float]
        The lower and upper ranges of the bins.

    Returns
    -------
    hist: np.ndarray
        Histogram data

    Raises
    ------
    ValueError
        If the positions array is not two-dimensional.
    """
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim).")

    rad_distances = np.linalg.norm(positions, axis=1)
    hist, _ = np.histogram(rad_distances, bins=edges, range=bin_range)
    return hist


def radial_cyl_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    bin_range: Tuple[float, float],
    dim: Literal[0, 1, 2]
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
        edge.
    bin_range : Tuple[float, float]
        The lower and upper ranges of the bins.
    dim : {0, 1, 2}
        The longitudinal axis (0=x, 1=y, 2=z).

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
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim).")

    trans_axes = np.roll(np.arange(3), -dim)[1:]  # selecting transverse axes
    trans_distances = np.linalg.norm(positions[:, trans_axes], axis=1)
    hist, _ = np.histogram(trans_distances, bins=edges, range=bin_range)
    return hist


def axial_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    bin_range: Tuple[float, float],
    dim: Literal[0, 1, 2]
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
        edge.
    bin_range : Tuple[float, float]
        The lower and upper ranges of the bins.
        The lower and upper ranges of the bins.
    dim : {0, 1, 2}
        Axis direction (0=x, 1=y, 2=z).

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
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim).")

    hist, _ = np.histogram(positions[:, dim], bins=edges, range=bin_range)
    return hist


def azimuth_cyl_histogram(
    positions: np.ndarray,
    edges: np.ndarray,
    bin_range: Tuple[float, float],
    dim: Literal[0, 1, 2]
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
        edge.
    bin_range : Tuple[float, float]
        The lower and upper ranges of the bins.
    dim : {0, 1, 2}
        The longitudinal axis (0=x, 1=y, 2=z).

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
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim).")

    transverse_axes = np.roll(np.arange(3), -dim)[1:]
    azimuthal_angles = np.arctan2(
        positions[:, transverse_axes[1]], positions[:, transverse_axes[0]]
        )
    hist, _ = np.histogram(azimuthal_angles, bins=edges, range=bin_range)
    return hist


def planar_cartesian_histogram(
    positions: np.ndarray,
    edges: List[np.ndarray],
    bin_ranges: List[Tuple[int, int]],
    dim: Literal[0, 1, 2]
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
        A monotonically increasing array of bin edges, including the rightmost
        edge.
    bin_range : List[Tuple[float, float]]
        The list of the lower and upper ranges of the bins in the transverse
        directions within the plane.
    dim : {0, 1, 2}
        Cartesian axis (0=x, 1=y, 2=z). The right-hand rule is used to pass the
        planar axes to the `np.histogram2d`: When `dim=0` (x), `dim=1` (y) and
        `dim=2` (z) values are passed respectively. When `dim=1` (y), `dim=2`
        (z) and `dim=0` (x) values are passed respectively. When `dim=2` (z),
        `dim=0` (x) and `dim=1` (y) values are passed respectively.

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
        If the length of `edges` or `bin_ranges` is not two.
    """
    if dim not in {0, 1, 2}:
        raise ValueError("'dim' must be one of {0, 1, 2}.")
    if positions.ndim != 2:
        raise ValueError(
            "'positions' must be a 2D array with shape (n_atoms, n_dim).")
    if len(edges) != 2:
        raise ValueError("'edges' must contain two arrays, one for each axis.")
    if len(bin_ranges) != 2:
        raise ValueError(
            "'bin_edges' must contain two tuples, one for each axis."
        )

    t_dims = np.roll(np.arange(3), -dim)[1:]  # selecting transverse axes
    hist, _, _ = np.histogram2d(
        positions[:, t_dims[0]],
        positions[:, t_dims[1]],
        bins=edges,
        range=bin_ranges)
    return hist
