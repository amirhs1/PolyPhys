"""\
==========================================================
:mod:`polyphys.analyze.measurer`
==========================================================
The :mod:`polyphys.analyze.measurer` module provides a suite of functions for
performing a varietyof measurements and analyses on molecular or particle-based
simulation data. This module is designed to assist with structural and
geometric calculations, statistical analysis, and handling of periodic boundary
conditions (PBCs) in orthogonal or confined simulation environments.

The module is particularly useful in molecular dynamics (MD) simulations, where
quantifying characteristics like end-to-end distance, transverse size, number
density, and volume fraction can help describe the structure and behavior of
polymers, macromolecules, or particle assemblies. Measurements span from single
particle interactions to bulk properties in both cubic and cylindrical
confinement geometries.

Functions
=========
.. autofunction:: apply_pbc_orthogonal
.. autofunction:: pair_distance
.. autofunction:: end_to_end
.. autofunction:: transverse_size
.. autofunction:: max_distance
.. autofunction:: fsd
.. autofunction:: simple_stats
.. autofunction:: sem
.. autofunction:: number_density_cube
.. autofunction:: volume_fraction_cube
.. autofunction:: number_density_cylinder
.. autofunction:: volume_fraction_cylinder

Dependencies
============
- `numpy`: For numerical operations on arrays, such as distances, means, and
  statistical calculations.
- `typing`: For type hinting, particularly with `Optional` and `Literal` types
  for enhanced function clarity.

Usage
=====
The functions in this module allow users to perform essential measurements on
spatial data, compute densities, apply boundary conditions, and gather
statistical properties. Many functions are designed with flexibility in mind,
allowing them to be applied across a range of MD simulation setups.

Examples
========
Example of computing the end-to-end distance for a polymer chain:

>>> import numpy as np
>>> import measure
>>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
>>> distance = measure.end_to_end(positions)
>>> print(distance)
3.4641016151377544

To calculate the volume fraction of particles in a cubic box with periodic
boundary conditions:

>>> n_atom = 100
>>> d_atom = 1.0
>>> l_cube = 10.0
>>> fraction = measure.volume_fraction_cube(n_atom, d_atom, l_cube, pbc=True)
>>> print(fraction)

Notes
=====
Functions in this module frequently rely on assumptions about input data:
- Atomic or particle positions are expected to be ordered and centered,
  typically around a center of geometry (COG) for accurate geometric
  calculations.
- When using functions that apply periodic boundary conditions, ensure that the
  PBC lengths are correctly specified to avoid calculation errors.

References
==========
For Feret's statistical diameter:
- Wang Y, Teraoka I, Hansen FY, Peters GH, Ole H. "A Theoretical Study of
  the Separation Principle in Size Exclusion Chromatography." Macromolecules
  2010, 43, 3, 1651-1659. https://doi.org/10.1021/ma902377g
"""
from typing import Dict, Tuple, Optional, Literal
import numpy as np


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


def end_to_end(positions: np.ndarray) -> float:
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
        perpendicualr to `axis`.

    Raises
    ------
    IndexError
        If `axis` is out of bounds for `positions`.

    Notes
    -----
    The distance is calculated in a reference frame centered at the center of
    geometry of atoms in the plane perpendicaulr to a given axis.
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
    Dict[str, flaot]
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
    pbc: Optional[bool] = False
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
    pbc : bool, optional
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
    pbc: Optional[bool] = False
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
    pbc : bool, optional
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
    pbc: Optional[bool] = False
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
    pbc : bool, optional
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
    pbc: Optional[bool] = False
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
    pbc : bool, optional
        Periodic boundary conditions along the longitudinal axis. If `True`,
        :math:`v_{avail} = \\pi * l_{cyl} * (d_{cyl} - d_{atom})^2 / 4`;
        otherwise, :math:`v_{avail} = \\pi * (l_cyl - d_{atom}) * (d_{cyl}
        - d_{atom})^2 / 4`. Defaults to `True`.

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
