"""
Contains tools for performing various measurements, such as averaging, on data.
"""
from typing import Any, Dict, Literal, Tuple, Optional
import numpy as np


def apply_pbc_orthogonal(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: Dict[int, float]
) -> Tuple[np.ndarray, np.ndarray]:
    """Updates the periodic boundary condition (PBC) information in the
    dimensions (keys) by the length (values) passed by `pbc`. The key
    are indexes of the `pbc_lengths` and `pbc_lengths_inverse` arrays.

    Parameters
    ----------
    pbc_lengths: numpy.ndarray
        An array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        An array of the inverses of lengths/sizes in different directions.
    pbc: dict
        A dictionary of dimensions (keys) and lengths (values) in which
        the pbc exists.

    Return
    ------
    pbc_lengths: numpy.ndarray
        The updated array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        The updated array of the inverses of lengths/sizes in different
        directions.
    """
    for dim, length in pbc.items():
        pbc_lengths[dim] = length
        pbc_lengths_inverse[dim] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def pair_distance(
    positions: np.ndarray,
    pbc: Optional[Dict[int, float]] = None
) -> np.ndarray:
    """
    Computes the pair distance between two particles along each axis in the
    Cartesian coordinate system, applying the minimum image convention if
    periodic boundary conditions (PBC) are provided.

    Note
    ----
    The distance along each axis is the difference between the coordinates of
    the second and the first particle. This function does not compute the
    Euclidean (absolute) distance but instead returns the distance vector
    components.

    Parameters
    ----------
    positions : numpy.ndarray
        A (2, n_dim) array containing the coordinates of the two atoms.
        The array should be sorted by atom number, with the first row
        corresponding to the first atom and the second row to the second atom.

    pbc : Optional[Dict[int, float]], default None
        A dictionary where the keys are the dimensions (0 for x, 1 for y,
        and 2 for z) and the values are the lengths of the simulation box
        in those dimensions. If provided, the minimum image convention
        will be applied.

    Returns
    -------
    dr_ij: np.ndarray
        A 1D numpy array containing the pair distances along each axis.
    """
    n_atoms, n_dims = positions.shape
    if n_atoms != 2:
        raise ValueError("'pair_distance' only works for two atoms.")
    # Calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
        )
        positions -= pbc_lengths * np.around(pbc_lengths_inverse * positions)
    # Using map to find the r_ij values
    dr_ij = np.array(
        list(map(lambda i: positions[1, i] - positions[0, i], range(n_dims))))
    return dr_ij


def simple_stats(property_name, array) -> Dict[str, Any]:
    """
    Measures the mean, standard deviation, variance, and standard error
    of the mean (sem) for an array.

    Parameters
    ----------
    property_name : str
        Name of physical property.
    array: numpy array of float
        Values of `property_name`.

    Return
    ------
    stats: Dict[str, Any]
        an dict of mean, std, var, and sem of `property_name`
    """
    # Unbiased std, var, and sem.
    stat = {
        property_name + '_mean': np.mean(array),
        property_name + '_var': np.var(array, ddof=1),
        property_name + '_sem': np.std(array, ddof=1) / np.sqrt(len(array))
    }
    return stat


def end_to_end(positions) -> Any:
    """
    Measures the end-to-end distance of a linear polymer, in the frame of
    reference located at the polymer's center of geometry.

    `positions` is sorted by atom id, so the end-to-end distance is simply
    the difference between last and first items in `positions`

    Parameters
    ----------
    positions : numpy array of float
        Positions (an n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.

    Return
    ------
    end_to_end: numpy array of float
    """
    # calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    positions = positions[-1] - positions[0]
    return np.linalg.norm(positions)


def transverse_size(positions: np.ndarray) -> np.floating:
    """
    Measures the mean transverse size of a linear polymer, in the frame of
    reference located at the polymer's center of geometry.

    Parameters
    ----------
    atomgroup : atomgroup
        An MDAnalysis AtomGroup for which the transverse size is measured.

    Return
    ------
    float: Twice of the maximum of the transverse (radial) distances of a group
    of atoms.
    """
    # calculation in the center of geometry of the atom group.
    # coordinates change for each frame
    transverse_pos = positions[:, :2]  # x,y not z
    transverse_pos = transverse_pos - np.mean(transverse_pos, axis=0)
    # times by 2, so we have the diameter not radius:
    return 2 * np.mean(np.linalg.norm(transverse_pos, axis=1))


def max_distance(positions: np.ndarray) -> np.ndarray:
    """
    Measures the maximum distance in each of the three Cartesian direction,
    in the frame of reference located at the polymer's center of geometry.

    The maximum distance is computed by subtracting the max and min of all
    the particle in an atom group in a given frame or snapshot or time step.

    Parameters
    ----------
    positions: numpy array of dtype float
        Positions (an n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.

    Return
    ------
    [xmax, ymax, zmax]: numpy array of  float
    """
    # calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    xmax = np.abs(np.amax(positions[:, 0]) - np.amin(positions[:, 0]))
    ymax = np.abs(np.amax(positions[:, 1]) - np.amin(positions[:, 1]))
    zmax = np.abs(np.amax(positions[:, 2]) - np.amin(positions[:, 2]))
    return np.array([xmax, ymax, zmax])


def fsd(positions: np.ndarray, axis: int = 2) -> np.ndarray:
    """
    Calculates the average size/diameter of a polymer confined in a
    cylindrical geometry based on the farthermost distance concept.

    fsd stands for Feret's statistical diameter: other names are the mean
    span dimension, the farthermost distance, or the mean caliper diameter.

    Parameters
    ----------
    positions: numpy array of  float
        Positions (an n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.
    axis: int or tuple of int, default 2
        The index of the axis of the cylinder; by default it is the z axis.
        It can be any integer in the range of the spatial dimension of the
        system.

    Return
    ------
    fsd: numpy array of  float

    References
    ----------
    "A Theoretical Study of the Separation Principle in Size Exclusion
    Chromatography", Wang Y Teraoka
    I Hansen FY Peters GH Ole H. Macromolecules 2010, 43, 3, 1651-1659
    https://doi.org/10.1021/ma902377g
    """
    # calculation in the center of geometry of the atom group:
    positions = positions - np.mean(positions, axis=0)
    positions = np.ptp(positions[:, axis])
    return positions


def sem(data: np.ndarray) -> float:
    """Calculates the standard error of mean for a given array-like sample
    data.

    ddof is set to 1 since the standard error of mean is measured for the
    sample data not the population.

    Parameters
    ----------
    data: array-like
        A array-like (or numpy 1D ndarray) data

    Return
    ------
    float:
        The standard error of mean of `data`

    Requirements
    ------------
    Numpy
    """
    return np.std(data, ddof=1) / len(data) ** 0.5


def size_ratio(
    dcrowd,
    dmon=1
) -> Literal['$a_c < 1$'] | Literal['$a_c = 1$'] | Literal['$a_c > 1$']:
    """Set the type of size ratio of crowders to monomers."""
    if dcrowd < dmon:
        return r"$a_c < 1$"
    if dcrowd == dmon:
        return r"$a_c = 1$"
    if dcrowd > dmon: 
        return r"$a_c > 1$"    


def size_ratio_equal(
    dcrowd,
    dmon=1
) -> Literal['$a_c \\leq 1$'] | Literal['$a_c > 1$']:
    """Set the type of size ratio of crowders to monomers."""
    if dcrowd < dmon:
        ratio = r"$a_c \leq 1$"
    else:
        ratio = r"$a_c > 1$"
    return ratio


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
    The bulk number density is calculated as :math:`n_{atom} / v_{avail`,
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
    Volume fraction is computed in the volume available to the center of
    geometry of each particle. For point-like particles, the notion of the
    volume fraction is meaningless. For finit-size particles, the available
    volume depends on whether periodic boundary conditions (PBCs) are applied.
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
        Number of particles in the cylindrical confinement.
    d_atom : float
        Diameter of the particles of the species.
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
    The bulk number density is calculated as :math:`n_{atom} / v_{avail`,
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
        Number of particles in the cylindrical confinement.
    d_atom : float
        Diameter of the particleof the species.
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
    Volume fraction is computed in the volume available to the center of
    geometry of each particle. For point-like particles, the notion of the
    volume fraction is meaningless. For finit-size particles, the available
    volume depends on whether periodic boundary conditions (PBCs) are applied
    """
    rho = number_density_cylinder(n_atom, d_atom, l_cyl, d_cyl, pbc)
    return rho * np.pi * d_atom ** 3 / 6
