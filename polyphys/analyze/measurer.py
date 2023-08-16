"""
Contains tools for performing statistical measurements such as averaging,
standard deviation, or the like on a given time series, correlation data,
or histogram.
"""
from typing import Any, Dict
import numpy as np
import MDAnalysis as mda


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


def transverse_size(atomgroup: mda.AtomGroup) -> np.floating:
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
    transverse_pos = atomgroup.positions[:, :2]  # x,y not z
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


def fsd(
    positions: np.ndarray,
    axis: int = 2
) -> np.ndarray:
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


def size_ratio(dcrowd, dmon=1):
    """Set the type of size ratio of crowders to monomers."""
    if dcrowd < dmon:
        ratio = r"$a_c < 1$"
    elif dcrowd == dmon:
        ratio = r"$a_c = 1$"
    else:
        ratio = r"$a_c > 1$"
    return ratio


def size_ratio_equal(dcrowd, dmon=1):
    """Set the type of size ratio of crowders to monomers."""
    if dcrowd < dmon:
        ratio = r"$a_c \leq 1$"
    else:
        ratio = r"$a_c > 1$"
    return ratio
