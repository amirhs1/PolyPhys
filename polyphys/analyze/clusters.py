"""A submodule to detect and analyze particle clusters in a given system. It is
worth noting there is a difference analyzing particle clusters and
trajectory (or ensemble) clusters.
Below, different algorithms found in literature are implemented in scientific
Python.
"""
from typing import Dict, List, Tuple, Optional, Union
from itertools import combinations
from collections import defaultdict
import numpy as np
import pandas as pd
import numpy.linalg as npla

from ..manage.parser import (TransFociCub, TransFociCyl,
                             SumRuleCubHeteroLinear,
                             SumRuleCubHeteroRing)
from ..manage.typer import AxisT, TopologyT
from .measurer import apply_pbc_orthogonal


def self_dist_array(
    positions: np.ndarray,
    pbc: Optional[Dict[AxisT, float]] = None
) -> np.ndarray:
    """
    Compute the pairwise Euclidean distances between atoms.

    Parameters
    ----------
    positions : np.ndarray, shape (n_atoms, n_dims)
        Cartesian coordinates of atoms. Each row corresponds to an atom, and
        each column corresponds to a spatial dimension.
    pbc : Dict[AxisT, float]
        Periodic boundary conditions (PBC) as a dictionary where keys are
        dimensions (0 for 'x', 1 for 'y', 2 for 'z'), and values are the box
        lengths along those dimensions. If None, PBC is not applied. Default is
        None.

    Returns
    -------
    np.ndarray, shape (n_atoms, n_atoms)
        A symmetric matrix where each element (i, j) is the center-to-center
        distance between atoms i and j.

    Examples
    --------
    >>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    >>> self_dist_array(positions)
    array([[0.        , 1.73205081, 3.46410162],
           [1.73205081, 0.        , 1.73205081],
           [3.46410162, 1.73205081, 0.        ]])

    >>> pbc = {0: 10.0, 1: 10.0, 2: 10.0}
    >>> self_dist_array(positions, pbc)
    array([[0.        , 1.73205081, 3.46410162],
           [1.73205081, 0.        , 1.73205081],
           [3.46410162, 1.73205081, 0.        ]])

    """
    n_atoms, n_dims = positions.shape
    # Calculate pairwise distance matrix
    pos_t = positions.T
    pos_i = np.reshape(pos_t, (n_dims, 1, n_atoms))
    pos_j = np.reshape(pos_t, (n_dims, n_atoms, 1))
    dist_ji = pos_j - pos_i

    # Calculate PBC lengths and inverses, and apply PBCs if PBCs are present
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
            )
        pbc_lengths = pbc_lengths.reshape(n_dims, 1, 1)
        pbc_lengths_inverse = pbc_lengths_inverse.reshape(n_dims, 1, 1)
        dist_ji -= pbc_lengths * np.around(pbc_lengths_inverse * dist_ji)

    # Calculate squared distances and return square root of matrix
    dist_ji = np.sqrt(np.sum(np.square(dist_ji), axis=0))
    np.fill_diagonal(dist_ji, 0)
    return dist_ji


def self_dist_array_opt(
    positions: np.ndarray,
    pbc: Optional[Dict[AxisT, float]] = None
) -> np.ndarray:
    """
    Compute the pairwise Euclidean distances between atoms.

    Parameters
    ----------
    positions : np.ndarray
        Cartesian coordinates of atoms. Each row corresponds to an atom, and
        each column corresponds to a spatial dimension.
    pbc : Dict[AxisT, float]
        Periodic boundary conditions (PBC) as a dictionary where keys are
        dimensions (0 for 'x', 1 for 'y', 2 for 'z'), and values are the box
        lengths along those dimensions. If None, PBC is not applied. Default is
        None.

    Return
    ------
    np.ndarray
        A square matrix of size n_atoms * n_atoms in which each element are
        the center-to-center distances between different pairs of atoms. All
        the diagonal elements of this matrix is 0.

    Notes
    -----
    In 3 dimensions, this method is not much faster the `self_dist_array`.

    Examples
    --------
    >>> positions = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    >>> self_dist_array_opt(positions)
    array([[0.        , 1.73205081, 3.46410162],
           [1.73205081, 0.        , 1.73205081],
           [3.46410162, 1.73205081, 0.        ]])

    >>> pbc = {0: 10.0, 1: 10.0, 2: 10.0}
    >>> self_dist_array_opt(positions, pbc)
    array([[0.        , 1.73205081, 3.46410162],
           [1.73205081, 0.        , 1.73205081],
           [3.46410162, 1.73205081, 0.        ]])

    References
    ----------
    Albanie S. (2019) Euclidean Distance Matrix Trick. In: Advances in Computer
    Vision. EECV 2018. Lecture Notes in Computer Science, vol 11131. Springer,
    Cham. https://doi.org/10.1007/978-3-030-01446-4_2
    """
    n_atoms, n_dims = positions.shape
    # Redefining the shapes of different arrays to combine linear algebra
    # with numpy broadcasting.
    gram = np.matmul(positions, positions.T)
    dist_sq = np.diag(gram).reshape(n_atoms, 1) + \
        np.diag(gram).reshape(1, n_atoms) - 2 * gram

    # Calculate PBC lengths and inverses, and apply PBCs if PBCs are present
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc_orthogonal(
            np.zeros(n_dims), np.zeros(n_dims), pbc
        )
        pbc_lengths = pbc_lengths**2
        pbc_lengths_inverse = pbc_lengths_inverse**2
        dist_sq = dist_sq - pbc_lengths * np.around(
            pbc_lengths_inverse * dist_sq
            )
    return np.sqrt(dist_sq)


def find_direct_contacts(
    dist: np.ndarray,
    cutoff: float,
    inclusive: bool = True
) -> np.ndarray:
    """
    Identify direct contacts between atoms based on a cutoff distance.

    Two atoms are in direct contact if the center-to-center distance between
    two atoms being within the specified cutoff distance. The inclusivity of
    the cutoff can be controlled using the `inclusive` parameter.

    Parameters
    ----------
    dist : np.ndarray, shape (n_atoms, m_atoms)
        Pairwise distance matrix where element (i, j) is the distance between
        atoms `i` and `j`.
    cutoff : float
        The threshold distance for defining direct contact. Must be a positive
        number.
    inclusive : bool, optional
        If True, includes pairs of atoms where the distance equals the cutoff.
        If False, only considers distances strictly less than the cutoff.
        Default is True.

    Returns
    -------
    np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where element (i, j) is 1 if the distance between atoms
        `i` and `j` is within the cutoff, and 0 otherwise. If the input `dist`
        matrix is symmetric, the output will also be symmetric, with all
        diagonal elements set to 1.

    Raises
    ------
    ValueError
        If the `cutoff` is not a positive number.

    Examples
    --------
    >>> distances = np.array([[0.0, 1.2, 3.4], [1.2, 0.0, 2.3],
                             [3.4, 2.3, 0.0]])
    >>> find_direct_contacts(distances, cutoff=2.0)
    array([[1, 1, 0],
           [1, 1, 1],
           [0, 1, 1]])

    >>> find_direct_contacts(distances, cutoff=2.0, inclusive=False)
    array([[1, 1, 0],
           [1, 1, 0],
           [0, 0, 1]])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://aip.scitation.org/doi/10.1063/1.454720 .
    """
    if cutoff <= 0:
        raise ValueError("'cutoff' must be a positive number.")

    direct_contacts = dist <= cutoff if inclusive else dist < cutoff
    return np.asarray(direct_contacts, dtype=int)


def generate_contact_matrix(direct_contacts: np.ndarray) -> np.ndarray:
    """
    Generate a symmetric contact matrix that includes all direct and
    indirect contacts.

    This function generates a symmetric binary matrix representing all direct
    and indirect contacts between atoms. Starting from the `direct_contacts`
    matrix, it iteratively checks pairs of columns to identify and propagate
    connections.

    Parameters
    ----------
    direct_contacts : np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where element (i, j) is 1 if atoms `i` and `j` have a
        direct contact, and 0 otherwise.

    Returns
    -------
    np.ndarray, shape (n_atoms, m_atoms)
        A binary matrix where element (i, j) is 1 if atoms `i` and `j` are
        directly or indirectly connected, and 0 otherwise.

    Examples
    --------
    >>> direct_contacts = np.array([[1, 1, 0],
    ...                             [1, 1, 1],
    ...                             [0, 1, 1]])
    >>> generate_contact_matrix(direct_contacts)
    array([[1, 1, 1],
           [1, 1, 1],
           [1, 1, 1]])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720

    """
    _, n_cols = direct_contacts.shape
    contact_matrix = direct_contacts.copy()

    # Iteratively compare pairs of columns and propagate connections
    for i, j in combinations(range(n_cols), 2):
        # If columns i and j share any common contacts, merge them
        if np.any(contact_matrix[:, i] & contact_matrix[:, j]):
            merged_column = contact_matrix[:, i] | contact_matrix[:, j]
            contact_matrix[:, i] = merged_column
            contact_matrix[:, j] = merged_column

    return contact_matrix


def random_full_contact_matrix(
    n_atoms: int,
    seed: int = 0,
    p: float = 0.2,
    save_to: Optional[str] = None
) -> np.ndarray:
    """
    Generate a random direct-contact matrix for testing clustering algorithms.

    This function creates a symmetric binary matrix where each element
    indicates whether a direct contact exists between a pair of atoms. It is
    useful for testing clustering algorithms such as `count_bonds` and
    `count_clusters`.

    Parameters
    ----------
    n_atoms : int
        Number of atoms or nodes in the system.
    seed : int, optional
        Seed for the random number generator. Default is 0.
    p : float, optional
        Probability of a contact existing between any two atoms. Must be
        between 0 and 1. Default is 0.2.
    save_to : str, optional
        The filename or filepath to which the contact martix is saved. If None,
        the matrix is not saved. Default is None.

    Returns
    -------
    np.ndarray, shape (n_atoms, n_atoms)
        A symmetric binary matrix where element `(i, j)` is 1 if a contact
        exists between atoms `i` and `j`, and 0 otherwise. All diagonal
        elements are 1.

    See Also
    --------
    `find_direct_contact`: Identify direct contacts between atoms based on a
    cutoff distance.

    Examples
    --------
    >>> random_direct_contacts(4, seed=42)
    array([[1, 1, 0, 1],
           [1, 1, 1, 0],
           [0, 1, 1, 0],
           [1, 0, 0, 1]])
    """
    # Creating direct-contacts matrix
    np.random.seed(seed)
    dir_contacts = np.random.binomial(n=1, p=p, size=(n_atoms, n_atoms))
    dir_contacts_up = np.triu(dir_contacts, k=1)  # Exclude diag
    # Create a symmetric matrix and add diagonal elements
    dir_contacts = dir_contacts_up + dir_contacts_up.T
    np.fill_diagonal(dir_contacts, 1)
    # Creating the full direct contacts matrix
    full_contacts = generate_contact_matrix(dir_contacts)
    if save_to is not None:
        np.save(save_to, full_contacts)
    return full_contacts


def count_foci_bonds(dir_contacts: np.ndarray) -> np.ndarray:
    """
    Compute the number and frequency of bonds in a direct contact matrix.

    This function analyzes the direct contact matrix to infer the distribution
    of bonds per atom. Each atom's number of bonds is determined by summing
    the contact matrix along its rows or columns, subtracting the diagonal
    self-bonds. Thus, number of bonds ranges from 0 to :math:`n_{atoms}-1`

    Parameters
    ----------
    dir_contacts : np.ndarray, shape (n_atoms, n_atoms)
        A symmetric binary matrix where element `(i, j)` is 1 if a contact
        exists between atoms `i` and `j`, and 0 otherwise. Diagonal elements
        must be 1, indicating self-bonds.

    Returns
    -------
    bind_list: np.ndarray, shape (n_atoms-1,)
        A 1D array where the index represents bond size (number of bonds per
        atom), and the value is the frequency of that bond size in the system.

    Raises
    ------
    ValueError
         If the input `dir_contacts` matrix is not square.

    Examples
    --------
    >>> dir_contacts = np.array([[1, 1, 0],
    ...                          [1, 1, 1],
    ...                          [0, 1, 1]])
    >>> count_foci_bonds(dir_contacts)
    array([0, 1, 2])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720
    """
    # Ensure the input is a square matrix
    n_atoms, n_cols = dir_contacts.shape

    if n_atoms != n_cols:
        raise ValueError("The direct contact matrix must be square.")

    # Sum along rows to count bonds for each atom and exclude self-bonds
    bonds = dir_contacts.sum(axis=1) - 1  # Removing diagonal self-bonds

    # Create a histogram of bond sizes
    bond_sizes, bond_counts = np.unique(bonds, return_counts=True)
    bond_list = np.zeros(n_atoms, dtype=int)  # Array for bond frequencies
    np.put(bond_list, bond_sizes, bond_counts)
    return bond_list


def is_symmetric(matrix: np.ndarray) -> bool:
    """
    Check if a given square matrix is symmetric.

    A matrix is symmetric if it is equal to its transpose.

    Parameters
    ----------
    matrix : np.ndarray
        The matrix to be checked for symmetry.

    Returns
    -------
    bool
        True if the matrix is symmetric, False otherwise.

    Raises
    ------
    ValueError
        If the input matrix is not square.

    Notes
    -----
    This function uses NumPy's `allclose` method to account for numerical
    tolerances in floating-point  comparisons.

    Examples
    --------
    >>> import numpy as np
    >>> mat = np.array([[1, 2], [2, 1]])
    >>> is_symmetric(mat)
    True

    >>> mat = np.array([[1, 0], [2, 1]])
    >>> is_symmetric(mat)
    False
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("The input matrix must be square.")
    return np.allclose(matrix, matrix.T)


def is_positive_semi_definite(matrix: np.ndarray) -> bool:
    """
    Check if a given square matrix is positive semi-definite.

    A matrix is positive semi-definite if all its eigenvalues are non-negative.

    Parameters
    ----------
    matrix : np.ndarray
        The matrix to be checked for positive semi-definiteness.

    Returns
    -------
    bool
        True if the matrix is positive semi-definite, False otherwise.

    Raises
    ------
    ValueError
        If the input matrix is not square.

    Notes
    -----
    This function attempts a Cholesky decomposition, which only succeeds for
    positive semi-definite matrices.

    Examples
    --------
    >>> import numpy as np
    >>> mat = np.array([[2, -1], [-1, 2]])
    >>> is_positive_semi_definite(mat)
    True

    >>> mat = np.array([[1, 2], [2, -3]])
    >>> is_positive_semi_definite(mat)
    False

    References
    ----------
    Higham, N. J. (2002). Accuracy and Stability of Numerical Algorithms
    (2nd ed.). SIAM.
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("The input matrix must be square.")
    try:
        _ = np.linalg.cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False


def count_foci_clusters(contacts: np.ndarray) -> np.ndarray:
    """
    Infer the number and size of clusters from a contact matrix.

    This function calculates the eigenvalues of the provided contact matrix.
    The number of non-zero eigenvalues corresponds to the number of clusters,
    and their values indicate the sizes of these clusters.

    Before computing eigenvalues, the function checks if the contact matrix is
    symmetric and positive semi-definite. Eigenvalues with negligible imaginary
    parts are treated as real.

    Parameters
    ----------
    contacts : np.ndarray
        A binary square matrix representing contacts between atoms. Each
        element is either 0 or 1, indicating the absence or presence of a
        contact, respectively. Diagonal elements are assumed to be 1,

    Returns
    -------
    cluster_list: np.ndarray, shape (cluster_size-1,)
        A 1D array where each index represents :math:`cluster_size - 1`, and
        the values are the frequency of clusters of that size.

    Raises
    ------
    ValueError
        If the input matrix is not symmetric, not positive semi-definite,
        or contains invalid cluster sizes.

    Notes
    -----
    An eigenvalue of 0 does not have a direct physical meaning; however, the
    frequency of this eigenvalue corresponds to the total number of repeated
    clusters of size 2 or larger. In terms of linear algebra, this is
    equivalent to the number of dependent columns in the contact matrix.

    Examples
    --------
    >>> contacts = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])
    >>> count_foci_clusters(contacts)
    array([1, 1, 1])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720
    """
    n_atoms, _ = contacts.shape
    # Compute eigenvalues and round for stability
    clusters_row = npla.eigvalsh(contacts)
    clusters = np.asarray(np.round(clusters_row), dtype=int)

    # Validate eigenvalues for physical consistency
    if np.any(clusters < 0):
        raise ValueError("A cluster with size smaller than 0 found!")

    # Count occurrences of each cluster size
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    cluster_list = np.zeros(n_atoms + 1, dtype=int)
    np.put(cluster_list, cluster_sizes, cluster_counts)

    return cluster_list


def count_foci_clusters_old(contacts: np.ndarray) -> np.ndarray:
    """
    Infer the number and size of clusters from a contact matrix.

    This function calculates eigenvalues of the contact matrix, where the
    number of non-zero eigenvalues corresponds to the number of clusters. The
    values of the eigenvalues represent the sizes of clusters.

    Parameters
    ----------
    contacts : np.ndarray
        A binary square matrix representing contacts between atoms. Diagonal
        elements must be 1.

    Returns
    -------
    cluster_list: np.ndarray, shape (cluster_size-1,)
        A 1D array where each index represents :math:`cluster_size - 1`, and
        the values are the frequency of clusters of that size.
    Raises
    ------
    ValueError
        If the input matrix contains invalid eigenvalues.

    Notes
    -----
    An eigenvalue of 0 does not have a direct physical meaning; however, the
    frequency of this eigenvalue corresponds to the total number of repeated
    clusters of size 2 or larger. In terms of linear algebra, this is
    equivalent to the number of dependent columns in the contact matrix.

    Examples
    --------
    >>> contacts = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])
    >>> count_foci_clusters_old(contacts)
    array([1, 1, 1])

    References
    ----------
    Sevick, E. M., Monson, P. A., & Ottino, J. M. (1988). Monte Carlo
    calculations of cluster statistics in continuum models of composite
    morphology. The Journal of Chemical Physics, 89(1), 668-676.
    https://doi.org/10.1063/1.454720
    """
    n_atoms, _ = contacts.shape

    # Compute eigenvalues and round for numerical stability
    clusters = npla.eigvalsh(contacts)

    # Validate eigenvalues
    if np.any(clusters < 0):
        raise ValueError("A cluster with size smaller than 0 found!")

    # Count occurrences of each cluster size
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    cluster_list = np.zeros(n_atoms + 1, dtype=int)
    np.put(cluster_list, cluster_sizes, cluster_counts)

    return cluster_list


def foci_info(
    genomic_pos: np.ndarray,
    nmon_large: int,
    nmon_small: int
) -> np.ndarray:
    """
    Generate information about all possible pairs of large monomers in a
    circular polymer.

    This function calculates the genomic indices and distances for all pairs of
    large monomers in a circular heterogeneous polymer composed of two types of
    monomers: small and large. The genomic distance is defined as the minimum
    distance along the ring topology.

    Parameters
    ----------
    genomic_pos : np.ndarray
        A 1D array of size `nmon_large` containing the genomic positions of
        large monomers.
    nmon_large : int
        Number of large monomers.
    nmon_small : int
        Number of small monomers.

    Returns
    -------
    pair_info: np.ndarray
        An array of shape `(n_pairs, 3)` where:
        - The first two columns represent the indices of monomer pairs along
          the backbone.
        - The third column contains the genomic distance between the pair.

    Notes
    -----
    The genomic distance accounts for the circular topology of the polymer. It
    is the minimum of the direct distance and the distance around the circle.

    Examples
    --------
    >>> genomic_pos = np.array([0, 5, 10, 15])
    >>> foci_info(genomic_pos, nmon_large=4, nmon_small=6)
    array([[ 0,  5,  4],
           [ 0, 10,  3],
           [ 0, 15,  2],
           [ 5, 10,  4],
           [ 5, 15,  3],
           [10, 15,  4]])
    """
    nmon = nmon_small + nmon_large

    # Generate all combinations of large monomer pairs
    pairs = np.array(list(combinations(range(nmon_large), 2)), dtype=np.int_)

    # Retrieve genomic positions for each pair
    # genomic_idx = np.choose(pairs, genomic_pos)
    # to-do: check wither choose is better or take.
    genomic_idx = np.take(genomic_pos, pairs)

    # Calculate direct genomic distances
    genomic_dist = np.abs(np.diff(genomic_idx, axis=1))
    # Applying the topological contain on the genomic distance, i.e. finding
    # minimum distance between a pair on a circle. The pair genomic distance is
    # inclusive, so their difference is one unit more than the actual number of
    # genomic units between the pair, but it is equal to the number of bonds
    # (not the total length of the bonds).

    # Adjust distances for circular topology
    # Number of small monomers between to large monomers:
    genomic_dist = np.minimum.reduce((genomic_dist, nmon - genomic_dist)) - 1

    # Combine pair indices and genomic distances
    # pair_info = np.concatenate((genomic_idx, genomic_dist), axis=1)
    pair_info = np.hstack((genomic_idx, genomic_dist))

    return pair_info


def v_histogram(
    data: np.ndarray,
    bins: int,
    bin_range: Tuple[float, float],
    **kwargs
) -> np.ndarray:
    """
    Compute the histogram of input data with specified bins and range.

    This function wraps NumPy's `np.histogram` to compute histograms of data
    with user-defined parameters. It is designed for vectorized operations
    across multiple datasets.

    Parameters
    ----------
    data : np.ndarray
        1D array-like input data to compute the histogram for.
    bins : int
        Number of bins to divide the data into.
    bin_range : Tuple[float, float]
        The lower and upper range of bins.

    Returns
    -------
    np.ndarray
        An array containing the frequencies of `data` in each bin.

    Examples
    --------
    >>> data = np.array([0.5, 1.5, 2.5, 3.5])
    >>> v_histogram(data, bins=3, bin_range=(0, 3))
    array([1, 1, 2])
    """
    hist, _ = np.histogram(data, bins=bins, range=bin_range, **kwargs)
    return hist


def foci_rdf(
    hist: np.ndarray,
    bin_centers: np.ndarray,
    binsize: float
) -> np.ndarray:
    """
    Compute the radial distribution function (RDF).

    The RDF provides the probability of finding a particle at a given radial
    distance relative to another particle. It is calculated by normalizing
    particle frequencies in each radial bin with respect to the shell volume.

    Parameters
    ----------
    hist : np.ndarray
        An array containing the frequencies of particles in each radial bin.
    bin_centers : np.ndarray
        An array containing the positions of bin centers along the radial
        direction.
    binsize : float
        The width of each radial bin.

    Returns
    -------
    np.ndarray
        The computed RDF values for each radial bin.

    Notes
    -----
    The RDF is normalized such that the integral over all bins equals 1.

    Examples
    --------
    >>> hist = np.array([10, 20, 30])
    >>> bin_centers = np.array([1, 2, 3])
    >>> binsize = 0.5
    >>> foci_rdf(hist, bin_centers, binsize)
    array([0.31830989, 0.15915494, 0.1061033 ])
    """
    shell_volumes = 4 * np.pi * bin_centers**2 * binsize
    rdf = hist / shell_volumes
    rdf /= np.sum(rdf)  # Normalize RDF
    return rdf


def foci_histogram(
    tags: List[str],
    pairs_dist: np.ndarray,
    binsize: float
) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    Generate histograms and RDFs for pairwise distances.

    This function computes histograms and radial distribution functions (RDFs)
    for all pairwise distances of large monomers along a heterogeneous polymer.
    It uses the same number of bins and bin range for all pairs. The lower bin
    range is 0 while the upper one is the maximum of physical distances between
    any pairs.

    Parameters
    ----------
    tags : List[str]
        A list of strings containing pair names.
    pairs_dist : np.ndarray
        An array of shape `(n_pairs, n_frames)` containing the time series of
        pairwise distances for each pair.
    binsize : float
        The width of each histogram bin.

    Returns
    -------
    Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]
        - `pair_hists`: A dictionary where keys are pair tags, and values are
          histograms of pair distances. Includes a special key 'bin_center'
          for bin centers.
        - `pair_rdfs`: A dictionary where keys are pair tags, and values are
          RDFs of pair distances. Includes a special key 'bin_center' for bin
          centers.

    Examples
    --------
    >>> tags = ['pair1', 'pair2']
    >>> pairs_info = np.array([[0, 1, 2], [1, 2, 3]])
    >>> pairs_dist = np.array([[1.0, 1.5], [2.0, 2.5]])
    >>> binsize = 0.5
    >>> foci_histogram(tags, pairs_info, pairs_dist, binsize)
    ({'pairDistHistFoci-pair1': array([...]), ...},
     {'pairDistRdfFoci-pair1': array([...]), ...})
    """
    # Determine histogram range and bin centers
    # the maximum pair distance used as upper limit of histogram range in all
    # pair histogram. This value is rounded to be a multiple of binsize.
    upper_range = np.ceil(np.max(pairs_dist) / binsize) * binsize
    max_range = (0, upper_range)
    max_bins = int(np.ceil((max_range[1] - max_range[0]) / binsize))
    bin_edges = np.arange(max_range[0], max_range[1] + binsize, binsize)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

    # Compute histograms for all pairs
    hists = [
        np.histogram(data, bins=max_bins, bin_range=max_range)
        for data in pairs_dist
    ]
    hist_tags = [f"pairDistHistFoci-{tag}" for tag in tags]
    pair_hists = dict(zip(hist_tags, hists))
    pair_hists['bin_center'] = bin_centers

    # Compute RDFs for all pairs
    rdfs = [
        foci_rdf(hist, bin_centers, binsize)
        for hist in hists
    ]
    rdf_tags = [f"pairDistRdfFoci-{tag}" for tag in tags]
    pair_rdfs = dict(zip(rdf_tags, rdfs))
    pair_rdfs['bin_center'] = bin_centers

    return pair_hists, pair_rdfs


def whole_dist_mat_foci(
    whole_path: str,
    whole_info: Union[TransFociCub, TransFociCyl, 
                      SumRuleCubHeteroRing, SumRuleCubHeteroLinear]
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate time series, histograms, and RDFs of foci distances from a
    distance matrix.

    This function processes the distance matrix of large monomers (foci) in a
    polymer system. It computes time series of pairwise distances, histograms,
    and radial distribution functions (RDFs) for all foci pairs.

    Parameters
    ----------
    whole_path : str
        Path to the distance matrix file. The matrix should be a 3D array where
        the first axis represents time frames.
    whole_info : Union[TransFociCub, TransFociCyl,
                       SumRuleCubHeteroRing, SumRuleCubHeteroLinear]
        Object containing metadata about the simulation, including the number
        of large and small monomers and monomer properties.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        - `pair_hists`: DataFrame where columns are pair tags, and values are
          pair histograms.
        - `pair_rdfs`: DataFrame where columns are pair tags, and values are
          pair RDFs.
        - `pair_tseries`: DataFrame where columns are pair tags, and values are
          time series of pair distances.

    Examples
    --------
    >>> whole_path = "distance_matrix.npy"
    >>> whole_info = TransFociCub(dmon_small=1.0, nmon_large=4, nmon_small=6)
    >>> pair_hists, pair_rdfs, pair_tseries = \
            whole_dist_mat_foci(whole_path, whole_info)

    Notes
    -----
    The distance matrix is assumed to be upper-triangular for efficiency. Tags
    are generated in the format "locXlocYgenDistZ" where X and Y are monomer
    indices, and Z is their genomic distance.
    """
    # Load the distance matrix
    whole = np.load(whole_path)

    # Extract parameters from whole_info
    dmon_small = whole_info.dmon_small
    nmon_large = whole_info.nmon_large
    nmon_small = whole_info.nmon_small

    # Extract genomic positions and compute pair information
    diagonal_idxs = np.diag_indices(nmon_large)
    genomic_pos = whole[0, diagonal_idxs[0], diagonal_idxs[1]]
    pairs_info = foci_info(genomic_pos, nmon_large, nmon_small)

    # Generate tags for each pair
    tags = [
        f"loc{x[0]}loc{x[1]}genDist{x[2]}"
        for x in pairs_info[:, :3]
    ]

    # Extract upper-triangular indices for distances
    triu_indices = np.triu_indices(nmon_large, k=1)
    pairs_dist = whole[:, triu_indices[0], triu_indices[1]].T

    # Compute histograms and RDFs
    binsize = 0.5 * dmon_small
    pair_hists, pair_rdfs = foci_histogram(tags, pairs_dist, binsize)

    # Convert results to DataFrames
    pair_hists = pd.DataFrame.from_dict(pair_hists)
    pair_rdfs = pd.DataFrame.from_dict(pair_rdfs)

    # Time series DataFrame
    tseries_tags = [f"pairDistTFoci-{tag}" for tag in tags]
    pair_tseries = pd.DataFrame(dict(zip(tseries_tags, pairs_dist)))

    return pair_hists, pair_rdfs, pair_tseries


def enforce_single_patch_dir_contact(
    contact_matrix: np.ndarray
) -> np.ndarray:
    """
    Enforces single-patch contact for each column in a contact matrix.

    This function modifies a binary contact matrix to ensure that each patch
    (column) is only in contact with one monomer (row). If a patch is in
    contact with multiple monomers, only the first (topmost) contact is
    preserved, and the rest are set to zero.

    Parameters
    ----------
    contact_matrix : np.ndarray
        Binary contact matrix of shape `(n_monomers, n_patches)`, where 1
        indicates a direct contact between a monomer and a patch.

    Returns
    -------
    contact_matrix: np.ndarray
        Modified contact matrix where all values after the first contact in
        each column are set to 0.

    Raises
    ------
    ValueError
        If the input matrix is not 2D.

    Notes
    -----
    Each H-NS patch can be in contact with one monomer. However, visual
    inspection has revealed that a single patch can be in contact with more
    than one monomers, which is an artifacts and does not have physical
    grounds. The inspection has revealed these monomers are adjacent monomers
    along the chain backbone. Therefore, it is decided to only consider the
    contact of the H-NS patch with the first of these monomers as a real
    contact and conisder the rest of contacts as artifacts. This means setting
    the non-zero elements below the first non-zero elements along a column
    equal to zero.

    Examples
    --------
    >>> contact_matrix = np.array([[1, 0, 1],
    ...                            [1, 1, 0],
    ...                            [0, 1, 0]])
    >>> enforce_single_patch_dir_contact(contact_matrix)
    array([[1, 0, 1],
           [0, 1, 0],
           [0, 0, 0]])
    """
    if contact_matrix.ndim != 2:
        raise ValueError("Input contact matrix must be a 2D array.")

    # Create a copy to avoid modifying the input matrix in place
    contact_matrix = contact_matrix.copy()

    # Iterate over columns to enforce single-patch contact
    for col_idx in range(contact_matrix.shape[1]):
        # Find the first contact index in the column
        first_contact_idx = np.argmax(contact_matrix[:, col_idx] != 0)
        # Zero out all contacts below the first contact
        contact_matrix[first_contact_idx + 1:, col_idx] = 0

    return contact_matrix


def genomic_distance(
    atoms_a: np.ndarray,
    atoms_b: np.ndarray,
    topology: TopologyT,
    num_backbone_atoms: int
) -> np.ndarray:
    """
    Calculate the genomic distance between pairs of atoms along a polymer
    backbone.

    Genomic distance is the intra-monomer distance between atoms along the
    backbone of a polymer. It is defined based on the polymer topology
    ('linear' or 'ring').

    Parameters
    ----------
    atoms_a : np.ndarray
        Array of atom indices in the first group along the backbone, with shape
        `(num_backbone_atoms,)` or `(num_backbone_atoms, 1)`.
    atoms_b : np.ndarray
        Array of atom indices in the second group along the backbone, with
        shape `(num_backbone_atoms,)` or `(num_backbone_atoms, 1)`.
    topology : TopologyT
        Topology of the polymer ('linear' or 'ring').
    num_backbone_atoms : int
        Total number of atoms along the backbone.

    Returns
    -------
    distances: np.ndarray
        Array of element-wise genomic distances between corresponding pairs of
        atoms along the backbone.

    Raises
    ------
    ValueError
        If `atoms_a` and `atoms_b` have different shapes.
        If `topology` is not 'linear' or 'ring'.

    Notes
    -----
    - For 'linear' topology, the genomic distance is the absolute difference
      between atom indices.
    - For 'ring' topology, the genomic distance is the minimum of the direct
      distance and the distance around the circle.
    - Genomic distance is measured in terms of the number of bonds, not atoms.
      It is 1 unit more than the number of atoms between a pair of atoms.

    Examples
    --------
    >>> atoms_a = np.array([0, 2, 4])
    >>> atoms_b = np.array([1, 3, 5])
    >>> genomic_distance(atoms_a, atoms_b, topology='linear',
                         num_backbone_atoms=6)
    array([1, 1, 1])
    >>> genomic_distance(atoms_a, atoms_b, topology='ring',
                         num_backbone_atoms=6)
    array([1, 1, 1])
    """
    if atoms_a.shape != atoms_b.shape:
        raise ValueError("'atoms_a' and 'atoms_b' must have the same shape.")
    if topology not in ('linear', 'ring'):
        raise ValueError(
            f"Invalid topology: '{topology}'. Must be 'linear' or 'ring'.")

    # Define topology-specific distance calculation
    topology_functions = {
        'linear': lambda a, b: np.abs(a - b),
        'ring': lambda a, b, n: np.minimum(np.abs(a - b), n - np.abs(a - b)),
    }
    distances: np.ndarray = \
        topology_functions[topology](atoms_a, atoms_b, num_backbone_atoms)
    return distances


def hns_genomic_distance(
    contact_nonzeros: np.ndarray,
    topology: TopologyT,
    n_mons: int
) -> np.ndarray:
    """
    Calculate genomic distances for monomer pairs bridged by an H-NS protein.

    This function computes the genomic distances for monomer pairs that are
    connected via an H-NS protein in a polymer system. The distances are
    determined based on the polymer's topology ('linear' or 'ring').

    Parameters
    ----------
    contact_nonzeros : np.ndarray
        Array of shape `(n_pairs, 2)` containing indices of non-zero elements
        in the monomer-monomer contact matrix.
    topology : TopologyT
        Topology of the polymer ('linear' or 'ring').
    n_mons : int
        Total number of monomers in the polymer.

    Returns
    -------
    np.ndarray
        A 2D array with shape `(n_pairs, 3)`:
        - Column 0: Index of the first monomer in the pair.
        - Column 1: Index of the second monomer in the pair.
        - Column 2: Genomic distance between the two monomers.

    Raises
    ------
    ValueError
        If `topology` is not 'linear' or 'ring'.

    Examples
    --------
    >>> contact_nonzeros = np.array([[0, 2], [1, 3], [2, 4]])
    >>> hns_genomic_distance(contact_nonzeros, topology='linear', n_mons=5)
    array([[0, 2, 2],
           [1, 3, 2],
           [2, 4, 2]])
    >>> hns_genomic_distance(contact_nonzeros, topology='ring', n_mons=5)
    array([[0, 2, 2],
           [1, 3, 2],
           [2, 4, 2]])
    """
    # Calculate genomic distances for the provided pairs
    genomic_distances = genomic_distance(
        contact_nonzeros[:, 0],
        contact_nonzeros[:, 1],
        topology,
        n_mons
    )

    # Combine pair indices with genomic distances
    contact_m_m_gen_dist = \
        np.column_stack((contact_nonzeros, genomic_distances))

    return contact_m_m_gen_dist


def hns_binding(
    direct_contact: np.ndarray,
    topology: str,
    cis_threshold: int = 4,
    binding_stats: Optional[Dict[str, List[int]]] = None,
    loop_length_hist: Optional[np.ndarray] = None
) -> Tuple[Dict[str, List[int]], np.ndarray]:
    """
    Calculate the binding statistics of H-NS proteins to monomers.

    This function evaluates the binding behavior of H-NS proteins in a polymer
    system. Each H-NS protein has a core and two patches, with patches able to
    bind monomers. The function calculates statistics for free, bridged, and
    dangled H-NS proteins and determines genomic distances for binding.

    Parameters
    ----------
    direct_contact : np.ndarray
        Binary matrix of shape `(n_mon, n_hpatch)`, where each element `(i, j)`
        is 1 if monomer `i` is in contact with patch `j`, and 0 otherwise.
    topology : str
        Topology of the polymer ('linear' or 'ring').
    cis_threshold : int, optional
        Genomic distance threshold (in bonds) below which an H-NS binding is
        classified as cis-binding. Default is 4.
    binding_stats : dict, optional
        A dictionary to store binding statistics. Keys include:
        - 'n_m_hpatch_bound': Number of monomer-patch contacts.
        - 'n_hpatch_engaged': Number of engaged patches.
        - 'n_hpatch_free': Number of free patches.
        - 'n_hcore_free': Number of free H-NS cores.
        - 'n_hcore_bridge': Number of bridging H-NS cores.
        - 'n_hcore_dangle': Number of dangling H-NS cores.
        If None, a new dictionary is initialized.
    loop_length_hist : np.ndarray, optional
        Histogram of loop sizes (genomic distances). If None, a new array is
        initialized based on `topology`.

    Returns
    -------
    Tuple[Dict[str, List[int]], np.ndarray]
        - Updated `binding_stats` dictionary with binding data.
        - Updated `loop_length_hist` array containing the histogram of loop
          sizes.

    Raises
    ------
    ValueError
        If the topology is invalid ('linear' or 'ring' are the only accepted
        values).

    Notes
    -----
    - The function handles rare cases where patches are in contact with
      multiple monomers by using `enforce_single_patch_dir_contact`.
    - Genomic distances are calculated using the H-NS core's bridging matrix.
    - Loop sizes exceeding the `cis_threshold` are classified as trans-binding.

    Examples
    --------
    >>> direct_contact = np.array([[1, 0, 1], [0, 1, 0], [1, 0, 0]])
    >>> topology = 'linear'
    >>> binding_stats, loop_hist = hns_binding(direct_contact, topology)
    """
    # Initialize results if not provided
    if binding_stats is None:
        binding_stats = defaultdict(list)
    if loop_length_hist is None:
        max_gen_distance = {'linear': direct_contact.shape[0], 'ring': direct_contact.shape[0] // 2 + 1}
        loop_length_hist = np.zeros(max_gen_distance[topology], dtype=int)

    n_mon, n_hpatch = direct_contact.shape
    n_hcore = n_hpatch // 2

    if topology not in ('linear', 'ring'):
        raise ValueError(f"Invalid topology: '{topology}'. Must be 'linear' or 'ring'.")

    # Compute basic binding statistics
    binding_stats['n_m_hpatch_bound'].append(np.sum(direct_contact))
    d_cont_per_hpatch = np.sum(direct_contact, axis=0)
    binding_stats['n_hpatch_engaged'].append(np.count_nonzero(d_cont_per_hpatch))
    binding_stats['n_hpatch_free'].append(n_hpatch - binding_stats['n_hpatch_engaged'][-1])

    # Compute core-level binding statistics
    binding_stats['n_hcore_free'].append(np.sum(
        (d_cont_per_hpatch[0::2] == 0) & (d_cont_per_hpatch[1::2] == 0)
    ))
    binding_stats['n_hcore_bridge'].append(np.sum(
        (d_cont_per_hpatch[0::2] > 0) & (d_cont_per_hpatch[1::2] > 0)
    ))
    binding_stats['n_hcore_dangle'].append(
        n_hcore - binding_stats['n_hcore_free'][-1] - binding_stats['n_hcore_bridge'][-1]
    )

    # Enforce single monomer-patch contact and compute monomer-monomer binding
    single_patch_dir_contact = enforce_single_patch_dir_contact(direct_contact)
    cont_m_hpatch = generate_contact_matrix(single_patch_dir_contact)
    cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
    cont_m_m = np.matmul(cont_m_hcore, cont_m_hcore.T)
    cont_m_m_triu = np.triu(cont_m_m, 1)

    # Extract non-zero indices and calculate genomic distances
    cont_m_m_nonzeros = np.array(np.where(cont_m_m_triu > 0)).T
    m_m_gen_dist = hns_genomic_distance(cont_m_m_nonzeros, topology, n_mon)

    # Classify bindings and update loop length histogram
    binding_stats['n_hcore_cis'].append(np.count_nonzero(
        (m_m_gen_dist[:, 2] > 0) & (m_m_gen_dist[:, 2] <= cis_threshold)
    ))
    binding_stats['n_hcore_trans'].append(np.count_nonzero(
        m_m_gen_dist[:, 2] > cis_threshold
    ))
    np.add.at(loop_length_hist, m_m_gen_dist[:, 2], 1)

    return binding_stats, loop_length_hist


def generate_mon_bind_direct(
    mon_patch_dir: np.ndarray,
    patch_per_binder: int
) -> np.ndarray:
    """Creates monomer-binder direct contact matrix from monomer-patch direct
    contact matrix.

    Parameters
    ----------
    mon_patch_direcy : np.ndarray
        Monomer-patch direct contact matrix of size n_monomer *
        (patch_per_binder)*n_binder

    Returns
    -------
    np.ndarray
        Monomer-binder direct contact matrix of size n_monomer * n_binder
    """
    if len(mon_patch_dir.shape) != 2:
        raise ValueError(
            "The monomer-patch direct contact matrix is not a matrix!")
    if mon_patch_dir.shape[1] % patch_per_binder != 0:
        raise ValueError(
            "The number of columns is nut a multiplier of"
            f"'n={patch_per_binder}'"
            )

    n_bind = mon_patch_dir.shape[1] // patch_per_binder
    mon_bind_dir = np.zeros((mon_patch_dir.shape[0], n_bind), dtype=np.int_)

    for i in range(n_bind):
        start_col = patch_per_binder * i
        end_col = start_col + patch_per_binder
        binder_patches = mon_patch_dir[:, start_col:end_col]
        binder = np.any(binder_patches, axis=1).astype(int)
        mon_bind_dir[:, i] = binder

    return mon_bind_dir


def split_binder_matrix(M_ij: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Splits the Monomer-Binder direct contact matrix into two matrices based
    on whether a binder is a bridger or a dangler.

    Parameters
    ----------
    M_ij (numpy.ndarray):
        The monomer-binder direct contact matrix, where each row represents a
        monomer and each column represents a binder. The matrix is binary,
        with '1' indicating direct contact.

    Returns
    -------
    Tuple of two numpy.ndarrays:
        D_ij: The matrix representing dangling binders, where each binder is
        attached to the polymer chain through only one patch.
        B_ij: The matrix representing bridging binders, where each binder
        connects two monomers using both patches.

    The function calculates the sum of each column in M_ij to classify each
    binder as dangling (sum = 1) or bridging (sum = 2) and populates the D_ij
    and B_ij matrices accordingly.
    """
    n_monomers, n_binders = M_ij.shape
    D_ij = np.zeros_like(M_ij)  # Dangling binders matrix
    B_ij = np.zeros_like(M_ij)  # Bridging binders matrix

    # Calculate the sum of each column to classify binders
    column_sums = np.sum(M_ij, axis=0)

    # Populate the D_ij and B_ij matrices
    for binder in range(n_binders):
        if column_sums[binder] == 1:
            D_ij[:, binder] = M_ij[:, binder]
        elif column_sums[binder] == 2:
            B_ij[:, binder] = M_ij[:, binder]

    return D_ij, B_ij


def find_binder_clusters_need_attention(M_ij: np.ndarray, topology) -> np.ndarray:
    """Identifies clusters of binders in a polymer system based on their
    binding to the same or adjacent monomers.

    Parameters
    ----------
    M_ij (numpy.ndarray)
        Monomer-binder direct contact matrix, where rows represent monomers and
        columns represent binders. It is a binary matrix with '1' indicating
        direct contact.

    Returns:
    numpy.ndarray:
        A square binary matrix where each element (i, j) is 1 if binders i and
        j are part of the same cluster, and 0 otherwise. This matrix is
        symmetric.

    The function sets all diagonal elements of the output matrix to 1,
    indicating each binder forms a cluster with itself. It then iterates over
    all pairs of binders, checking for direct connections to the same monomer
    or connections via adjacent monomers.
    """
    n_binders = M_ij.shape[1]
    n_monomers = M_ij.shape[0]
    B_ij = np.zeros((n_binders, n_binders), dtype=int)
    print('need attention: see comments')
    # diagonal elements have meaning here; if it is 0 means unbound hns/binder
    # but if 1 means it bridged two momnomers or dangle from on monomer.
    np.fill_diagonal(B_ij, 1)
    # Check for binders connected to the same monomer
    for i in range(n_binders):
        for j in range(i + 1, n_binders):
            if np.any(M_ij[:, i] & M_ij[:, j]):
                B_ij[i, j] = B_ij[j, i] = 1

    # Check for binders connected to adjacent monomers
    for monomer in range(n_monomers - 1):
        adjacent_binders = np.where(M_ij[monomer] | M_ij[monomer + 1])[0]
        for binder in adjacent_binders:
            B_ij[binder, adjacent_binders] = 1

    if topology == 'ring':
        adjacent_binders = np.where(M_ij[0] | M_ij[-1])[0]
        for binder in adjacent_binders:
            B_ij[binder, adjacent_binders] = 1

    return B_ij
