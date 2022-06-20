"""A sub-module to detect and analyze particle clusters in a given system. It
worth noting there is a difference analyzing particle clusters and
trajectory (or ensemble) clusters.
Below, different algorithms found in literature are implemented in scientific
Python.
"""
from typing import Tuple, Optional
from itertools import combinations
import numpy as np
import numpy.linalg as npla


def apply_pbc(
    pbc_lengths: np.ndarray,
    pbc_lengths_inverse: np.ndarray,
    pbc: dict
) -> Tuple[np.ndarray, np.ndarray]:
    """Updates the periodic boundary condition (PBC) information in the
    dimentions (keys) by the length (values) passed by `pbc`. The key
    are indexes of the `pbc_lengths` and `pbc_lengths_inverse` arrays.

    Parameters
    ----------
    pbc_lengths: numpy.ndarray
        An array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        An array of the invesres of lengths/sizes in different directions.
    pbc: dict
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Return
    ------
    pbc_lengths: numpy.ndarray
        The updated array of lengths/sizes in different directions.
    pbc_lengths_inverse: numpy.ndarray
        The updated array of the invesres of lengths/sizes in different
        directions.
    """
    for dim, length in pbc.items():
        pbc_lengths[dim] = length
        pbc_lengths_inverse[dim] = 1 / length
    return pbc_lengths, pbc_lengths_inverse


def dist_sq_matrix(
    positions: np.ndarray,
    pbc: Optional[dict] = None
) -> np.ndarray:
    """Finds the squares of distances between all the pair of atoms
    (even self-distances where i=j) for a given list of atoms' `positions`.

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different dimentions
        or components.

    pbc: dict, defualt None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Rutern
    ------
    dist_sq: np.ndarray
        A square matrix of size n_atoms * n_atoms in which each elements are
        the squares of the center-to-center distances between different pairs
        of atoms. All the diagonal elements of this materix is 0.
    """
    n_atoms, n_dims = positions.shape
    # Defining pbc lengths with value zero in all directions:
    # The approach allows us to combine linear agebra and numpy broadcasting
    # and efficently apply pbc.
    pbc_lengths = np.zeros(n_dims)
    pbc_lengths_inverse = np.zeros(n_dims)
    if pbc is not None:
        pbc_lengths, pbc_lengths_inverse = apply_pbc(
            pbc_lengths, pbc_lengths_inverse, pbc
        )
    # Redefining the shapres ofdifferent arrays to combine linear algebra
    # with numpy broadcasting.
    pbc_lengths = np.reshape(pbc_lengths, (n_dims, 1, 1))
    pbc_lengths_inverse = np.reshape(pbc_lengths_inverse, (n_dims, 1, 1))
    pos_T = positions.T
    pos_j_element = np.reshape(pos_T, (n_dims, n_atoms, 1))
    pos_i_element = np.reshape(pos_T, (n_dims, 1, n_atoms))
    dist_sq = pos_j_element - pos_i_element
    dist_sq = dist_sq - pbc_lengths * np.around(pbc_lengths_inverse * dist_sq)
    dist_sq = dist_sq ** 2
    dist_sq = np.sum(dist_sq, axis=0)
    return dist_sq


def direct_connectivity(dist_sq: np.ndarray, cutoff: float) -> np.ndarray:
    """Finds direct contancts between pairs of atoms.

    Direct connectivity refers to a conditon in which the center-to-center
    distance between two atoms is equal to or less than the `cutoff` distance.
    This definition means that an atom is in direct connectivity with itself.
    Such a definition is correct since clusters with size 1 exist (i. e. a
    signle atom that is not clustered with other atoms); moreover, this
    definition esnures there is a correct number of atoms in a cluster.

    Parameters
    ----------
    dist_sq: np.ndarray
        A square matrix of size n_atoms * n_atoms in which each elements are
        the squares of the center-to-center distances between different pairs
        of atoms. All the diagonal elements of this materix is 0.

    cutoff: float
        The cutoff center-to-center distance for defining a direct
        connectivity.

    Return
    ------
    np.ndarray:
        A symmeteric matrix of size n_atoms * n_atoms in which each element is
        either 0 or 1 and indicates whether an direct connectivity is between a
        pair of atoms or not. 1 means direct connectivity exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition given above.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    return np.asarray(dist_sq <= (cutoff * cutoff), dtype=int)


def connectivity_matrix(dir_connectivites: np.ndarray) -> np.ndarray:
    """Finds the direct and indirect contants between atoms in a system
    from the direct connectivites between these atoms `dir_connectivites`.

    Parameters
    ----------
    dir_connectivites: np.ndarray
        A symmeteric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether an direct connectivity is
        between a pair of atoms or not. 1 means direct connectivity exists
        while 0 means it does not. All the diagonal elements of this matrix
        are 1 by the definition given above.

    Return
    ------
    connectivites: np.ndarray
        A block diagonal matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a connectivity is between a pair
        of atoms or not. 1 means a connectivity exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the definition
        in `direct_connectivites` function.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = dir_connectivites.shape
    col_pairs = list(combinations(range(n_atoms), 2))
    connectivites = dir_connectivites.copy()
    # TODO Is this approach good or converting dir_connectivites to
    # connectivity_matrix?
    for (i, j) in col_pairs:
        if np.any(connectivites[:, i] & connectivites[:, j]):
            transient = connectivites[:, i] | connectivites[:, j]
            connectivites[:, i] = transient
            connectivites[:, j] = transient
    return connectivites


def direct_connectivites_from_sq_dists(
    positions: np.ndarray,
    cutoff: float,
    pbc: Optional[dict] = None,
) -> np.ndarray:
    """Creates the direct connectivity matrix for a group of atoms using a
    cutoff for the center-to-center distance as the ceritrion for clustering.

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different
        dimentions or components.

    cutoff: float
        The cutoff center-to-center distance for defining a direct
        connectivity.

    pbc: dict, default None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Return
    ------
    connectivites: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a connectivity is between a
        pair of atoms or not. 1 means a connectivity exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_connectivites` function.
    """
    dist_sq_mat = dist_sq_matrix(positions, pbc)
    dir_connectivites = direct_connectivity(dist_sq_mat, cutoff)
    connectivites = connectivity_matrix(dir_connectivites)
    return connectivites


def random_direct_connectivites(
    n_atoms: int,
    seed: Optional[int] = 0,
    save_to: Optional[Tuple[str, str]] = None
) -> np.ndarray:
    """Generates a direct-connectivity matrix with size `n_atoms` which can be
    used to test a clustering algorithm. See `direct_connectivity`function
    description for more info.

    A connectivity matrix is symmetric and all its diagonal elements are 1.

    Parameters
    ----------
    n_atoms: int
        Number of atoms or nodes.
    seed: int, default 0
        The seed for the random number generator.
    save_to: tuple of str, default None
        A tumple of two strings in which the first element is an/a
        absolute/relative path of a directory to which the random direct
        connectivity matrix is saved and the second one is the matrix name.

    Return
    ------
    connectivites: np.ndarray
        A symmeteric matrix of size `n_atoms` * `n_atoms` in which each
        element is either 0 or 1 and indicates whether an direct connectivity
        is between a pair of atoms or not. 1 means direct connectivity exists
        while 0 means it does not. All the diagonal elements of this matrix
        are 1 by the definition given in `direct_connectivity` function
        description.
    """
    n_atoms = 5
    # Creating direct-connectivites matrix
    np.random(seed)
    connectivites = np.random.binomial(n=1, p=0.35, size=[n_atoms * n_atoms])
    connectivites = np.reshape(connectivites, (n_atoms, n_atoms))
    connectivites_up = np.triu(connectivites)  # upper triangle
    connectivites_diag = np.identity(n_atoms, dtype=int)
    # dropping diagnal elements:
    connectivites_up = connectivites_up - connectivites_diag
    # ensuring all diagonal elements are either 0 or 1:
    connectivites_up[connectivites_up < 0] = 0
    connectivites_low = connectivites_up.T  # lower traingle: symmetry!
    connectivites = connectivites_up + connectivites_low + connectivites_diag
    # Creating full-connectivites matrix
    connectivites = connectivity_matrix(connectivites)
    if save_to is not None:
        output = save_to[0] + \
            f"random_direct_connectivity-natoms{n_atoms}" + save_to[1]
        np.save(output, connectivites)
    return connectivites


def count_clusters(connectivities: np.ndarray):
    """Infers the number and size of clusters from the connectivity matrix.

    The connectivity matrix is a symmetric and its diagonal elements are
    all 1. While the number of non-zero eigenvalues of the connectivity
    matrix is equla to the number of clusters, the value of a eigenvalue
    shows the size of a cluster. The frequency of eigenvalues is equal to
    the frequency of clusters with different sizes.

    Parameters
    ----------
    connectivities: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a connectivity is between a
        pair of atoms or not. 1 means a connectivity exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_connectivites` function.

    Return
    ------
    cluster_list: np.ndarray
        A 1D array in which the indexes are "cluster_size - 1" and the values
        are frequency of clusters with that size.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = connectivities.shape
    cluster_list = np.zeros(n_atoms, dtype=int)
    clusters, _ = npla.eig(connectivities)
    # numpy returns infinitesimal numbers instead of zero when an eigenvalue
    # is 0. Below, these eigenvalues are rounded to 0. The non-zero
    # values are already integer, so rounding does not affect them.
    clusters = np.asarray(np.round(clusters), dtype=int)
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    # Finding the index of eigenvalue 0 in the cluster_sizes and deleting
    # it from cluster_sizes and its frequency from cluster_counts (
    # Accroding to numpy doc, the indexes are the same in both arrays)
    index_of_zero = np.where(cluster_sizes == 0)
    index_of_zero = index_of_zero[0][0]
    cluster_sizes = np.delete(cluster_sizes, index_of_zero)
    cluster_counts = np.delete(cluster_counts, index_of_zero)
    for size, count in zip(cluster_sizes, cluster_counts):
        cluster_list[size-1] = count
    return cluster_list
