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


def direct_contact(dist_sq: np.ndarray, cutoff: float) -> np.ndarray:
    """Finds direct contancts between pairs of atoms.

    Direct contact refers to a conditon in which the center-to-center
    distance between two atoms is equal to or less than the `cutoff` distance.
    This definition means that an atom is in direct contact with itself.
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
        contact.

    Return
    ------
    np.ndarray:
        A symmeteric matrix of size n_atoms * n_atoms in which each element is
        either 0 or 1 and indicates whether an direct contact is between a
        pair of atoms or not. 1 means direct contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition given above.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    return np.asarray(dist_sq <= (cutoff * cutoff), dtype=int)


def contact_matrix(dir_contacts: np.ndarray) -> np.ndarray:
    """Finds the direct and indirect contants between atoms in a system
    from the direct contacts between these atoms `dir_contacts`.

    Parameters
    ----------
    dir_contacts: np.ndarray
        A symmeteric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether an direct contact is
        between a pair of atoms or not. 1 means direct contact exists
        while 0 means it does not. All the diagonal elements of this matrix
        are 1 by the definition given above.

    Return
    ------
    contacts: np.ndarray
        A block diagonal matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a pair
        of atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the definition
        in `direct_contacts` function.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = dir_contacts.shape
    col_pairs = list(combinations(range(n_atoms), 2))
    contacts = dir_contacts.copy()
    # TODO Is this approach good or converting dir_contacts to
    # contact_matrix?
    for (i, j) in col_pairs:
        if np.any(contacts[:, i] & contacts[:, j]):
            transient = contacts[:, i] | contacts[:, j]
            contacts[:, i] = transient
            contacts[:, j] = transient
    return contacts


def contact_from_sq_dists(
    positions: np.ndarray,
    cutoff: float,
    pbc: Optional[dict] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Creates the direct contact matrix for a group of atoms using a
    cutoff for the center-to-center distance as the ceritrion for clustering.

    All the diagonal elemnts of the contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many
    repeated columns.

    Parameters
    ----------
    positions: np.ndarray
        a matrix with size n_atoms * n_dims of atoms' positions in which
        the rows are atoms' positions and the columns are different
        dimentions or components.

    cutoff: float
        The cutoff center-to-center distance for defining a direct
        contact.

    pbc: dict, default None
        A dictionary of dimensions (keys) and lenghts (values) in which
        the pbc exists.

    Return
    ------
    dir_contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    contacts: np.ndarray
        A matrix of size n_atoms * n_atoms in which each element is either
        0 or 1 and indicates whether a contact is between a pair of
        atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    foci_pair_dist: np.ndarray
        An upper triangluar matrix with 0 zero diagonal elements that
        contains the center-to-center distance between all the pairs of
        atoms.
    """
    dist_sq_mat = dist_sq_matrix(positions, pbc)
    dir_contacts = direct_contact(dist_sq_mat, cutoff)
    contacts = contact_matrix(dir_contacts)
    return dir_contacts, contacts


def random_direct_contacts(
    n_atoms: int,
    seed: Optional[int] = 0,
    save_to: Optional[Tuple[str, str]] = None
) -> np.ndarray:
    """Generates a direct-contact matrix with size `n_atoms` which can be
    used to test a clustering algorithm. See `direct_contact`function
    description for more info.

    All the diagonal elemnts of the direct contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many
    repeated columns.

    Note
    ----
    This function can be used for testing the `count_bonds` and
    `count_clusters` methods.

    Parameters
    ----------
    n_atoms: int
        Number of atoms or nodes.
    seed: int, default 0
        The seed for the random number generator.
    save_to: tuple of str, default None
        A tumple of two strings in which the first element is an/a
        absolute/relative path of a directory to which the random direct
        contact matrix is saved and the second one is the matrix name.

    Return
    ------
    dir_contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.
    """
    n_atoms = 5
    # Creating direct-contacts matrix
    np.random(seed)
    dir_contacts = np.random.binomial(n=1, p=0.35, size=[n_atoms * n_atoms])
    dir_contacts = np.reshape(dir_contacts, (n_atoms, n_atoms))
    dir_contacts_up = np.triu(dir_contacts)  # upper triangle
    dir_contacts_diag = np.identity(n_atoms, dtype=int)
    # dropping diagnal elements:
    dir_contacts_up = dir_contacts_up - dir_contacts_diag
    # ensuring all diagonal elements are either 0 or 1:
    dir_contacts_up[dir_contacts_up < 0] = 0
    dir_contacts_low = dir_contacts_up.T  # lower traingle: symmetry!
    dir_contacts = dir_contacts_up + dir_contacts_low + dir_contacts_diag
    # Creating the full direct contacts matrix
    dir_contacts = contact_matrix(dir_contacts)
    if save_to is not None:
        output = save_to[0] + \
            f"random_direct_contact-natoms{n_atoms}" + save_to[1]
        np.save(output, dir_contacts)
    return dir_contacts


def count_bonds(dir_contacts: np.ndarray):
    """Infers the number and size of bonds from the direct contact matrix
    `dir_contacts`.

    The direct contact matrix is a symmetric and its diagonal elements are
    all 1. The column- or row-wise sum of the matrix shows the number of bonds
    an atom has with itself and other atoms. A self-bond happens due to the
    definition of the direct contact matirx. With this definition, the
    number of bonds per atom ranges from 1 to the number of atoms.

    Parameters
    ----------
    contacts: np.ndarray
        A symmetric matrix of size n_atoms * n_atoms in which each element
        is either 0 or 1 and indicates whether a contact is between a
        pair of atoms or not. 1 means a contact exists while 0 means
        it does not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    Return
    ------
    cluster_list: array-like
        A 1D array in which the indexes are "bond" sizes and the values
        are frequency of bonds with that size.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = dir_contacts.shape
    bond_list = np.zeros(n_atoms, dtype=int)
    # removing self-bonds:
    bonds = dir_contacts.sum(axis=1) - 1
    bond_sizes, bond_counts = np.unique(bonds, return_counts=True)
    np.put(bond_list, bond_sizes, bond_counts)
    return bond_list


def count_clusters(contacts: np.ndarray) -> np.ndarray:
    """Infers the number and size of clusters from the contact matrix
    `contacts`.

    All the diagonal elemnts of the direct contact matrix are 1; moreover,
    depending on the size and number of clusters, the matrix has many repeated
    columns. While the number of non-zero eigenvalues of the contact
    matrix is equal to the number of clusters, the values of eigenvalues
    show the sizes of clusters. The frequencies of eigenvalues are equal to
    the frequencies of clusters with different sizes. An eigenvalue with value
    0 does not have any physical meaning, but the frequency of this eigenvalue
    shows the total number of repeated clusters with sizes equal to or larger
    than 2 (that is equal to the number of dependant columns in the
    contact matrix in the linear algebra lingo).

    Parameters
    ----------
    contacts: np.ndarray
        A matrix of size n_atoms * n_atoms in which each element is either
        0 or 1 and indicates whether a contact is between a pair of
        atoms or not. 1 means a contact exists while 0 means it does
        not. All the diagonal elements of this matrix are 1 by the
        definition in `direct_contacts` function.

    Return
    ------
    cluster_list: array-like
        A 1D array in which the indexes are "cluster_size - 1" sizes and the
        values are frequency of clusters with that size.

    Reference
    ---------
    https://aip.scitation.org/doi/10.1063/1.454720
    """
    n_atoms, _ = contacts.shape
    # While Tte eigenvalue 0 is not physically important, the eigenvaue n_atoms
    # is. Here, a cluster list with values ranging from 0 to n_atoms is created
    # to properly collect all the eigenvalues.
    cluster_list = np.zeros(n_atoms + 1, dtype=int)
    clusters, _ = npla.eig(contacts)
    # Sanity check for having meaning cluster sizes:
    if np.any(clusters < 0):
        raise ValueError("A cluster with size smaller than 0 found!")
    # numpy returns infinitesimal numbers instead of zero when an eigenvalue
    # is 0. Below, these eigenvalues are rounded to 0. The non-zero
    # values are already integer, so rounding does not affect them.
    clusters = np.asarray(np.round(clusters), dtype=int)
    cluster_sizes, cluster_counts = np.unique(clusters, return_counts=True)
    # Above, the eigenvalue with value 0 can be interpreted in this way:
    # The system has clusters with size equal to or larger than 1
    # The frequency of eigenvalue 0 is the total number of clusters (or
    # dependant columns in the contact matrix) with size equal to
    # or larger than 2.
    np.put(cluster_list, cluster_sizes, cluster_counts)
    return cluster_list
