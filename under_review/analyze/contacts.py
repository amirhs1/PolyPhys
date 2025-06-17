"""
Contact analysis --- :mod:`polyphys.analyze.contacts`
=====================================================


This module contains classes to analyze contacts between atoms over a
trajectory.


Functions
---------

.. autofunction:: find_direct_contacts
.. autofunction:: generate_contact_matrix
.. autofunction:: random_full_contact_matrix
"""
from typing import Optional
from itertools import combinations
import numpy as np


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
        The filename or filepath to which the contact matrix is saved. If None,
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
