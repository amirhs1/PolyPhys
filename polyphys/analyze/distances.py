"""
Distance Computations --- :mod:`polyphys.analyze.distances`
===========================================================

The :mod:`polyphys.analyze.distances` module provides functions for calculating
pairwise distances from molecular coordinates.
Functions
=========
.. autofunction:: self_dist_array
.. autofunction:: self_dist_array_opt

Notes
=====
- Distance calculations can account for periodic boundary conditions, ensuring
  correct measurements in periodic simulation boxes.
- Optimized versions of functions are available for improved performance in
  certain dimensions and scenarios.

References
==========
- Albanie, S. (2019). Euclidean Distance Matrix Trick. In Advances in Computer
  Vision. Lecture Notes in Computer Science, vol 11131. Springer, Cham.
  https://doi.org/10.1007/978-3-030-01446-4_2
"""
from typing import Optional, Dict
import numpy as np
from .measurer import apply_pbc_orthogonal
from ..manage.types import AxisT


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
