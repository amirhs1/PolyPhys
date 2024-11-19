"""
Helper functions --- :mod:`polyphys.analyze.helper`
===================================================


Helper functions for mathematical operations.

Functions
---------

.. autofunction:: is_symmetric
.. autofunction:: is_positive_semi_definite
"""
import numpy as np


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
