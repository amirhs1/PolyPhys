"""
Contains tools fro perfroming statistical measurements such as avaerging,
standard deviation, or the like on a given time series, correlation data,
or histogram.
"""

import numpy as np


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
