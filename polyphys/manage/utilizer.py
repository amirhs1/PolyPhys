"""Miscellaneous tools used by different modules.
"""
import numpy as np


def round_down_first_non_zero(x: float) -> float:
    """rounds down a number to its first non-zero digit.

    Parameters:
    x (float or int): number which is rounded.

    Returns:
    a float or integer number to which to which x is rounded down to its
        first non-zero digit.
    """
    if x == 0:
        return x
    else:
        exponent = np.floor(np.log10(abs(x)))
        non_zero = 10 ** exponent
        return round(np.floor(x/non_zero)*non_zero, int(abs(exponent)))


def round_up_nearest(dividend: float, diviser: float) -> float:
    """rounds up the floating-point number dividend by the diviser.

    Parameters
    ----------
    dividend: float
        The number should be rounded up.
    diviser: float
        The number used as the diviser.

    Return
    ------
    float:
        The rounded number which is divisable by the divisor.
    """
    return np.round(dividend / diviser) * diviser
