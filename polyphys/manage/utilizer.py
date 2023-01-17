"""Miscellaneous tools used by different modules.
"""
import numpy as np
import gzip
from gzip import GzipFile
import os
from typing import Generator, Optional, List, IO, TextIO, Any, Union
from contextlib import contextmanager

InputT = Union[GzipFile, TextIO, IO[Any]]


def openany(
    filepath: str,
    mode: str = 'r'
) -> InputT:
    """
    Open a regular or gzipped file.

    Parameters
    ----------
    filepath : str
        Path to the file.
    mode : str, optional
        The mode by the file is opened, by default 'r'

    Yields
    ------
    Generator[IO, None, None]
        A file object.
    """
    _, ext = os.path.splitext(filepath)
    if ext == ".gz":
        file = gzip.open(filepath, mode=mode)
    else:
        file = open(filepath, mode=mode)
    return file


@contextmanager
def openany_context(
    filepath: str,
    mode: str = 'r'
) -> Generator[InputT, None, None]:
    """
    Open a regular or gzipped file.

    Parameters
    ----------
    filepath : str
        Path to the file.
    mode : str, optional
        The mode by the file is opened, by default 'r'

    Yields
    ------
    Generator:
        A file object.
    """
    _, ext = os.path.splitext(filepath)
    if ext == ".gz":
        file = gzip.open(filepath, mode=mode)
    else:
        file = open(filepath, mode=mode)
    try:
        yield file
    except Exception as exp:
        print(exp)
    file.close()


def round_down_first_non_zero(x: float) -> float:
    """rounds down a number to its first non-zero digit.

    Parameters
    ----------
    x : float
        number which is rounded.

    Returns
    -------
    float:
        a float number to which x is rounded down to its first non-zero digit.
    """
    if x == 0:
        return x
    else:
        exponent = np.floor(np.log10(abs(x)))
        non_zero = 10 ** exponent
        return round(np.floor(x/non_zero)*non_zero, int(abs(exponent)))


def round_up_nearest(
    dividend: float,
    divider: float,
    round_to: int
) -> float:
    """rounds up `dividend` by `divider` up to `round_to` significant digits.

    Parameters
    ----------
    dividend: float
        The number should be rounded up.
    divider: float
        The number used as the divider.
    round_to: int
        The number of the significant digits.

    Return
    ------
    float:
        The rounded number which is divisible by the divisor.
    """
    return round(round(dividend / divider) * divider, round_to)


def invalid_keyword(
    keyword: str,
    valid_keywords: List[str],
    message: Optional[str] = None
) -> None:
    """
    Raises an error if `keyword` is not in `valid_keywords`.

    Parameters
    ----------
    keyword: str
        Name of the `keyword`
    valid_keywords: array of str
        Array of valid keywords
    message: str
        Message to be printed.
    """
    if message is None:
        message = " is an invalid option. Please select one of " + \
            f"{valid_keywords} options."
    if keyword not in valid_keywords:
        raise ValueError(f"'{keyword}'" + message)
