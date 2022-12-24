"""Miscellaneous tools used by different modules.
"""
import numpy as np
import gzip
import os
from typing import List, Generator, IO, Any
from contextlib import contextmanager


def openany(filepath: str, mode: str = 'r') -> IO[Any]:
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
def openany_context(filepath: str, mode: str = 'r'
                    ) -> Generator[IO[Any], None, None]:
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
    try:
        yield file
    except Exception as exp:
        print(exp)
    file.close()


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


def round_up_nearest(dividend: float, diviser: float, round_to: int) -> float:
    """rounds up `dividend` by `diviser` up to `round_to` significatn digits.

    Parameters
    ----------
    dividend: float
        The number should be rounded up.
    diviser: float
        The number used as the diviser.
    round_to: int
        The nuumber of the significatn digits.

    Return
    ------
    float:
        The rounded number which is divisable by the divisor.
    """
    # pylance: disable-next=reportGeneralTypeIssues
    return np.around(np.around(dividend / diviser) * diviser, round_to)


def invalid_keyword(
    keyword: str,
    valid_keywords: List[str],
    message: str = ''
) -> None:
    """
    Raises an error if `keyowrd` is not in `valid_keywords`.

    Parameters
    ----------
    keyword: str
        Name of the `keyword`
    valid_keywords: array of str
        Array of valid keywords
    """
    if message is None:
        message = " is an invalid option. Please select one of " + \
            f"{valid_keywords} options."
    if keyword not in valid_keywords:
        raise ValueError(f"'{keyword}'" + message)
