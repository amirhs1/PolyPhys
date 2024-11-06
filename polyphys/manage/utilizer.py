"""Miscellaneous tools used by different modules.
"""
import numpy as np
import gzip
from gzip import GzipFile
import os
from typing import Generator, Optional, List, IO, TextIO, Any, Union
from contextlib import contextmanager

InputT = Union[GzipFile, TextIO, IO[Any]]


def read_camel_case(word: str) -> List[Union[str, Tuple[str]]]:
    """
    Splits a camelCase or CamelCase string into its component words.

    Parameters
    ----------
    word : str
        The camelCase or CamelCase string to be split.

    Returns
    -------
    str
        The split string with spaces inserted between words.

    Examples
    --------
    >>> read_camel_case("camelCase")
    ['camel', 'Case']
    >>> read_camel_case("CamelCaseString")
    ['Camel', 'Case', 'String']
    """
    return re.findall(r"[A-Z]?[a-z]+|[A-Z]+(?=[A-Z]|$)", word)


def to_float_if_possible(value: str) -> Union[float, str]:
    """
    Attempts to convert a string to a float. If conversion fails,
    returns the original string.

    Parameters
    ----------
    value : str
        The string to attempt to convert to a float.

    Returns
    -------
    Union[float, str]
        The converted float if the input can be converted; otherwise,
        the original string.

    Examples
    --------
    >>> to_float_if_possible("3.14")
    3.14
    >>> to_float_if_possible("not_a_float")
    'not_a_float'
    >>> to_float_if_possible("42")
    42.0
    """
    try:
        return float(value)
    except ValueError:
        return value


def split_alphanumeric(alphanumeric: str) -> List[Union[int, str, float]]:
    """
    Splits an alphanumeric string into a list of strings, integers, and floats.

    This function identifies contiguous sections of digits, letters, and
    decimal numbers, returning them as separate elements in a list, making the
    string suitable for alphanumeric sorting.

    Parameters
    ----------
    alphanumeric : str
        An alphanumeric string to be split.

    Returns
    -------
    List[Union[int, str, float]]
        A list of components including integers, floats, and strings.

    Examples
    --------
    >>> split_alphanumeric("file20.5name10")
    ['file', 20.5, 'name', 10]
    """
    number_pattern = re.compile(r"(\d+\.*\d*)")
    parts = number_pattern.split(alphanumeric)
    parts = [
        int(part) if part.isdigit() else to_float_if_possible(part)
        for part in parts
        if part
    ]
    return parts


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
