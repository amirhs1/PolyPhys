"""Miscellaneous tools used by different modules.
"""
import re
import gzip
from typing import Generator, Optional, List, Union, Tuple
from contextlib import contextmanager
import numpy as np
from .typer import InputType


def read_camel_case(word: str) -> List[Union[str, Tuple[str]]]:
    """
    Split a camelCase or CamelCase string into its component words.

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
    Attempt to convert a string to a float. If conversion fails,
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
    Split an alphanumeric string into a list of strings, integers, and floats.

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
    return [int(part) if part.isdigit() else
            to_float_if_possible(part) for part in parts if part]


def sort_filenames(
    filenames: List[str], fmts: List[str]
) -> List[Tuple[str, ...]]:
    """
    Group and alphanumerically sort filenames by specified formats.

    This function groups `filenames` by the extensions in `formats`, sorting
    each group alphanumerically. It then returns a list of tuples where each
    tuple contains filenames with matching base names across the specified
    formats.

    Parameters
    ----------
    filenames : List[str]
        A list of filenames to sort and group.
    fmts : List[str]
        A list specifying the file formats. Each item can be a single extension
        (e.g., 'data') or a tuple of alternative extensions (e.g., ('trj',
        'lammpstrj')).

    Returns
    -------
    List[Tuple[str, ...]]
        A list of tuples where each tuple contains filenames grouped and
        sorted by the specified formats.

    Examples
    --------
    >>> sort_filenames(['file1.data', 'file2.trj', 'file1.trj', 'file2.data'],
                       ['data', ('lammpstrj', 'trj')])
    [('file1.data', 'file1.trj'), ('file2.data', 'file2.trj')]
    """
    grouped_filenames = []
    for exts in fmts:
        grouped_filenames.append([f for f in filenames if f.endswith(exts)])
    for idx, filenames_group in enumerate(grouped_filenames):
        grouped_filenames[idx] = sorted(
            filenames_group, key=split_alphanumeric
        )
    return list(zip(*grouped_filenames))


def openany(filepath: str, mode: str = 'r') -> InputType:
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
    open_func = gzip.open if filepath.endswith(".gz") else open
    return open_func(filepath, mode=mode)


@contextmanager
def openany_context(
    filepath: str,
    mode: str = 'r'
) -> Generator[InputType, None, None]:
    """
    Open a regular or gzipped file, providing a file-like object.

    Parameters
    ----------
    filepath : str
        Path to the file.
    mode : str, optional
        The mode by the file is opened, by default 'r'

    Yields
    ------
    Generator[InputT, None, None]
        A file object.
    """
    open_func = gzip.open if filepath.endswith(".gz") else open
    file = open_func(filepath, mode=mode)
    try:
        yield file
    except Exception as exp:
        print(f"Error opening {filepath}: {exp}")
        raise
    finally:
        file.close()


def round_down_first_non_zero(x: float) -> float:
    """
    Round down a number to its first non-zero digit.

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
    exponent = np.floor(np.log10(abs(x)))
    non_zero = 10 ** exponent
    return round(np.floor(x/non_zero)*non_zero, int(abs(exponent)))


def round_up_nearest(
    dividend: float,
    divider: float,
    round_to: int
) -> float:
    """
    Round up `dividend` by `divider` up to `round_to` significant digits.

    Parameters
    ----------
    dividend: float
        The number should be rounded up.
    divider: float
        The number used as the divider.
    round_to: int
        The number of the significant digits.

    Returns
    -------
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
    Raise an error if `keyword` is not in `valid_keywords`.

    Parameters
    ----------
    keyword: str
        Name of the `keyword`
    valid_keywords: array of str
        Array of valid keywords
    message: str
        Message to be printed.

    Raises
    ------
    ValueError
        If `keyword` is not in `valid_keywords`.
    """
    if message is None:
        message = " is an invalid option. Please select one of " + \
            f"{valid_keywords} options."
    if keyword not in valid_keywords:
        raise ValueError(f"'{keyword}'" + message)
