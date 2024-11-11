"""Miscellaneous tools used by different modules.
"""
import re
import numpy as np
import gzip
from gzip import GzipFile
import os
from typing import Generator, Optional, List, IO, TextIO, Any, Union, Literal, Tuple
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


def sort_filenames(
    filenames: List[str], fmts: List[str]
) -> List[Tuple[str, ...]]:
    """
    Groups and alphanumerically sorts filenames by specified formats.

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


def create_fullname(name: str, group: str, prop: str) -> str:
    """
    Creates a structured filename based on a base name, particle group,
    and property.

    Parameters
    ----------
    name : str
        The base name for the output file.
    group : str
        Type of the particle group (e.g., 'bug' or 'nucleoid').
    prop : str
        The physical property name for the data (e.g., 'density').

    Returns
    -------
    str
        A string that combines `name`, `group`, and `prop` in a
        hyphen-separated format, suitable for use as a filename.

    Examples
    --------
    >>> create_fullname("sample", "bug", "density")
    'sample-bug-density'
    """
    return "-".join([name, group, prop])


def save_parent(
    filename: str,
    data: Union[np.ndarray, pd.DataFrame, Dict[str, np.ndarray]],
    save_to: str,
) -> None:
    """
    Saves the `data` to a specified file format, allowing structured file
    naming.

    Parameters
    ----------
    filename : str
        The output file name.
    data : Union[np.ndarray, pd.DataFrame, Dict[str, np.ndarray]]
        Data to be saved. Accepted types:

        - `np.ndarray`: Saved as a .npy file.
        - `pd.DataFrame`: Saved as a .csv file.
        - `dict` of `np.ndarray`: Each entry saved as a separate .npy file
        with suffix.

    save_to : str
        Path to the directory where the file will be saved.

    Raises
    ------
    TypeError
        If `data` type does not match any of the supported formats.

    Examples
    --------
    Save a numpy array to a `.npy` file in the specified directory:

    >>> save_parent("output", np.array([1, 2, 3]), "density", "all", "/data/")

    Save a DataFrame to a `.csv` file:

    >>> import pandas as pd
    >>> df = pd.DataFrame({"A": [1, 2, 3]})
    >>> save_parent("output", df, "density", "all", "/data/")
    """
    file_path = os.path.join(save_to, filename)

    # Save data based on type
    if isinstance(data, pd.DataFrame):
        data.to_csv(f"{file_path}.csv", index=False)
    elif isinstance(data, np.ndarray):
        np.save(f"{file_path}.npy", data)
    elif isinstance(data, dict):
        for prop_key, prop_data in data.items():
            _, prop_measure = prop_key.split("-")
            np.save(f"{file_path}-{prop_measure}.npy", prop_data)
    else:
        raise TypeError(
            f"Unsupported data type {type(data).__name__}."
            + "Expected pd.DataFrame, np.ndarray, or dict of np.ndarray."
        )


def make_database(
    old_database: str,
    phase: Literal[
        "simAll", "simCont", "log", "trj", "probe", "analysis", "viz", "galaxy"
    ],
    stage: Literal["segment", "wholeSim", "ens", "ensAvg", "space", "galaxy"],
    group: Literal["bug", "nucleoid", "all"],
) -> str:
    """
    Create a new directory path based on the provided `old_database` path and
    specified parameters (`phase`, `group`, and `stage`). If the directory does
    not already exist, it will be created.

    The `old_database` is expected to follow a structured naming convention:
    `prefix-old_phase-old_group-old_stage`, where each part represents an
    aspect of the directory's role. This base structure helps generate the
    new path.

    The newly constructed directory name follows:
        `prefix-phase-group-stage`

    Parameters
    ----------
    old_database : str
        Path to the reference directory whose structure will be used as a base.
    phase : {"simAll", "simCont", "log", "trj", "probe", "analysis", "viz",
            "galaxy"}
        The new phase name for the directory, specifying its role in
        the workflow.
    stage : {"segment", "wholeSim", "ens", "ensAvg", "space", "galaxy"}
        The stage of processing or data type represented in the new directory.
    group : {"bug", "nucleoid", "all"}
        Type of particle group related to the data in the directory.

    Returns
    -------
    str
        The path to the newly created directory or the existing path if it
        already exists.

    Raises
    ------
    ValueError
        If any parameter is not in its list of accepted values.

    Examples
    --------
    Examples of directory transformation:
        - Input: `old_database = root/parent1/.../parentN/old_phase/old_dir`
        - Result: `new_database = root/parent1/.../parentN/phase/new_dir`

    Construct a new directory path based on an existing database path:

    >>> make_database('/root/data/old_analysis/simulations', 'analysis',
                     'wholeSim', 'bug')
    '/root/data/analysis/prefix-bug-wholeSim/'
    """
    invalid_keyword(
        phase,
        ["simAll", "simCont", "log", "trj", "probe", "analysis", "viz",
         "galaxy"],
    )
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(
        stage, ["segment", "wholeSim", "ens", "ensAvg", "space", "galaxy"]
    )

    old_path = pathlib.Path(old_database)

    # Common prefix derived from old_database to create new directory name
    prefix = old_path.parts[-1].split("*")[0].split("-")[0]
    new_directory = "-".join([part for part in [prefix, group, stage] if part])

    # Construct the new database path in the same parent directory level as
    # the old one
    new_database_parts = list(old_path.parts[:-2]) + [phase, new_directory]
    new_database = pathlib.Path(*new_database_parts)

    # Avoid double slashes at the start of the path
    if str(new_database).startswith("//"):
        new_database = pathlib.Path(str(new_database)[1:])

    # Create the new directory if it doesn't already exist
    try:
        new_database.mkdir(parents=True, exist_ok=False)
    except FileExistsError as error:
        print(error)
        print("Files are saved/overwritten in an existing directory.")

    return str(new_database) + "/"


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
