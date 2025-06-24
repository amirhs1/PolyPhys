"""\
===============================================================
Miscellaneous utility codes --- :mod:`polyphys.manage.utils`
===============================================================

The :mod:`~polyphys.manage.utilizer` module provides shared utility functions
used across various parts of the PolyPhys package. These functions handle
string processing, file handling, numeric rounding, and error checking.


Functions
=========

.. autofunction:: read_camel_case
.. autofunction:: to_float_if_possible
.. autofunction:: split_alphanumeric
.. autofunction:: sort_filenames
.. autofunction:: openany
.. autofunction:: openany_context
.. autofunction:: round_down_first_non_zero
.. autofunction:: round_up_nearest
.. autofunction:: invalid_keyword

Dependencies
============
- Python standard library: `re`, `gzip`, `contextlib`, `typing`
- Third-party: `numpy`

See also: :mod:`~polyphys.manage.types` for the definition of
:data:`InputType`.
"""
import re
import gzip
from typing import Generator, Optional, List, Union, Tuple, Sequence
from contextlib import contextmanager
import numpy as np
from .types import InputType


def read_camel_case(string: str) -> List[Union[str, Tuple[str]]]:
    """
    Split a *camelCase* or *CamelCase* `string` into its component substrings.

    Parameters
    ----------
    string : str
        The *camelCase* or *CamelCase* `string` to be split.

    Returns
    -------
    List[str]
        The split `string` with substrings as elements.

    Examples
    --------
    >>> read_camel_case("camelCase")
    ['camel', 'Case']
    >>> read_camel_case("CamelCaseString")
    ['Camel', 'Case', 'String']
    """
    return re.findall(r"[A-Z]?[a-z]+|[A-Z]+(?=[A-Z]|$)", string)


def to_float_if_possible(string: str) -> Union[float, str]:
    """
    Attempt to convert a `string` to a float. If conversion fails,
    returns the original `string`.

    Parameters
    ----------
    string : str
        The `string` to attempt to convert to a float.

    Returns
    -------
    Union[float, str]
        The converted float if the input can be converted; otherwise,
        the original `string`.

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
        return float(string)
    except ValueError:
        return string


def split_alphanumeric(alphanumeric: str) -> List[Union[int, str, float]]:
    """
    Split an *alphanumeric* string into a list of strings, integers, and
    floats.

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
    filenames: List[str],
    formats: List[Union[str, Tuple[str, ...]]]
) -> List[Tuple[str, ...]]:
    """
    Group and alphanumerically sort `filenames` by specified formats.

    This function groups `filenames` by the extensions in `formats`, sorting
    each group alphanumerically. It then returns a list of tuples where each
    tuple contains filenames with matching base names across the specified
    formats.

    Parameters
    ----------
    filenames : List[str]
        A list of filenames to sort and group.
    formats : List[Union[str, Tuple[str, ...]]]
        A list specifying the file formats. Each item in `formats` can be a
        string representing a file extension (e.g., 'data'), or a tuple of
        alternative extensions considered equivalent (e.g.,
        ('trj', 'lammpstrj')).

    Returns
    -------
    List[Tuple[str, ...]]
        A list of tuples where each tuple contains filenames grouped and
        sorted by the specified formats.

    Notes
    -----
    If any of the specified formats match no files, the returned list may be
    empty. This is because Python's built-in `zip()` truncates to the shortest
    input sequence. It is the caller's responsibility to ensure that each
    format group matches at least one file.

    Consider validating inputs beforehand or adding error handling if such
    cases should be treated as exceptional.

    Examples
    --------
    >>> sort_filenames(['file1.data', 'file2.trj', 'file1.trj', 'file2.data'],
    ...                ['data', ('lammpstrj', 'trj')])
    [('file1.data', 'file1.trj'), ('file2.data', 'file2.trj')]

    >>> sort_filenames(['file1.data', 'file2.data'], ['data', ('trj',)])
    []
    """
    grouped_filenames = []
    for exts in formats:
        grouped_filenames.append([f for f in filenames if f.endswith(exts)])

    print(grouped_filenames)
    for idx, filenames_group in enumerate(grouped_filenames):
        grouped_filenames[idx] = sorted(
            filenames_group, key=split_alphanumeric
        )
    return list(zip(*grouped_filenames))


def openany(filepath: str, mode: str = 'rt') -> InputType:
    """
    Open a regular or gzipped file.

    Parameters
    ----------
    filepath : str
        Path to the file.
    mode : str, optional
        The mode by the file is opened, by default 'r'.

    Returns
    -------
    InputType
        A file-like object (e.g., regular file or GzipFile).
    """
    open_func = gzip.open if filepath.endswith(".gz") else open
    return open_func(filepath, mode=mode)


@contextmanager
def openany_context(
    filepath: str,
    mode: str = 'rt'
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
        A file-like object (e.g., regular file or GzipFile) opened in the
        specified mode.
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


def round_down_first_non_zero(num: float) -> float:
    """
    Round down number `num` to its first non-zero digit.

    Parameters
    ----------
    num : float
        Number which is rounded.

    Returns
    -------
    float:
        The number `num` rounded down to its first non-zero digit.
    """
    if num == 0:
        return num
    exponent = np.floor(np.log10(abs(num)))
    non_zero = 10 ** exponent
    return round(np.floor(num/non_zero)*non_zero, int(abs(exponent)))


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
        The number to round up.
    divider: float
        The number used as the divider.
    round_to: int
        The number of the significant digits.

    Returns
    -------
    float:
        The smallest number greater than or equal to `dividend` that is
        divisible by `divider`, rounded to `round_to` significant digits.
    """
    return round(round(dividend / divider) * divider, round_to)


def invalid_keyword(
    keyword: str,
    valid_keywords: Sequence[Optional[str]],
    message: Optional[str] = None
) -> None:
    """
    Raise an error if `keyword` is not in `valid_keywords`.

    Parameters
    ----------
    keyword: str
        Name of the `keyword`
    valid_keywords: List[Union[str,None]]
        Array of valid keywords
    message: str
        Message to be printed.

    Raises
    ------
    ValueError
        If `keyword` is not in `valid_keywords`.
    """
    if message is None:
        message = " is an invalid option. Please select one of: " + \
            f"{', '.join(map(str, valid_keywords))}."
    if keyword not in valid_keywords:
        raise ValueError(f"'{keyword}'" + message)


def number_density_cube(
    n_atom: float, d_atom: float, l_cube: float, pbc: bool = False
) -> float:
    """
    Compute the bulk number density of a species in a cubic box.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cube : float
        Length of one side of the cubic box.
    pbc : bool
        Periodic boundary conditions along all box sides. If `True`,
        :math:`v_{avail} = l_{cube}^3`; otherwise,
        :math:`v_{avail} = (l_{cube} - d_{atom})^3`. Defaults to `False`.

    Returns
    -------
    float
        Bulk number density of the species in the cubic box.

    Notes
    -----
    The bulk number density is calculated as :math:`n_{atom} / v_{avail}`,
    where `v_{avail}` is the available volume to the center of geometry of a
    particle based on the presence of periodic boundary conditions.

    Examples
    --------
    >>> number_density_cube(1000, 1.0, 10.0)
    1.3717421124828531
    >>> number_density_cube(1000, 1.0, 10.0, pbc=True)
    1.0
    """
    v_avail = (l_cube - int(not pbc) * d_atom) ** 3
    return n_atom / v_avail


def volume_fraction_cube(
    n_atom: float, d_atom: float, l_cube: float, pbc: bool = False
) -> float:
    """
    Compute the volume fraction of a species in a cubic box.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cube : float
        Length of one side of the cubic box.
    pbc : bool
        Periodic boundary conditions along all box sides. If `True`,
        :math:`v_{avail} = l_{cube}^3`; otherwise,
        :math:`v_{avail} = (l_{cube} - d_{atom})^3`. Defaults to `False`.

    Returns
    -------
    float
        Volume fraction of the species in the cubic box.

    Raises
    ------
    ValueError
        If the computed volume fraction exceeds 1.0, which is physically
        invalid for spherical particles.

    Notes
    -----
    The available volume accounts for the excluded-volume boundary effect if
    `pbc=False`, meaning particles cannot cross the box boundary. If
    `pbc=True`, particles are assumed to wrap around, and the full box volume
    is used.

    Examples
    --------
    >>> volume_fraction_cube(1000, 1.0, 10.0)
    0.7189611722461486

    >>> volume_fraction_cube(1000, 1.0, 10.0, pbc=True)
    0.5235987755982988
    """
    rho = number_density_cube(n_atom, d_atom, l_cube, pbc)
    phi = rho * np.pi * d_atom**3 / 6
    if phi > 1.0:
        raise ValueError(
            "Volume fraction exceeds 1.0, which is physically invalid."
        )
    return phi


def number_density_cylinder(
    n_atom: float, d_atom: float, l_cyl: float, d_cyl: float, pbc: bool = False
) -> float:
    """
    Compute the bulk number density of a species in a cylindrical confinement.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of the particle of the species.
    l_cyl : float
        Length of the cylindrical confinement.
    d_cyl : float
        Diameter of the cylindrical confinement.
    pbc : bool
        Periodic boundary conditions along the longitudinal axis. If `True`,
        :math:`v_{avail} = \\pi * l_{cyl} * (d_{cyl} - d_{atom})^2 / 4`;
        otherwise, :math:`v_{avail} = \\pi * (l_{cyl} - d_{atom}) * (d_{cyl}
        - d_{atom})^2 / 4`. Defaults to `False`.

    Returns
    -------
    float
        Bulk number density of the species in the cylindrical confinement.

    Notes
    -----
    The bulk number density is calculated as :math:`n_{atom} / v_{avail}`,
    where `v_{avail}` is the available volume to the center of geometry of a
    particle based on the presence of periodic boundary conditions.

    Examples
    --------
    >>> number_density_cylinder(100, 1.0, 10.0, 5.0)
    0.8841941282883075
    >>> number_density_cylinder(100, 1.0, 10.0, 5.0, pbc=True)
    0.7957747154594768
    """
    v_avail = (
        np.pi * (l_cyl - int(not pbc) * d_atom) * (d_cyl - d_atom) ** 2 / 4
    )
    return n_atom / v_avail


def volume_fraction_cylinder(
    n_atom: float, d_atom: float, l_cyl: float, d_cyl: float, pbc: bool = False
) -> float:
    """
    Compute the volume fraction of spherical particles in a cylindrical
    confinement.

    Parameters
    ----------
    n_atom : float
        Number of particles.
    d_atom : float
        Diameter of a single particle.
    l_cyl : float
        Length of the cylindrical confinement (along the longitudinal axis).
    d_cyl : float
        Diameter of the cylindrical confinement (transverse direction).
    pbc : bool, default False
        Whether periodic boundary conditions are applied along the cylinder
        axis. If `True`, the available volume is:

        .. math::

            V = \\frac{\\pi}{4} l_{cyl} (d_{cyl} - d_{atom})^2

        If `False`, the available volume is:

        .. math::

            V = \\frac{\\pi}{4} (l_{cyl} - d_{atom}) (d_{cyl} - d_{atom})^2

    Returns
    -------
    float
        Volume fraction of the particles, defined as:

        .. math::

            \\phi = \\rho \\cdot \\frac{\\pi d^3}{6}

        where :math:`\\rho` is the number density of particles in the
        available cylindrical volume.

    Raises
    ------
    ValueError
        If the computed volume fraction exceeds 1.0, which is physically
        invalid.

    Notes
    -----
    The available volume excludes regions where the particle centers cannot
    access due to finite size. When `pbc=False`, particles are excluded from
    the cylinder ends; when `pbc=True`, they are not.

    Examples
    --------
    >>> volume_fraction_cylinder(100, 1.0, 10.0, 5.0)
    0.46296296296296297

    >>> volume_fraction_cylinder(100, 1.0, 10.0, 5.0, pbc=True)
    0.4166666666666667
    """
    rho = number_density_cylinder(n_atom, d_atom, l_cyl, d_cyl, pbc)
    phi = rho * np.pi * d_atom**3 / 6
    if phi > 1.0:
        raise ValueError(
            "Volume fraction exceeds 1.0, which is physically invalid."
        )
    return phi
