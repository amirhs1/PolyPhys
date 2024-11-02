"""This oduel organize fa hierarchy of files.
"""

import warnings
from typing import List, Dict, Tuple, Optional, Union, Callable
import os
import pathlib
from glob import glob
import re
import numpy as np
import pandas as pd

from ..manage.typer import ParserT, TransFociT
from ..analyze.clusters import whole_dist_mat_foci
from .utilizer import round_up_nearest, invalid_keyword


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
    filenames: List[str],
    fmts: List[str]
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


def database_path(
    old_database: str,
    phase: str,
    stage: Optional[str] = None,
    group: Optional[str] = None,
) -> str:
    """
    Constructs and creates a directory path based on the specified
    `old_database` path and new parameters `phase`, `group`, and `stage`. The
    resulting directory is created at the same level as `old_database`.

    The `old_database` is expected to follow a structured naming pattern:
    `prefix-old_phase-old_group-old_stage`, which helps infer the base
    structure.

    The new directory will be named as:
        `prefix-phase-group-stage`

    Given a directory structure:
        old_database = root/parent1/.../parentN/old_phase/old_directory
    This function would create:
        new_database = root/parent1/.../parentN/phase/new_directory

    Parameters
    ----------
    old_database : str
        The path to the existing directory, which serves as a reference point.
    phase : str
        A new phase name for the directory.
    stage : Optional[str], default=None
        The stage of the new directory.
    group : Optional[str], default=None
        The particle group type.

    Returns
    -------
    str
        The path to the newly created directory, or an existing directory if it
        already exists.

    Raises
    ------
    ValueError
        If any parameter is not in its list of accepted values.

    Examples
    --------
    Construct a new directory path based on an existing database path:

    >>> database_path('/root/data/old_analysis/simulations', 'analysis',
                     'wholeSim', 'bug')
    '/root/data/analysis/prefix-bug-wholeSim/'
    """
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


def whole_from_segment(
    property_: str,
    segments: List[Tuple[str]],
    parser: ParserT,
    geometry: str,
    group: str,
    topology: str,
    relation: str,
    save_to: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """
    Generates `whole` array for `property_` of the particle `group` in the
    `geometry` of interest from its `segments`.

    Parameters
    ----------
    property_: str
        The physical property.
    segments: list of tuples
        List of tuples where each tuple at least has one member (the path to
        a csv file for the `property_`).
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    topology: str
        Topology of the polymer.
    relation : {'histogram', 'tseries', 'bin_edges'}
        Relation between segments and wholes:

        'histogram'
            Child is an N-dimensional histogram-like 'segment' file, so the
            children of a parent should be sum along "axis=0" in numpy lingo.
            In 2-dimensional space, a histogram has a (x_nbins, y_nbins) shape.

        'tseries'
            Child is a time-series-like 'segment' file, so the children of
            a parent should be concatenated vertically (along "axis=0" in
            pandas lingo).

        'bin_edge'
            Child is a bin edges 'segment' file. All the siblings' 'segments'
            are the same, so one parent's 'bin_edges' is created by rewriting
            similar segments on each other.

    save_to : str, defualt None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    siblings: dict of np.ndarray
        Dict of siblings where keys are 'whole' names (str) and values
        are whole data (arrays).

    Notes
    -----
    Please see the 'organizer' documentation for definitions
    of 'geometry', and 'group', and the definitions of their keywords.
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(relation, ["histogram", "tseries", "bin_edge"])
    mapping_func = {
        "histogram": lambda whole: (whole[0], np.sum(whole[1], axis=0)),
        "tseries": lambda whole: (whole[0], np.concatenate(whole[1])),
        # 'segments' of a 'whole' have the same bin_edge, so np.unique used
        # to pick one.
        "bin_edge": lambda whole: (whole[0], np.unique(np.array(whole[1]))),
    }
    wholes = {}
    for segment in segments:
        segment_info = parser(segment[0], "segment", geometry, group, topology)
        whole_name = getattr(segment_info, "whole")
        child_arr = np.load(segment[0])
        # if not bool(wholes):  # is ens_names empty or not?
        #    wholes[whole_name] = [child_arr]
        # elif whole_name not in wholes.keys():
        #    wholes[whole_name] = [child_arr]
        # a better approach 2023-08-15
        if whole_name not in wholes:
            wholes[whole_name] = [child_arr]
        else:
            wholes[whole_name].append(child_arr)
    wholes = dict(map(mapping_func[relation], wholes.items()))
    if save_to is not None:
        _ = dict(
            map(
                lambda whole: (
                    whole[0],
                    save_parent(whole[0], whole[1], property_, save_to, group)
                ),
                wholes.items(),
            )
        )
    return wholes


def whole_from_file(
    whole_paths: List[Tuple[str]],
    parser: ParserT,
    geometry: str,
    group: str,
    topology: str,
) -> Dict[str, np.ndarray]:
    """Loads `whole` numpy arrays for a given physical property of the
    particle `group` in the `geometry` of interest from their paths
    `whole_paths`.

    Parameters
    ----------
    whole_paths : list of tuples
        List of tuples where each tuple at least has one member (the path to
        a csv file for the `property_`).
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    topology: str
        Topology of the polymer.

    Return
    ------
    siblings: dict of np.ndarray
        Dict of siblings where keys are 'whole' names (str) and values
        are whole data (arrays).

    Notes
    -----
    Please see the 'organizer' documentation for definitions
    of 'geometry', and 'group', and the definitions of their keywords.
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    wholes = {}
    for whole_path in whole_paths:
        whole_info = parser(whole_path[0], "whole", geometry, group, topology)
        whole_name = getattr(whole_info, "whole")
        wholes[whole_name] = np.load(whole_path[0])
    return wholes


def whole_from_dist_mat_t(
    whole_paths: List[Tuple[str]],
    parser: TransFociT,
    geometry: str,
    group: str,
    topology: str,
) -> Tuple[Dict[str, pd.DataFrame], ...]:
    """
    Loads `whole` 2D numpy arrays for a given physical property of the
    particle `group` in the `geometry` of interest from their paths
    `whole_paths`.

    Parameters
    ----------
    whole_paths : list of tuples
        List of tuples where each tuple at least has one member (the path to
        a csv file for the `property_`).
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    topology: str
        Topology of the polymer.

    Return
    ------
    wholes_freqs: dict of np.ndarray
    wholes_rdfs: dict of np.ndarray
    wholes_tseries: dict of np.ndarray
        In each of these three dicts, keys are 'whole' names (str) and values
        are whole data (arrays).

    Notes
    -----
    Please see the 'organizer' documentation for definitions
    of 'geometry', and 'group', and the definitions of their keywords.
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    wholes_freqs = {}
    wholes_rdfs = {}
    wholes_tseries = {}
    for whole_path in whole_paths:
        whole_info: TransFociT = parser(
            whole_path[0], "whole", geometry, group, topology
        )
        whole_freqs, whole_rdfs, whole_tseries = whole_dist_mat_foci(
            whole_path[0], whole_info
        )
        whole_name = getattr(whole_info, "whole")
        wholes_freqs[whole_name] = whole_freqs
        wholes_rdfs[whole_name] = whole_rdfs
        wholes_tseries[whole_name] = whole_tseries
    return wholes_freqs, wholes_rdfs, wholes_tseries


def ens_from_bin_edge(ens) -> Tuple[str, np.ndarray]:
    """
    Convert an "ensemble" bin edge from several "whole" bin edges where all
    the "whole" bin edges are similar but named differently because they
    belong to different wholes.

    Parameters
    ----------
    ens: np.ndarray
        All the "whole" bin edges in a given `ens`.

    Returns
    -------
    tuple of str and np.ndarray
        The "ensemble" name  and its bin edge values.
    """
    return ens[0], np.unique(list(ens[1].values()))


def ens_from_vec(
    ens: Tuple[str, Dict[str, np.ndarray]]
) -> Tuple[str, pd.DataFrame]:
    """
    Creates an "ensemble" dataframe from a dictionary of wholes where
    each "whole" is a numpy vector or 1D array.

    Parameters
    ----------
    ens: tuple of np.ndarray
        A tuple in which the first element is an ensemble name and the second
        one is a dictionary in which  keys are whole names and values are
        whole-type data.

    Return
    ------
    A tuple of ensemble name and its associated dataframe.
    """
    return ens[0], pd.DataFrame.from_dict(ens[1], orient="columns")


def ens_from_mat_t(
    ens: Tuple[str, Dict[str, np.ndarray]]
) -> Tuple[str, np.ndarray]:
    """
    creates an "ensemble" dataframe from a dictionary of wholes where
    each "whole" is a numpy 2D array.

    Parameters
    ----------
    ens: tuple of np.ndarray
        A tuple in which the first element is an ensemble name and the second
        one is a dictionary in which  keys are whole names and values are
        whole-type data.
    Return
    ------
    A tuple of ensemble name and its associated dataframe.
    """
    return ens[0], np.stack(list(ens[1].values()), axis=0)


def ens_from_df(
    ens: Tuple[str, Dict[str, pd.DataFrame]]
) -> Tuple[str, pd.DataFrame]:
    """
    creates an "ensemble" dataframe from a dictionary of wholes where
    each "whole" is a pandas dataframe.

    In each "whole" dataframe, headers are "elements" of a "matrix" or 2D
    quantity and columns are the values of a given property.

    A "whole" dataframe is of two types: "histogram" or "timeseries". A "whole"
    "histogram" dataframe has an additional header "bin_center" that contains
    the values of bin centers.

    There is a stark difference between ensembles of "vectors" and "matrices"
    whole types and ensembles of "dataframe" whole type. In the former, the
    headers in an "ensemble" dataframe are "whole" names and the columns are
    the values of a given properties. In the latter, the headers in an
    "ensemble" dataframe are "elements" of a matrix and the columns are the
    values of "elements", each averaged over all the "whole" dataframes.

    Note
    ----
    It is assumed that performing "mean" (see below) over other headers than
    "element" headers does not change those headers. "bin_center" header is an
    example of such mean-invariant headers.

    Parameters
    ----------
    ens: tuple of pd.DataFrame
        A tuple in which the first element is an ensemble name and the second
        one is a dictionary in which keys are whole names and values are
        whole-type data.

    Return
    ------
    A tuple of ensemble name and its associated dataframe.
    """
    return ens[0], pd.concat(list(ens[1].values())).groupby(level=0).mean()


def ensemble(
    property_: str,
    wholes: Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]],
    parser: ParserT,
    geometry: str,
    group: str,
    topology: str,
    whole_type: str,
    edge_wholes: Optional[Dict[str, np.ndarray]] = None,
    save_to: Optional[str] = None,
) -> Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]]:
    """Generates ensembles from `wholes` for the physical property `property_`
    of a particle `group` in a `geometry` of interest.

    The `whole_type` can be "vector" or "1D" numpy array, "matrix" or "2D"
    numpy array, or "dataframe". The "dataframe" `whole_type` are "whole"
    dataframe in which the headers are elements of a matrix and the columns
    are either the values of that elements over time (a "timeseries"
    "dataframe" `whole_type`) or a measurement on values; for instance, the
    measurement can be counting (histogram). For such a "histogram"
    "dataframe" `whole_type`, there is an additional header that is
    "bin_center". The length of "timeseries" "whole" "dataframe" is equal to
    the number of time frames while the length of "histogram" "whole"
    "dataframe" is equal to the number of bins.

    There is a stark difference between ensembles of "vectors" and "matrices"
    whole types and ensembles of "dataframe" whole type. In the former, the
    headers in an "ensemble" dataframe are "whole" names and the columns are
    the values of a given properties. In the latter, the headers in an
    "ensemble" dataframe are "elements" of a matrix and the columns are the
    values of "elements", each averaged over all the "whole" dataframes; as
    a result, ensembles of "dataframe" `whole-type` are ensemble-averaged
    and do not any further steps.

    Parameters
    ----------
    property_ : str
        The physical property.
    wholes : WholeT
        A dictionary in which keys are 'whole' names and values are data
        (1D array).
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    geometry : {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        The type of the particle group.
    topology: str
        Topology of the polymer.
    whole_type: {'vector', 'matrix', 'dataframe', 'bin_edge'}
        The type of "whole" files:

        'vector':
            A numpy 1D array; for example, a time series or a histogram.
        'matrix':
            A numpy 2D array; for example, the gyration matrix.
        'dataframe':
            A pandas dataframe; for example, the pair distances of a group of
            monomers.
        'bin_edge':
            A bin_edge 1D array.

    edge_wholes:  dict of np.ndarray, default None
        A dictionary in which keys are 'whole' names and values are bin_edges.
        This option is used if `wholes` are histograms. Since the bin edges (
        and thus bin centers) are the same for all the "whole" histogram in an
        ensemble, it is written just once as a column to the final "ensemble"
        dataframe. `edge_wholes` is meaningful if "vector" `whole_type` is
        used.
    save_to : str, default None
        Absolute or relative path of a directory to which outputs are
        saved.

    Return
    ------
    ensembles : EnsembleT
        A dictionary for `property_` in which keys are ensemble names and
        values are dataframes. In each dataframe, the columns are wholes of
        that ensemble.
    """
    # Averaging over ensembles with similar initial parameters
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(whole_type, ["vector", "matrix", "dataframe", "bin_edge"])
    ensembles: Dict[str, Dict[str, Union[pd.DataFrame, np.ndarray]]] = {}
    bin_centers = {}
    for w_name, w_arr in wholes.items():
        w_info = parser(
            w_name, "whole", geometry, group, topology, ispath=False
        )
        ens_name = getattr(w_info, "ensemble_long")
        if not bool(ensembles):  # is ens_names empty or not?
            ensembles[ens_name] = {w_name: w_arr}
        elif ens_name not in ensembles.keys():
            ensembles[ens_name] = {w_name: w_arr}
        else:
            ensembles[ens_name].update({w_name: w_arr})
        if edge_wholes is not None:
            bin_centers[ens_name] = 0.5 * (
                edge_wholes[w_name][:-1] + edge_wholes[w_name][1:]
            )
    whole_types = {
        "vector": {"mapping_func": ens_from_vec, "ext": "csv"},
        "matrix": {"mapping_func": ens_from_mat_t, "ext": "npy"},
        "dataframe": {"mapping_func": ens_from_df, "ext": "csv"},
        "bin_edge": {"mapping_func": ens_from_bin_edge, "ext": "npy"},
    }
    ensembles_merged = dict(
        map(whole_types[whole_type]["mapping_func"], ensembles.items())
    )
    if edge_wholes is not None:
        for ens_name in ensembles_merged.keys():
            ensembles_merged[ens_name]["bin_center"] = bin_centers[ens_name]
    if save_to is not None:
        _ = dict(
            map(
                lambda ens: (
                    ens[0],
                    save_parent(
                        ens[0],
                        ens[1],
                        property_,
                        save_to,
                        group
                    ),
                ),
                ensembles_merged.items(),
            )
        )
    return ensembles_merged


def ens_avg_from_df(
    ens_prop: str, ens_data: pd.DataFrame, exclude: List[str]
) -> pd.DataFrame:
    """
    Creates an "ensAvg" dataframe from an "ensemble" dataframe. The columns
    in the "ensemble" dataframe `ens_data` are the "whole" data and any other
    variable/data given by `exclude`; for instance, if the ensembles is a
    histogram, then there is "bin_center" column in addition to the "whole"
    columns.

    Parameters
    ----------
    ens_prop: str
        The property name to which theis "ensemble" data belongs.
    ens_data: pd.DataFrame
        The dataframe of "ensemble"data.
    exclude: list of str
        The list of columns other than "whole" columns.

    Return
    ------
    ens_data: pd.DataFrame
        Update the `ens_data` with averaging statistics.
    """
    wholes = [col for col in ens_data.columns if col not in exclude]
    ens_data[ens_prop + "-mean"] = ens_data[wholes].mean(axis=1)
    ens_data[ens_prop + "-var"] = ens_data[wholes].var(axis=1, ddof=1)
    ens_data[ens_prop + "-sem"] = ens_data[wholes].sem(axis=1, ddof=1)
    ens_data.drop(columns=wholes, inplace=True)
    return ens_data


def ens_avg_from_bin_edge(
    ens_prop: str, ens_data: np.ndarray, exclude: List[str]
) -> np.ndarray:
    """
    Not written yet
    """
    return np.unique(ens_data)


def ens_avg_from_ndarray(
    ens_prop: str, ens_data: np.ndarray, exclude: List[str]
) -> Dict[str, np.ndarray]:
    """
    Creates an "ensAvg" matrix from an "ensemble" matrix. The first axis of
    the "ensemble" matrix (axis=0 in numpy lingo) contains the whole matrices
    and statistical measurements are performed on this axis.

    Parameters
    ----------
    ens_prop: str
        The property name to which theis "ensemble" data belongs.
    ens_data: pd.DataFrame
        The matrix of "ensemble" data.
    exclude: list of str
        The list of columns other than "whole" columns. It does not have any
        application but is defined for consistency with `ens_avg_from_df`.

    Return
    ------
    ens_data: dict
        A dictionary in which values of various statistical measures and
        values are the ensemble-averages.
    """
    size_sqrt = ens_data.shape[0] ** 0.5
    ens_avg = {
        ens_prop + "-mean": np.mean(ens_data, axis=0),
        ens_prop + "-var": np.var(ens_data, axis=0, ddof=1),
        ens_prop + "-sem": np.std(ens_data, axis=0, ddof=1) / size_sqrt,
    }
    return ens_avg


def ensemble_avg(
    property_: str,
    ensembles: Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]],
    geometry: str,
    group: str,
    ens_type: str,
    exclude: Optional[list] = None,
    save_to: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """
    Performa averaging over the "whole" data in each "ensemble" data in the
    `ensembles` of the physical property `property_` of a particle `group`, if
    that columns is a valid 'whole' simulation name, not a `exclude` columns.

    The typ

    Parameters
    ----------
    property_ : str
        The physical property.
    ensembles : dict of DataFrame
        Dictionary of siblings (dataframes) where keys are parent names and
        values are siblings (dataframes). In each siblings' dataframe, the
        number of columns is equal to or more than the number of
        children of that parent. Columns' names are the children' names and
        the names given by `exclude`.
    geometry : {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    ens_type: {'vector', 'matrix'}
        The type of "ens" values.

        'dataframe':
            A dataframe in which each column is either a "whole" or an item
            from `exclude` list.

        'ndarray':
            A ndarray in which the elements along the first axis (axis=0 in
            numpy lingo) are "whole" matrices. The length of ndarray along the
            first axis is equal to the number of different wholes.

    exclude: list of str, default ['bin_center']
        List of columns that are not 'whole' or 'segment' simulation names.
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    ens_avgs : dict of pd.DataFrame
        Dict of  on `property_` where keys are ensemble names and
        values are dataframes of ensemble-averaged measurements.
    """
    # Averaging over ensembles with similar initial parameters
    if exclude is None:
        exclude = ["bin_center"]
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(ens_type, ["dataframe", "ndarray", "bin_edge"])
    ens_avgs = {}
    ens_types = {
        "dataframe": {"ens_avg_func": ens_avg_from_df, "ext": "csv"},
        "ndarray": {
            "ens_avg_func": ens_avg_from_ndarray,
            "ext": "dict_of_npy",
        },
        "bin_edge": {"ens_avg_func": ens_avg_from_bin_edge, "ext": "npy"},
    }
    for ens, ens_data in ensembles.items():
        ens_avg = ens_types[ens_type]["ens_avg_func"](
            property_, ens_data, exclude
        )
        ens_avgs[ens] = ens_avg
    if save_to is not None:
        property_ = property_ + "-ensAvg"
        _ = dict(
            map(
                lambda ens_dummy: (
                    ens_dummy[0],
                    save_parent(
                        ens_dummy[0],
                        ens_dummy[1],
                        property_,
                        save_to,
                        group
                    ),
                ),
                ens_avgs.items(),
            )
        )
    return ens_avgs


def children_stamps(
    stamps: List[Tuple[str]],
    group: str,
    lineage: str,
    save_to: Optional[str] = None,
) -> pd.DataFrame:
    """
    Generates a dataset of the physical attributes and equilibrium
    (time-averaged) physical properties of all the `lineage` simulations of a
    particle `group` in a 'space' in a `geometry`.

    The name of space is created from the first value of 'filename' in
    the generated dataset.

    Parameters
    ----------
    stamps: list of tuple
        List of tuples where each tuple has one member and that member is a
        filepath to the stamp of a 'segment' simulation in a space.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    lineage: {'segment', 'whole'}
        Lineage type of child stamps
    save_to : str, default None
        Absolute or relative path of a directory to which outputs are saved.

    Return
    ------
    space_stamps: pd.DataFrame
        All the child stamps in the `group` in a space.
    """
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(lineage, ["segment", "whole"])
    space_stamps = [pd.read_csv(stamp[0]) for stamp in stamps]
    space_stamps = pd.concat(space_stamps)
    space_stamps.reset_index(inplace=True, drop=True)
    cols_to_drop = ["segment", "segment_id"]
    if lineage == "whole":
        # Some older version of parsers use 'segment'= "N/A" and
        # 'segment_id'="N/A" in a "whole" stamp when "whole" linage
        # is used.
        try:
            space_stamps.drop(columns=cols_to_drop, inplace=True)
            warnings.warn(
                "'segment' and 'segment_id' columns are dropped when"
                " individual 'whole' stamps combined to create a single"
                " dataframe of 'whole' stamps by 'children_stamps'.",
                UserWarning,
            )
        except KeyError:
            print(f"'{cols_to_drop}' are not among columns.")
    if save_to is not None:
        space_name = str(space_stamps.loc[0, "space"])
        filename = "-".join([space_name, group, lineage, "stamps.csv"])
        space_stamps.to_csv(save_to + filename, index=False)
    return space_stamps


def parents_stamps(
    stamps: pd.DataFrame,
    geometry: str,
    group: str,
    lineage: str,
    properties: Optional[Dict[str, Callable]] = None,
    save_to: Optional[str] = None,
) -> pd.DataFrame:
    """
    Performs merging/ensemble-averaging over all the 'segment/'whole'
    simulation stamps in a 'space' in a given `geometry` for a given
    `group` based on the given `lineage`.

    Parameters
    ----------
    stamps: DataFrame
        Dataframe of all the simulation stamps in the `group` in a space.
    geometry : str in {'cylindrical', 'slit', 'cubic'}
        Shape of the simulation box
    group: str in {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    lineage: str in  {'segment', 'whole'}
        Lineage type of children's stamps.
    properties: dict of str
        A dictionary in which the keys are properties such as the
        time-averaged radius of gyration which are measured during the 'probe'
        phase and the values are user-defined or numpy functions which are
        used as the aggregation function by pandas.
    save_to : str, default None
        Absolute or relative path of a directory to which outputs are saved.

    Notes
    -----
    If `lineage='segment'`, then stamps are for 'segments' and they have only
    different 'segment_id' for a given 'whole' parent. If `lineage='whole'`,
    then stamps are for 'wholes' and they have only different 'ensemble_id'
    for a given 'ensemble' parent. In either scenarios, the 'stamp' files
    have 'segment' and 'segment_id' columns. If `lineage='whole'`, the
    values of these two columns are "N/A".

    Return
    ------
    stamps_avg: pd.DataFrame
        Dataframe of all the parents stamps in the `group` in a space.
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(lineage, ["segment", "whole"])
    # attributes, properties and genealogy:
    stamps_cols = list(stamps.columns)
    try:
        stamps_cols.remove("lineage_name")
        stamps_cols.remove(lineage)
    except ValueError:
        print(
            f"'lineage_name' and '{lineage}'"
            " columns are not among in stamps column:"
            f"'{stamps_cols}', they are probably removed in"
            " a previous call of 'parents_stamps' function."
        )
    # aggregation dictionary: See Note above.
    agg_funcs = dict()
    attr_agg_funcs = ["last"] * len(stamps_cols)
    agg_funcs.update(zip(stamps_cols, attr_agg_funcs))
    if properties is not None:  # add/update agg funcs for properties.
        agg_funcs.update(properties)
    # handling lineage-specific details:
    if lineage == "whole":
        parent_groupby = "ensemble_long"
        # aggregating functions for properties
        agg_funcs["ensemble_id"] = "count"
        agg_funcs["n_frames"] = "last"
        file_lastname = "ensAvg"
    else:
        parent_groupby = "whole"
        # aggregating functions for properties
        agg_funcs["segment_id"] = "count"
        agg_funcs["n_frames"] = "sum"
        file_lastname = "whole"
    agg_funcs.pop(parent_groupby)
    parent_stps = stamps.groupby([parent_groupby]).agg(agg_funcs)
    parent_stps.reset_index(inplace=True)
    if lineage == "whole":
        parent_stps.rename(
            columns={"ensemble_id": "n_ensembles"}, inplace=True
        )
        # If the 'whole' stamps are generated directly in the 'probe' phase,
        # then 'segment' and 'segment_id' columns are "N/A" and are removed
        # from the list of stamps columns that are added to the parents
        # stamps.
        # There is no need to have the "n_segment" column in the parent
        # stamps, so it is removed. The "whole" stamps directly generated
        # in the "probe" phase do not have such a column, but those generated
        # from "segment" stamps have.
        # Dropping redundant columns silently:
        parent_stps.drop(
            columns=["n_segments", "segment_id", "segment"],
            inplace=True,
            errors="ignore",
        )
    else:
        parent_stps.rename(columns={"segment_id": "n_segments"}, inplace=True)
    if save_to is not None:
        space_name = str(parent_stps.loc[0, "space"])
        filename = "-".join(
            [space_name, group, "stamps", file_lastname + ".csv"]
        )
        parent_stps.to_csv(save_to + filename, index=False)
    return parent_stps


def unique_property(
    filepath: str,
    prop_idx: int,
    extensions: List[str],
    drop_properties: Optional[List[str]] = None,
    sep: Optional[str] = "-",
) -> Tuple[List[str], List[str]]:
    """
    Finds unique physical properties and physical property-measures by
    splitting filenames given by the glob-friendly 'filepath'.

    A measure refers to some measurement done on a physical property.

    A physical property-measure is defined as a measurement done on a
    physical;
    for instance, "gyrT-acf" means the auto-correlation function of the radius
    of gyration.

    Parameters
    ----------
    filepath: str
        The globe-friendly filepath.
    prop_idx: int
        The index after which a property name or property-measure name starts.
    extensions: list of ext
        The extensions that comes after "property" or "property-measure" name
        such as "-ensAvg", "-ens", or "-whole"
        This is different from a file's extensions/format such as "csv" or
        "npy".
    drop_properties: list of str, default None
        The properties that should be ignored.
    sep: str, default "-"
        The seperator between a "property" and its "measure".

    Return
    ------
    uniq_props: list of str
        A sorted list of unique physical properties.
    uniq_prop_measures: list of str
        A sorted list of unique property-measures.
    """
    props_measures = glob(filepath)
    uniq_prop_measures = []
    for ext in extensions:
        prop_measure_per_ext = list(
            set(
                [
                    sep.join(
                        prop.split("/")[-1].split(ext)[0].split(sep)[prop_idx:]
                    )
                    for prop in props_measures
                ]
            )
        )
        uniq_prop_measures.extend(prop_measure_per_ext)
    if drop_properties is not None:
        for drop_property in drop_properties:
            try:
                uniq_prop_measures.remove(drop_property)
            except ValueError:
                print(f"'{drop_property}' is not among unique properties.")
    uniq_props = set(
        [property_.split(sep)[0] for property_ in uniq_prop_measures]
    )
    uniq_prop_measures = set(uniq_prop_measures)
    uniq_prop_measures = list(uniq_prop_measures.difference(uniq_props))
    uniq_props = list(uniq_props)
    uniq_prop_measures.sort()
    uniq_props.sort()
    return uniq_props, uniq_prop_measures


def space_tseries(
    input_database: str,
    property_: str,
    parser: ParserT,
    hierarchy: str,
    physical_attrs: List[str],
    group: str,
    geometry: str,
    topology: str,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False,
) -> pd.DataFrame:
    """
    Takes the `property_path` to 'ensAvg' time series of a given `property_`
    in a given space `input_database`,  adds the `physical_attrs` of interest
    as the new columns to each 'ensAvg' dataframe, and merges all the 'ensAvg'
    dataframes into one 'space' dataframe along the 0 (or 'row' or 'index')
    in pandas lingo,

    In each 'ensemble-averaged' dataframe, there are 3 columns with this name
    pattern:

    column name = '[long_ensemble]-[property_][-measure]-[stat]'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'.

    Parameters
    ----------
    input_database: str
        Path to the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    hierarchy: str
        The pattern by which the filenames of timeseries are started with; for
        instance, "N*" means files start with "N"
    physical_attrs: list of str
        The physical attributes that will be added as new columns to the
        concatenated timeseries.
    group: str in {'bug', 'nucleoid', 'all'}
        The type of the particle group.
    geometry : str in {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    topology: str
        Topology of the polymer.
    divisor: float, default 0.025
        The step by which the values of "phi_c_bulk" attribute are rounded.
    round_to: int, default 3
        The number of significant decimal digits in the values of "phi_c_bulk"
        attribute.
    is_save : bool, default False
        whether to save output to file or not.

    Return
    ------
    all_in_one: pandas.DataFrame
        a dataframe in which all the timeseries are concatenated along `orient`
        of interest, and "properties and attributes" of interest are added to
        it as the new columns.
    """
    property_ext = "-" + property_ + "-ensAvg.csv"
    ens_avg_csvs = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs, fmts=[property_ext])
    property_db = []
    parser_name = getattr(parser, "__name__", "unknown")
    if parser_name == "unknown":
        raise ValueError(f"'{parser}' does not have a name!")

    dumping_freq = {
        "TransFociCyl": "bdump",
        "TransFociCub": "bdump",
        "SumRuleCubHeteroLinear": "bdump",
        "SumRuleCubHeteroRing": "bdump",
        "SumRuleCyl": "bdump",
        "HnsCub": "ndump",
        "HnsCyl": "ndump",
    }

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        property_info = parser(
            ens_avg_csv[0], "ensemble_long", geometry, group, topology
        )
        ens_avg.reset_index(inplace=True)
        ens_avg.rename(columns={"index": "t_index"}, inplace=True)
        ens_avg["t_index"] = ens_avg["t_index"] * getattr(
            property_info, dumping_freq[parser_name]
        )
        ens_avg["time"] = ens_avg["t_index"] * property_info.dt
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)
        ens_avg["phi_c_bulk_round"] = ens_avg["phi_c_bulk"].apply(
            round_up_nearest, args=[divisor, round_to]
        )
        property_db.append(ens_avg)
    property_db = pd.concat(property_db, axis=0)
    property_db.reset_index(inplace=True, drop=True)
    if is_save is not False:
        save_to_space = database_path(
            input_database, "analysis", stage="space", group="bug"
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        output = "-".join([space, group, property_]) + "-space.csv"
        property_db.to_csv(save_to_space + output, index=False)
    return property_db


def space_hists(
    input_database: str,
    property_: str,
    parser: ParserT,
    hierarchy: str,
    physical_attrs: List[str],
    group: str,
    geometry: str,
    topology: str,
    bin_center: Optional[np.ndarray] = None,
    normalize: Optional[bool] = False,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False,
) -> pd.DataFrame:
    """
    Takes the `property_path` to 'ensAvg' time series of a given `property_`
    in a given space `input_database`,  adds the `physical_attrs` of interest
    as the new columns to each 'ensAvg' dataframe, and merges all the 'ensAvg'
    dataframes into one 'space' dataframe along the 0 (or 'row' or 'index')
    in pandas lingo,

    In each 'ensemble-averaged' dataframe, there are 4 columns with this name
    pattern:

    column name = '[long_ensemble]-[property_][-measure]-[stat]'

    , and sometimes

    column name = 'bin_center'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'. If the 'bin_center' presents as a column in a
    'ensemble_averaged' dataframe, then it is inferred; otherwise, it should be
    passed to the function. See `bin_center` kw argument below.

    Parameters
    ----------
    input_database: str
        Path to the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    hierarchy: str
        The pattern by which the filenames of timeseries are started with; for
        instance, "N*" means files start with "N"
    physical_attrs: list of str
        The physical attributes that will be added as new columns to the
        concatenated timeseries.
    group: {'bug', 'nucleoid', 'all'}
        The type of the particle group.
    geometry : {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    topology: str
        Topology of the polymer.
    bin_center: numpy array, default None
        The bin centers. The argument should be given if the 'bin_center' is
        not in the ensemble-averaged dataframe and the same array of bin
        centers is used in all different ensemble-averaged dataframes. This is
        the case for "clustersHistTFoci" or "bondsHistTFoci" properties, but
        not for "zHistMon" or "rHistCrd".
    normalize: bool, default False
        Whether normalize frequencies or not.
    divisor: float, default 0.025
        The step by which the values of "phi_c_bulk" attribute are rounded.
    round_to: int, default 3
        The number of significant decimal digits in the values of "phi_c_bulk"
        attribute.
    is_save : bool, default False
        whether to save output to file or not.

    Return
    ------
    all_in_one: pandas.DataFrame
        a dataframe in which all the timeseries are concatenated along `orient`
        of interest, and "properties and attributes" of interest are added to
        it as the new columns.

    Requirements
    ------------
    PolyPhys, Pandas
    """
    property_ext = "-" + property_ + "-ensAvg.csv"
    ens_avg_csvs = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs, fmts=[property_ext])
    property_db = []
    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        property_info = parser(
            ens_avg_csv[0], "ensemble_long", geometry, group, topology
        )
        if bin_center is not None:
            ens_avg["bin_center"] = bin_center.tolist()
            ens_avg["bin_center-norm"] = (
                ens_avg["bin_center"] / ens_avg["bin_center"].max()
            )
        if normalize is True:
            normalizer = ens_avg[property_ + "-mean"].sum()
            if normalizer != 0:
                ens_avg[property_ + "-norm"] = (
                    ens_avg[property_ + "-mean"] / normalizer
                )
            else:
                warnings.warn(
                    "All the frequencies are zero, so all the normalized"
                    " frequencies are set to zero.",
                    UserWarning,
                )
                ens_avg[property_ + "-norm"] = 0
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)
        ens_avg["phi_c_bulk_round"] = ens_avg["phi_c_bulk"].apply(
            round_up_nearest, args=[divisor, round_to]
        )
        property_db.append(ens_avg)
    property_db = pd.concat(property_db, axis=0)
    property_db.reset_index(inplace=True, drop=True)
    if is_save is not False:
        save_to_space = database_path(
            input_database, "analysis", stage="space", group=group
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        output = "-".join([space, group, property_]) + "-space.csv"
        property_db.to_csv(save_to_space + output, index=False)
    return property_db


def normalize_z(
    prop: str, ens_avg: pd.DataFrame, norm_direction: bool = True
) -> pd.DataFrame:
    """Normalizes the ensemble-average local distribution `ens_avg` of `prop`
    along z direction in cylindrical geometry by the maximum value of the
    `ens_avg`, and takes average over the absolute values of bin centers along
    z direction if `norm_direction` is `True`.

    Parameters
    ----------
    prop: str
        Name of the physical property.
    ens_avg: pd.DataFrame
        Ensemble-average local distribution.
    norm_direction: bool, default True
        Whether averaging over absolute values of bin_centers or not.

    Return
    ------
    ens_avg or esn_sum: pd.DataFrame
        Normalized ensemble-average local distribution.
    """
    ens_avg_max = ens_avg[prop + "-scale"].max()
    ens_avg[prop + "-normalizer"] = ens_avg_max
    if ens_avg_max != 0:
        ens_avg[prop + "-norm"] = ens_avg[prop + "-scale"] / ens_avg_max
    else:
        warnings.warn(
            "All the frequencies are zero, so all the normalized"
            " frequencies are set to zero.",
            UserWarning,
        )
        ens_avg[prop + "-norm"] = 0
    # If system is symmetric with respect to z=0, then an average can be
    # applied with respect to absolut size of bin centers.
    if norm_direction is True:
        df_len = ens_avg.shape[0]
        if df_len % 2 == 0:
            # index or bin_center or z >= 0:
            ens_pos = ens_avg.iloc[df_len // 2 :, :].copy()
            ens_pos.reset_index(inplace=True, drop=True)
            # index or bin_center or z < 0:
            ens_neg = ens_avg.iloc[: df_len // 2, :].copy()
            ens_neg.index = -1 * ens_neg.index
            ens_neg.sort_index(inplace=True)
            ens_neg.reset_index(inplace=True, drop=True)
            ens_neg["bin_center"] = -1 * ens_neg["bin_center"]
            # averaging over |z|>0:
            ens_sum = 0.5 * (ens_pos + ens_neg)
        else:
            # index or bin_center or z > 0
            ens_pos = ens_avg.iloc[df_len // 2 + 1 :, :].copy()
            ens_pos.reset_index(inplace=True, drop=True)
            # index or bin_center or z < 0
            ens_neg = ens_avg.iloc[: df_len // 2, :].copy()
            ens_neg.index = -1 * ens_neg.index
            ens_neg.sort_index(inplace=True)
            ens_neg.reset_index(inplace=True, drop=True)
            ens_neg["bin_center"] = -1 * ens_neg["bin_center"]
            # averaging over |z|>0:
            ens_sum = 0.5 * (ens_pos + ens_neg)
            # index or bin_center or z = 0
            ens_sum.set_index("bin_center", inplace=True)
            ens_sum.loc[0, :] = ens_avg.loc[df_len // 2, :]
            ens_sum.sort_index(inplace=True)
            ens_sum.reset_index(inplace=True)
        return ens_sum
    else:
        return ens_avg


def normalize_r(
    prop: str, ens_avg: pd.DataFrame, method: str = "first"
) -> pd.DataFrame:
    """Normalizes the ensemble-average local distribution `ens_avg` of `prop`
    along r direction in cylindrical geometry by the `method` chosen for the
    value of thr normalizer.

    Parameters
    ----------
    prop: str
        Name of the physical property.
    ens_avg: pd.DataFrame
        Ensemble-average local distribution.
    method: {'first', 'max'}, default 'first'
        Normalization method
        'first':
            Normalizing by the first value of `ens_avg`.
        'max':
            Normalizing by the maximum of `ens_avg`.

    Return
    ------
    ens_avg: pd.DataFrame
        Normalized ensemble-average local distribution.
    """
    if method == "first":
        ens_avg_max = ens_avg[prop + "-scale"].values[0]
    elif method == "max":
        ens_avg_max = ens_avg[prop + "-scale"].max()
    else:
        raise NotImplementedError("Choose either 'first' or 'max' method.")
    ens_avg[prop + "-normalizer"] = ens_avg_max
    if ens_avg_max != 0:
        ens_avg[prop + "-norm"] = ens_avg[prop + "-scale"] / ens_avg_max
    else:
        warnings.warn(
            "All the frequencies are zero, so all the normalized"
            " frequencies are set to zero.",
            UserWarning,
        )
        ens_avg[prop + "-norm"] = 0
    return ens_avg


def space_sum_rule(
    input_database: str,
    property_: str,
    parser: ParserT,
    hierarchy: str,
    physical_attrs: List[str],
    species: str,
    size_attr: str,
    group: str,
    geometry: str,
    topology: str,
    direction: str,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False,
) -> pd.DataFrame:
    """Takes the `property_path` to 'ensAvg' local distribution of a given
    `property_` in a given space `input_database`, normalize and scale that
    distribution, adds the `physical_attrs` of interest as the new columns to
    each 'ensAvg' distribution, and merges all the 'ensAvg' distributions into
    one 'space' dataframe along the 0 (or 'row' or 'index') in pandas lingo,

    In each 'ensemble-averaged' dataframe, there are 4 columns with this name
    pattern:

    column name = '[long_ensemble]-[property_][-measure]-[stat]'

    , and sometimes

    column name = 'bin_center'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'. If the 'bin_center' presents as a column in a
    'ensemble_averaged' dataframe, then it is inferred; otherwise, it should
    be passed to the function. See `bin_center` kw argument below.

    A few words on terminology:
    - 'scaling' means changing the range of data while 'normalizing' means
    changing the distribution of data.
    - For a given particle species, the column "prop + '-scaler'" is set to the
    particle size (particle size squared) for volume fraction (number density)
    profiles.


    Parameters
    ----------
    input_database: str
        Path to the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    hierarchy: str
        The pattern by which the filenames of timeseries are started with; for
        instance, "N*" means files start with "N"
    physical_attrs: list of str
        The physical attributes that will be added as new columns to the
        concatenated timeseries.
    species: {'Mon', 'Crd', 'Foci', 'Dna'}
        The species of the particles in a group.
            'Mon': Monomers or small monomers
            'Crd': Crowders
            'Foci': Large monomers
            'Dna': Small and large monomers
    size_attr: str
        The attribute of the `parser` object that is the size (diameter) of
        species.
    group: {'bug', 'nucleoid', 'all'}
        The type of the particle group.
    geometry : {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    topology: str
        Topology of the polymer.
    direction: {'r', 'z'}
        The direction along which operation is done.
    divisor: float, default 0.025
        The step by which the values of "phi_c_bulk" attribute are rounded.
    round_to: int, default 3
        The number of significant decimal digits in the values of "phi_c_bulk"
        attribute.
    is_save : bool, default False
        whether to save output to file or not.

    Return
    ------
    all_in_one: pandas.DataFrame
        a dataframe in which all the timeseries are concatenated along `orient`
        of interest, and "properties and attributes" of interest are added to
        it as the new columns.

    Requirements
    ------------
    PolyPhys, Pandas
    """
    normalizer = {"r": normalize_r, "z": normalize_z}
    property_ext = "-" + group + "-" + direction + property_ + species
    property_ext += "-ensAvg.csv"
    prop = direction + property_ + species  # full name of physical property
    ens_avg_csvs = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs, fmts=[property_ext])
    property_db = []
    # ens_csvs is a list of tuples, each has one member.
    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        property_info: ParserT = parser(
            ens_avg_csv[0], "ensemble_long", geometry, group, topology
        )
        if property_ == "Phi":
            scaler = getattr(property_info, size_attr)
            ens_avg[prop + "-scaler"] = scaler
            ens_avg[prop + "-scale"] = ens_avg[prop + "-mean"] / scaler
        elif property_ == "Rho":
            scaler = getattr(property_info, size_attr)
            ens_avg[prop + "-scaler"] = scaler**2
            ens_avg[prop + "-scale"] = ens_avg[prop + "-mean"] * scaler**2
        else:
            raise NotImplementedError(
                "Sum rule's scaler is only defined for "
                "'rho' (density) or 'phi' (volume fraction) properties."
            )
        ens_avg[prop + "-scale-normalized_curve"] = (
            ens_avg[prop + "-scale"] / ens_avg[prop + "-scale"].sum()
        )
        # this is the part the sumr-rule is calculated:
        if direction == "r" and geometry == "cubic":
            ens_avg = normalizer[direction](prop, ens_avg, method="max")
        else:
            ens_avg = normalizer[direction](prop, ens_avg)
        ens_avg[prop + "-sumrule_constant"] = ens_avg[prop + "-scaler"]
        ens_avg["bin_center-norm"] = (
            ens_avg["bin_center"] / ens_avg["bin_center"].max()
        )
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)
        ens_avg["bin_center-dcrowd"] = (
            2 * ens_avg["bin_center"] / ens_avg["dcrowd"]
        )
        ens_avg["phi_c_bulk_round"] = ens_avg["phi_c_bulk"].apply(
            round_up_nearest, args=[divisor, round_to]
        )
        if geometry == "cylindrical":
            ens_avg["temp"] = (ens_avg["dcyl"] % ens_avg["dcrowd"]) / (
                ens_avg["dcrowd"]
            )
            ens_avg["bin_center-dcrowd-recentered"] = (
                ens_avg["bin_center-dcrowd"] - ens_avg["temp"]
            )
            ens_avg["bin_center-recentered-norm"] = ens_avg["bin_center"] - (
                ens_avg["dcyl"] % ens_avg["dcrowd"]
            )
            ens_avg["bin_center-recentered-norm"] = (
                ens_avg["bin_center-recentered-norm"]
                / ens_avg["bin_center-recentered-norm"].max()
            )
            ens_avg.drop(columns=["temp"], inplace=True)
        property_db.append(ens_avg)
    property_db = pd.concat(property_db, axis=0)
    property_db.reset_index(inplace=True, drop=True)
    if is_save is not False:
        save_to_space = database_path(
            input_database, "analysis", stage="space", group=group
        )
        space = save_to_space.split("/")[-2].split("-")[0]
        output = "-".join([space, group, property_, species])
        output += "-normalizedRescaled-space.csv"
        property_db.to_csv(save_to_space + output, index=False)
        print("done")
    return property_db
