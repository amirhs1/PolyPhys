"""This oduel organize fa hierarchy of files.
"""

import warnings
from typing import List, Dict, Tuple, Optional, Union, Callable, Literal
import os
import pathlib
from glob import glob
import numpy as np
import pandas as pd

from ..manage.typer import ParserT, TransFociT
from ..analyze.clusters import whole_dist_mat_foci
from .utilizer import split_alphanumeric, round_up_nearest, invalid_keyword


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
        [
            "simAll",
            "simCont",
            "log",
            "trj",
            "probe",
            "analysis",
            "viz",
            "galaxy",
        ],
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


def whole_from_segments(
    property_: str,
    segments: List[Tuple[str]],
    parser: ParserT,
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    topology: Literal["ring", "linear", "branched"],
    relation: Literal["histogram", "tseries", "bin_edge"],
    save_to: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """
    Generate "whole" arrays by combining data segments for a specified
    property of a particle `group` within a defined `geometry`. The combination
    method is determined by the `relation` parameter.

    Parameters
    ----------
    property_ : str
        The physical property to process (e.g., "density").
    segments : List[Tuple[str]]
        List of tuples where each tuple contains at least one segment file path
        (e.g., to CSV files) for the `property_`.
    parser : ParserT
        A parser instance from `PolyPhys.manage.parser` module that interprets
        file paths or names and extracts attributes such as the `whole` name.
    geometry : str
        The simulation box geometry. Must be one of {'cylindrical', 'slit',
        'cubic'}.
    group : str
        Type of the particle group. Must be one of {'bug', 'nucleoid', 'all'}.
    topology : str
        The polymer topology. Must be one of {'linear', 'ring', 'branched'}.
    relation : str
        Specifies how to combine segments into a whole. Accepted values:

        - 'histogram'
            'segemnt' is an N-dimensional histogram-like file, so the segments
            of a 'whole' should be sum along "axis=0" in numpy lingo. In the
            two-dimensional space, a histogram has a (x_nbins, y_nbins) shape.

        - 'tseries'
            'segment' is a time-series-like file, so the segments of a 'whole'
            should be concatenated vertically (along "axis=0" in pandas lingo).

        - 'bin_edge'
            'segment' is a bin-edge file. All the the segments are similar,
            so one 'whole' bin-edge file is created by picking one of
            segements.

    save_to : str, optional
        Directory path where the output files are saved, if specified.

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary with whole names as keys and the combined whole data as
        numpy arrays.

    Raises
    ------
    ValueError
        If `geometry`, `group`, or `relation` have invalid values.

    Notes
    -----
    Refer to the 'organizer' module documentation for additional context on the
    terms and keywords such as "geometry" and "group".

    Examples
    --------
    >>> whole_data = whole_from_segment(
            property_='density',
            segments=[('path/to/segment1.npy',), ('path/to/segment2.npy',)],
            parser=my_parser_instance,
            geometry='cubic',
            group='bug',
            topology='linear',
            relation='histogram',
            save_to='/output/path'
        )
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(topology, ["linear", "ring", "branched"])
    invalid_keyword(relation, ["histogram", "tseries", "bin_edge"])

    # Define combination strategies for each relation type
    mapping_func = {
        "histogram": lambda whole: (whole[0], np.sum(whole[1], axis=0)),
        "tseries": lambda whole: (whole[0], np.concatenate(whole[1])),
        "bin_edge": lambda whole: (whole[0], np.unique(np.array(whole[1]))),
    }

    # Combine segments into wholes based on segment names and relation type
    grouped_segments: Dict[str, List[np.ndarray]] = {}
    for seg_path in segments:
        seg_info = parser(seg_path[0], "segment", geometry, group, topology)
        whole_name = getattr(seg_info, "whole")
        segment = np.load(seg_path[0])

        if whole_name not in grouped_segments:
            grouped_segments[whole_name] = [segment]
        else:
            grouped_segments[whole_name].append(segment)

    # Apply the appropriate mapping function to merge each group of segments
    # into a whole.
    wholes: Dict[str, np.ndarray] = dict(
        map(mapping_func[relation], grouped_segments.items())
    )

    if save_to:
        for whole_name, whole_data in wholes.items():
            whole_fullname = create_fullname(whole_name, group, property_)
            save_parent(whole_fullname, whole_data, save_to)

    return wholes


def whole_from_file(
    whole_paths: List[Tuple[str]],
    parser: ParserT,
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    topology: Literal["ring", "linear", "branched"],
) -> Dict[str, np.ndarray]:
    """
    Loads data from "whole" files for a specified physical property of a
    particle `group` within a given `geometry`. Each file path in `whole_paths`
    points to a single whole file.

    Parameters
    ----------
    whole_paths : List[Tuple[str]]
        List of tuples where each tuple contains the path to a single file
        representing a whole dataset for the `property_`. Each file is loaded
        and processed independently.
    parser : ParserT
        A parser instance from `PolyPhys.manage.parser` module that interprets
        file paths and extracts attributes such as the `whole` name.
    geometry : {"cylindrical", "slit", "cubic"}
        Shape of the simulation box, indicating spatial configuration.
    group : {"bug", "nucleoid", "all"}
        The particle group to which the data pertains.
    topology : {"ring", "linear", "branched"}
        The polymer topology associated with the data.

    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary where keys are the `whole` names (derived from each file
        path)and values are numpy arrays containing the loaded data.

    Raises
    ------
    ValueError
        If `geometry`, `group`, or `topology` contain invalid values.

    Notes
    -----
    This function assumes that each file in `whole_paths` contains data related
    to a complete simulation or experiment segment, referred to as a "whole."
    For further information on terminology, see the `organizer` module
    documentation.

    Examples
    --------
    Load and organize whole data for a specific property from a list of file
    paths:

    >>> whole_data = whole_from_file(
            whole_paths=[('path/to/whole1.npy',), ('path/to/whole2.npy',)],
            parser=my_parser_instance,
            geometry='cubic',
            group='bug',
            topology='linear'
        )
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(topology, ["linear", "ring", "branched"])

    wholes: Dict[str, np.ndarray] = {}
    for whole_path in whole_paths:
        whole_info = parser(whole_path[0], "whole", geometry, group, topology)
        whole_name = getattr(whole_info, "whole")
        wholes[whole_name] = np.load(whole_path[0])

    return wholes


def whole_from_dist_mat_t(
    whole_paths: List[Tuple[str]],
    parser: TransFociT,
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    topology: Literal["ring", "linear", "branched"],
) -> Tuple[
    Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]
]:
    """
    Loads pair distance matrix data for "whole" files from given file paths,
    processing specific distance-related physical properties (e.g.,
    frequencies, radial distribution functions, time series). Each data type
    is extracted and stored in separate dictionaries by whole name.

    Parameters
    ----------
    whole_paths : List[Tuple[str]]
        List of tuples where each tuple contains the path to a file
        representing a whole dataset.
    parser : TransFociT
        A parser instance that processes file paths to extract attributes such
        as the whole name and distance matrix information.
    geometry : {"cylindrical", "slit", "cubic"}
        Shape of the simulation box, indicating spatial configuration.
    group : {"bug", "nucleoid", "all"}
        The particle group to which the data pertains.
    topology : {"ring", "linear", "branched"}
        The polymer topology associated with the data.

    Returns
    -------
    Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame],
          Dict[str, pd.DataFrame]]
        Three dictionaries with whole names as keys and DataFrames as values:

        - First dictionary contains frequency data.
        - Second dictionary contains radial distribution functions (RDF).
        - Third dictionary contains time series data.

    Raises
    ------
    ValueError
        If `geometry`, `group`, or `topology` have invalid values.

    Notes
    -----
    This function assumes that each file in `whole_paths` contains data in a
    distance matrix format, which is typically used in molecular simulation
    analysis for spatial or time-based properties.

    Examples
    --------
    Load distance matrix data for multiple physical properties:

    >>> freqs, rdfs, tseries = whole_from_dist_mat_t(
            whole_paths=[('path/to/whole1.npy',), ('path/to/whole2.npy',)],
            parser=my_parser_instance,
            geometry='cubic',
            group='bug',
            topology='linear'
        )
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(topology, ["linear", "ring", "branched"])

    wholes_freqs: Dict[str, pd.DataFrame] = {}
    wholes_rdfs: Dict[str, pd.DataFrame] = {}
    wholes_tseries: Dict[str, pd.DataFrame] = {}
    for whole_path in whole_paths:
        whole_info = parser(whole_path[0], "whole", geometry, group, topology)
        whole_freqs, whole_rdfs, whole_tseries = whole_dist_mat_foci(
            whole_path[0], whole_info
        )
        whole_name = getattr(whole_info, "whole")
        wholes_freqs[whole_name] = whole_freqs
        wholes_rdfs[whole_name] = whole_rdfs
        wholes_tseries[whole_name] = whole_tseries
    return wholes_freqs, wholes_rdfs, wholes_tseries


def ens_from_bin_edge(
    ens_data: Tuple[str, Dict[str, np.ndarray]]
) -> Tuple[str, np.ndarray]:
    """
    Combines multiple bin edge arrays into a single unique array for an
    "ensemble" by removing duplicates. This is typically used when bin edges
    are identical across "whole" datasets in an ensemble.

    Parameters
    ----------
    ens_data : Tuple[str, Dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name (str) and the
        second element is a dictionary of bin edge arrays from different
        "whole" datasets within that ensemble. Whitin the dictionary, keys are
        whole names while values are whole bin edges.

    Returns
    -------
    Tuple[str, np.ndarray]
        A tuple where the first element is the ensemble name, and the second
        is a 1D numpy array of unique bin edges.

    Examples
    --------
    >>> ensemble_bin_edges = ens_from_bin_edge(
            ('ensemble1', {'whole1': np.array([0.1, 0.2, 0.2, 0.3],
             'whole2': np.array([0.1, 0.2, 0.2, 0.3]})
        )
    >>> print(ensemble_bin_edges)
    ('ensemble1', array([0.1, 0.2, 0.3]))
    """
    ensemble_name, whole_edges = ens_data
    unique_bin_edges = np.unique(list(whole_edges.values()))
    return ensemble_name, unique_bin_edges


def ens_from_vec(
    ens_data: Tuple[str, Dict[str, np.ndarray]]
) -> Tuple[str, pd.DataFrame]:
    """
    Converts a dictionary of 1D arrays (representing different "whole"
    datasets) into a single ensemble-level DataFrame. This is useful when each
    "whole" dataset in the ensemble contains a vector (1D array) of
    measurements.

    Parameters
    ----------
    ens_data : Tuple[str, Dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name, and the second
        is a dictionary mapping each "whole" dataset name to a 1D numpy array
        of measurements.

    Returns
    -------
    Tuple[str, pd.DataFrame]
        A tuple where the first element is the ensemble name, and the second
        is a pandas DataFrame. Each column represents a "whole" dataset, and
        rows contain corresponding vector values from each dataset.

    Examples
    --------
    >>> ens_name, ensemble_df = ens_from_vec(
            ('ensemble1', {'whole1': np.array([0.1, 0.2, 0.3]),
                           'whole2': np.array([0.4, 0.5, 0.6])})
        )
    >>> print(ens_name)
    'ensemble1'
    >>> print(ensemble_df)
       whole1  whole2
    0     0.1     0.4
    1     0.2     0.5
    2     0.3     0.6
    """
    ensemble_name, whole_vectors = ens_data
    ensemble_df = pd.DataFrame.from_dict(whole_vectors, orient="columns")
    return ensemble_name, ensemble_df


def ens_from_mat_t(
    ens_data: Tuple[str, Dict[str, np.ndarray]]
) -> Tuple[str, np.ndarray]:
    """
    Combines multiple 2D arrays into a single 3D array for a specified
    "ensemble", where each 2D array corresponds to a "whole" dataset in the
    ensemble.

    Parameters
    ----------
    ens_data : Tuple[str, Dict[str, np.ndarray]]
        A tuple where the first element is the ensemble name (str) and the
        second element is a dictionary with "whole" names as keys and 2D numpy
        arrays as values.

    Returns
    -------
    Tuple[str, np.ndarray]
        A tuple where the first element is the ensemble name and the second
        is a 3D numpy array formed by stacking the 2D arrays along the first
        axis.

    Examples
    --------
    >>> ens_data = (
            "ensemble1",
            {"whole1": np.array([[1, 2], [3, 4]]),
            "whole2": np.array([[5, 6], [7, 8]])}
        )
    >>> combined_array = ens_from_mat_t(ens_data)
    >>> print(combined_array)
    ('ensemble1', array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]))
    """
    ensemble_name, whole_matrices = ens_data
    combined_matrix = np.stack(list(whole_matrices.values()), axis=0)
    return ensemble_name, combined_matrix


def ens_from_df(
    ens_data: Tuple[str, Dict[str, pd.DataFrame]]
) -> Tuple[str, pd.DataFrame]:
    """
    Averages multiple DataFrames in an ensemble by grouping them into a single
    averaged DataFrame. Each DataFrame corresponds to a "whole" dataset in the
    ensemble.

    In a "whole" DataFrame, headers are "elements" of a "matrix" or 2D
    quantity and columns are the values of a given property.

    In each "whole" DataFrame, headers represent "elements" of a 2D quantity
    (matrix-like structure), and columns contain the values for a specified
    property. Each "whole" DataFrame may be either a "histogram" or a
    "time-series":

    - In a "histogram" DataFrame, an additional column, `bin_center`, holds
      bin center values that are invariant under averaging and are retained
      in the output DataFrame. The length of histogram DataFrame is equal to
      the number of bins.

    - In a "time-series" DataFrame, headers (columns) are elements of a matrix,
      and rows represent the values of each element over time. The length of
      timeseries DataFrame is equal to the number of time frames.


    Parameters
    ----------
    ens_data : Tuple[str, Dict[str, pd.DataFrame]]
        A tuple where the first element is the ensemble name and the second
        is a DataFrame containing the mean of all input DataFrames, grouped
        by index.

    Returns
    -------
    Tuple[str, pd.DataFrame]
        A tuple with:
            - The ensemble name (str).
            - A DataFrame with the averaged values across all "whole"
            DataFrames.

    Notes
    -----
    A clear distinction exists between ensembles of "vector" or "matrix" types
    and ensembles of "DataFrame" type. In vector or matrix types, each column
    in the resulting ensemble DataFrame is named after a "whole" and represents
    that whole's values for a specific property. In contrast, for DataFrame
    type ensembles, the columns are elements of a 2D quantity (e.g., matrix or
    time-series values) and hold averaged values across all whole DataFrames.

    This function assumes that averaging non-element headers, such as
    "bin_center," does not change their values, and these are carried forward
    without modification.

    Examples
    --------
    >>> df1 = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
    >>> df2 = pd.DataFrame({"A": [2, 3, 4], "B": [5, 6, 7]})
    >>> ens_data = ("ensemble1", {"whole1": df1, "whole2": df2})
    >>> averaged_df = ens_from_df(ens_data)
    >>> print(averaged_df)
    ('ensemble1',
       A    B
    0  1.5  4.5
    1  2.5  5.5
    2  3.5  6.5)
    """
    ensemble_name, whole_data = ens_data
    ensemble_df = pd.concat(list(whole_data.values())).groupby(level=0).mean()
    return ensemble_name, ensemble_df


def ensemble(
    property_: str,
    wholes: Dict[str, Union[np.ndarray, pd.DataFrame]],
    parser: ParserT,
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    topology: Literal["ring", "linear", "branched"],
    whole_type: Literal["vector", "matrix", "DataFrame", "bin_edge"],
    whole_edges: Optional[Dict[str, np.ndarray]] = None,
    save_to: Optional[str] = None,
) -> Dict[str, Union[np.ndarray, pd.DataFrame]]:
    """
    Generates an ensemble by  merging "whole" data arrays or DataFrames
    representing a specified `property_` of a particle `group` within a
    defined `geometry`. The `whole_type` determines the structure of each
    "whole" and the merging approach.

    Parameters
    ----------
    property_ : str
        The physical property for which the ensemble is generated (e.g.,
        "density").
    wholes : Dict[str, Union[np.ndarray, pd.DataFrame]]
        Dictionary containing "whole" names as keys and their corresponding
        data as values, represented as numpy arrays or DataFrames.
    parser : ParserT
        An instance of the parser class to infer attributes like ensemble name
        from each "whole".
    geometry : {"cylindrical", "slit", "cubic"}
        The geometry of the simulation box.
    group : {"bug", "nucleoid", "all"}
        Type of the particle group being processed.
    topology : {"ring", "linear", "branched"}
        The polymer topology.
    whole_type : {"vector", "matrix", "DataFrame", "bin_edge"}
        Specifies the type of each "whole" and defines the merging method:

        - "vector": A 1D numpy array, such as a time series or histogram.
        - "matrix": A 2D numpy array, such as a gyration tensor.
        - "DataFrame": A DataFrame representing either a time series (with
          headers as matrix elements and rows as values over time) or a
          histogram (headers as elements, rows as counts or frequencies,
          and an additional "bin_center" column for bin centers).
        - "bin_edge": A 1D array representing bin edges. Assumes all bin edges
          are identical and selects one to represent the ensemble.

    whole_edges : Dict[str, np.ndarray], optional
        Dictionary with "whole" names as keys and bin edge arrays as values,
        useful when `whole_type` is "vector" and bin edges are shared across
        wholes in the ensemble.
    save_to : str, optional
        Directory path where output files are saved, if specified.

    Returns
    -------
    Dict[str, Union[np.ndarray, pd.DataFrame]]
        Dictionary where keys are ensemble names and values are
        ensemble data as numpy arrays or DataFrames.

    Raises
    ------
    ValueError
        If any parameter has an invalid value.

    Notes
    -----
    Different averaging methods apply based on `whole_type`. When `whole_type`
    is "DataFrame," the resulting DataFrame averages the values of each matrix
    element across wholes, retaining headers that are mean-invariant (e.g.,
    "bin_center"). For "vector" and "matrix" types, columns represent
    individual wholes in the ensemble DataFrame.

    Examples
    --------
    >>> wholes = {
            "whole1": np.array([1, 2, 3]),
            "whole2": np.array([2, 3, 4])
        }
    >>> ensembles = ensemble(
            property_="density",
            wholes=wholes,
            parser=my_parser_instance,
            geometry="cubic",
            group="bug",
            topology="linear",
            whole_type="vector",
            save_to="/output/path"
        )
    """
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(whole_type, ["vector", "matrix", "DataFrame", "bin_edge"])

    merging_func: Dict[str, Callable] = {
        "vector": ens_from_vec,
        "matrix": ens_from_mat_t,
        "DataFrame": ens_from_df,
        "bin_edge": ens_from_bin_edge,
    }
    # Collecting wholes for each ensemble
    ensemble_wholes: Dict[str, Dict[str, Union[pd.DataFrame, np.ndarray]]] = {}
    bin_centers: Dict[str, np.ndarray] = {}
    if whole_edges is not None:
        for whole_name, whole_data in wholes.items():
            whole_info = parser(
                whole_name, "whole", geometry, group, topology, ispath=False
            )
            ensemble_name = getattr(whole_info, "ensemble_long")
            if ensemble_name not in ensemble_wholes:
                ensemble_wholes[ensemble_name] = {}
            ensemble_wholes[ensemble_name][whole_name] = whole_data
            bin_centers[ensemble_name] = 0.5 * (
                whole_edges[whole_name][:-1] + whole_edges[whole_name][1:]
            )
    else:
        for whole_name, whole_data in wholes.items():
            whole_info = parser(
                whole_name, "whole", geometry, group, topology, ispath=False
            )
            ensemble_name = getattr(whole_info, "ensemble_long")
            if ensemble_name not in ensemble_wholes:
                ensemble_wholes[ensemble_name] = {}
            ensemble_wholes[ensemble_name][whole_name] = whole_data

    # Merging wholesdata into a singel ensembel data object
    ensembles: Dict[str, Union[np.ndarray, pd.DataFrame]] = \
        dict(map(merging_func[whole_type], ensemble_wholes.items()))

    # Save ensemble data object if save_to is specified
    if save_to:
        for ensemble_name, data in ensembles.items():
            ensemble_fullname = create_fullname(
                ensemble_name, group, property_
            )
            save_parent(ensemble_fullname, data, save_to)

    return ensembles


def ens_avg_from_df(
    ens_prop: str, ens_data: pd.DataFrame, exclude: List[str]
) -> pd.DataFrame:
    """
    Creates an "ensAvg" DataFrame from an "ensemble" DataFrame. The columns
    in the "ensemble" DataFrame `ens_data` are the "whole" data and any other
    variable/data given by `exclude`; for instance, if the ensembles is a
    histogram, then there is "bin_center" column in addition to the "whole"
    columns.

    Parameters
    ----------
    ens_prop: str
        The property name to which theis "ensemble" data belongs.
    ens_data: pd.DataFrame
        The DataFrame of "ensemble"data.
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
        Dictionary of siblings (DataFrames) where keys are parent names and
        values are siblings (DataFrames). In each siblings' DataFrame, the
        number of columns is equal to or more than the number of
        children of that parent. Columns' names are the children' names and
        the names given by `exclude`.
    geometry : {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    group: {'bug', 'nucleoid', 'all'}
        Type of the particle group.
    ens_type: {'vector', 'matrix'}
        The type of "ens" values.

        'DataFrame':
            A DataFrame in which each column is either a "whole" or an item
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
        values are DataFrames of ensemble-averaged measurements.
    """
    # Averaging over ensembles with similar initial parameters
    if exclude is None:
        exclude = ["bin_center"]
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(ens_type, ["DataFrame", "ndarray", "bin_edge"])
    ens_avgs = {}
    ens_types = {
        "DataFrame": {"ens_avg_func": ens_avg_from_df, "ext": "csv"},
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
                        ens_dummy[0], ens_dummy[1], property_, save_to, group
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
                " DataFrame of 'whole' stamps by 'children_stamps'.",
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
    as the new columns to each 'ensAvg' DataFrame, and merges all the 'ensAvg'
    DataFrames into one 'space' DataFrame along the 0 (or 'row' or 'index')
    in pandas lingo,

    In each 'ensemble-averaged' DataFrame, there are 3 columns with this name
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
        a DataFrame in which all the timeseries are concatenated along `orient`
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
    as the new columns to each 'ensAvg' DataFrame, and merges all the 'ensAvg'
    DataFrames into one 'space' DataFrame along the 0 (or 'row' or 'index')
    in pandas lingo,

    In each 'ensemble-averaged' DataFrame, there are 4 columns with this name
    pattern:

    column name = '[long_ensemble]-[property_][-measure]-[stat]'

    , and sometimes

    column name = 'bin_center'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'. If the 'bin_center' presents as a column in a
    'ensemble_averaged' DataFrame, then it is inferred; otherwise, it should be
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
        not in the ensemble-averaged DataFrame and the same array of bin
        centers is used in all different ensemble-averaged DataFrames. This is
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
        a DataFrame in which all the timeseries are concatenated along `orient`
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
    one 'space' DataFrame along the 0 (or 'row' or 'index') in pandas lingo,

    In each 'ensemble-averaged' DataFrame, there are 4 columns with this name
    pattern:

    column name = '[long_ensemble]-[property_][-measure]-[stat]'

    , and sometimes

    column name = 'bin_center'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'. If the 'bin_center' presents as a column in a
    'ensemble_averaged' DataFrame, then it is inferred; otherwise, it should
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
        a DataFrame in which all the timeseries are concatenated along `orient`
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
