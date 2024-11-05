"""
The `organizer` module provides tools for organizing and combining data from
simulations or experiments, based on patterns in filenames and file structures.
It supports various stages of post-processing, including assembling segmented
data, grouping experiments into ensembles, consolidating ensembles into
spaces, and building galaxies for high-level analysis.

This module introduces several key concepts and terminology to structure
simulation output and maintain consistency across post-processing steps:

Terminology
-----------

- **Experiment**:
    A complete simulation or run using tools like LAMMPS. Each experiment may
    generate multiple data segments or complete datasets representing a
    particular physical property or system configuration.

- **Lineage**:
    A hierarchical categorization of data files and directories based on
    increasing levels of data aggregation. The main levels of lineage are:

    - **Segment**:
        Represents a chunk of data from a longer simulation. Segments are
        typically smaller time or spatial units that together make up a
        complete dataset or **whole**.

    - **Whole** (or "Instance"):
        A single, complete dataset from a simulation, created by combining
        multiple segments. Each whole represents a standalone data instance
        for a particular set of conditions.

    - **Ensemble**:
        A collection of multiple whole datasets with identical macroscopic
        properties but varying initial conditions (e.g., random seed values).
        These represent thermodynamic states of the same system, enabling
        ensemble averaging for statistical analysis.

    - **Space**:
        An aggregation of ensembles where all conditions are fixed except one
        macroscopic attribute (e.g., box size or temperature). Spaces allow
        for systematic variation of a single attribute and support comparative
        analysis across different conditions.

    - **Galaxy**:
        The broadest level, representing multiple spaces, each with varying
        attributes and system configurations. Galaxies enable high-level
        cross-attribute analyses, encompassing data from diverse experiments
        within a unified framework.

Attributes and Properties
-------------------------

- **Attribute**:
    A macroscopic characteristic of a system, such as box size, polymer
    topology, or particle type. Attributes remain fixed within each space,
    except for one that defines the specific space.

- **Property**:
    A physical feature measured during a simulation, such as density or
    radius of gyration. Properties may have associated measures (e.g.,
    autocorrelation or distribution) that detail the property's behavior or
    evolution.

- **Species**:
    Type of particle or molecular entity, such as "Mon" for monomers. The
    module defines specific behavior or attributes for each species as
    needed.

- **Direction**:
    Indicates a spatial direction, such as "r" for radial in spherical
    coordinates. Used to orient data appropriately based on coordinate
    system.

File Naming and Structure
-------------------------

Files and directories generated at different phases of post-processing are
organized according to the lineage and phase of the data they contain. This
module relies on consistent naming conventions to streamline data handling.

- **Phases**:
    Directories are organized into several phases to manage data from
    simulation to visualization:

    - `simulationsAll`:
        Contains all LAMMPS-related files for running simulations in a space.

    - `logs`, `trjs`:
        Store logs and trajectory data, respectively.

    - `probe`:
        Holds files created by directly probing simulation data, such as
        time-series or histograms.

    - `analyze`:
        Contains files resulting from analyzing or averaging data across
        ensembles, such as autocorrelation functions.

    - `viz`:
        Houses data files ready for visualization.

- **Stages**:
    Indicates the level of aggregation within each phase:

    - `wholeSim`:
        Combines all segments into a whole for each experiment.

    - `ens`:
        Collects whole files into an ensemble.

    - `ensAvg`:
        Holds ensemble-averaged files.

- **Naming Conventions**:
    Files and directories are named using consistent patterns, helping to
    identify data lineage, phase, and attributes. Naming is organized by
    the following patterns:

    - `lineage-phase-group-property_stage.extension`

    - `whole|segment.group.run.lammpstrj` (for trajectory files)

Each name pattern provides quick insight into the contents of the file or
directory, such as the type of data (e.g., segment or whole), its phase,
and physical properties.

Examples
--------

1. **Whole File**:
    A single dataset containing complete simulation results for a specific
    condition, stored under the `wholeSim` stage.

2. **Ensemble-Averaged Data**:
    An averaged dataset across all whole files in an ensemble, saved in
    `ensAvg`.

3. **Space Group**:
    A directory containing ensemble-averaged files with one varying attribute,
    useful for comparative analysis across different ensemble groups.

Each aggregation level (segment, whole, ensemble, space, galaxy) supports
specific functions in this module to streamline data organization, enabling
efficient post-processing of simulation results.

Attributes for Consistency
--------------------------

- `attribute`: A high-level system feature such as box size or topology.
- `property_`: A physical property measured during simulation, like density.
- `species`: Type of particle or entity in the simulation.
- `direction`: Spatial direction relevant to the property.
- `group`: A particle grouping in LAMMPS, like "bug".
- `phase`: Data processing stage, from "simulations" to "visualization".
- `stage`: Aggregation level, like "wholeSim" or "ensAvg".

Notes
-----
This module is optimized for use with LAMMPS outputs and general simulation
data. It includes functionalities to rename and organize files based on
metadata extracted from their paths or names, facilitating systematic data
exploration.
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
    and ensembles of "dataframe" type. In vector or matrix types, each column
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
    wholes: Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]],
    parser: ParserT,
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    topology: Literal["ring", "linear", "branched"],
    whole_type: Literal["vector", "matrix", "dataframe", "bin_edge"],
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
    wholes : Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]]
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
    whole_type : {"vector", "matrix", "dataframe", "bin_edge"}
        Specifies the type of each "whole" and defines the merging method:

        - "vector": A 1D numpy array, such as a time series or histogram.
        - "matrix": A 2D numpy array, such as a gyration tensor.
        - "dataframe": A DataFrame representing either a time series (with
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
        Dictionary where keys are ensemble names and values are ensemble data
        as numpy arrays or DataFrames.

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
    invalid_keyword(whole_type, ["vector", "matrix", "dataframe", "bin_edge"])

    merging_func: Dict[str, Callable] = {
        "vector": ens_from_vec,
        "matrix": ens_from_mat_t,
        "dataframe": ens_from_df,
        "bin_edge": ens_from_bin_edge,
    }
    # Collecting wholes for each ensemble
    ens_wholes: Dict[str, Dict[str, Union[pd.DataFrame, np.ndarray]]] = {}
    bin_centers: Dict[str, np.ndarray] = {}
    if whole_edges is not None:
        for whole_name, whole_data in wholes.items():
            whole_info = parser(
                whole_name, "whole", geometry, group, topology, ispath=False
            )
            ens_name = getattr(whole_info, "ensemble_long")
            if ens_name not in ens_wholes:
                ens_wholes[ens_name] = {}
            ens_wholes[ens_name][whole_name] = whole_data
            bin_centers[ens_name] = 0.5 * (
                whole_edges[whole_name][:-1] + whole_edges[whole_name][1:]
            )
    else:
        for whole_name, whole_data in wholes.items():
            whole_info = parser(
                whole_name, "whole", geometry, group, topology, ispath=False
            )
            ens_name = getattr(whole_info, "ensemble_long")
            if ens_name not in ens_wholes:
                ens_wholes[ens_name] = {}
            ens_wholes[ens_name][whole_name] = whole_data

    # Merging wholesdata into a singel ensembel data object
    ensembles: Dict[str, Union[np.ndarray, pd.DataFrame]] = \
        dict(map(merging_func[whole_type], ens_wholes.items()))

    # Save ensemble data object if save_to is specified
    if save_to:
        for ens_name, data in ensembles.items():
            ens_fullname = create_fullname(ens_name, group, property_)
            save_parent(ens_fullname, data, save_to)

    return ensembles


def ens_avg_from_bin_edge(ens_data: np.ndarray) -> np.ndarray:
    """
    Generates a unique set of bin edges for an "ensemble" dataset based on
    the bin edges in `ens_data`. This function assumes that the bin edges
    across the ensemble are identical and selects one representative set.

    Parameters
    ----------
    ens_data : List[np.ndarray]
        A 1D array of bin edges across multiple wholes in the ensemble. The
        function assumes these bin edges are identical across all wholes.

    Returns
    -------
    np.ndarray
        A unique, representative set of bin edges for the ensemble.

    Notes
    -----
    This function does not perform any averaging or statistical calculations.
    Instead, it identifies unique bin edges across all wholes, assuming they
    are identical and discards duplicates.

    Examples
    --------
    >>> ens_data = [np.array([0.0, 0.1, 0.2]), np.array([0.0, 0.1, 0.2])]
    >>> unique_bin_edges = ens_avg_from_bin_edge(ens_data)
    >>> print(unique_bin_edges)
    [0.  0.1 0.2 0.3]
    """
    return np.unique(ens_data)


def ens_avg_from_ndarray(
    ens_prop: str,
    ens_data: np.ndarray
) -> Dict[str, np.ndarray]:
    """
    Computes ensemble statistics (mean, variance, and standard error of the
    mean)across an "ensemble" 3D array where each slice along the first axis
    (axis=0) represents a "whole" data matrix. The statistics are calculated
    along this axis.

    Parameters
    ----------
    ens_prop : str
        The name of the property to which this "ensemble" data belongs. Used
        as a prefix for the computed statistics.
    ens_data : np.ndarray
        A 3D numpy array where each matrix slice along axis=0 represents a
        "whole" dataset.

    Returns
    -------
    Dict[str, np.ndarray]
        A dictionary where keys are statistical measures (mean, variance, SEM)
        with `ens_prop` as a prefix while values are numpy arrays of the
        computed statistics across "whole" matrices.

    Notes
    -----
    - The SEM is calculated with a degrees-of-freedom adjustment (`ddof=1`) for
      an unbiased estimate and is scaled by the square root of the number of
      wholes.
    - This function performs averaging along axis=0 for each element across
      the "whole" data matrices.

    Examples
    --------
    >>> ens_data = np.array([
            [[1, 2], [3, 4]],
            [[2, 3], [4, 5]]
        ])
    >>> ens_avg = ens_avg_from_ndarray("gyration_tensor", ens_data, exclude=[])
    >>> print(ens_avg)
    {
        'gyration_tensor-mean': array([[1.5, 2.5],
                               [3.5, 4.5]]),
        'gyration_tensor-var': array([[0.5, 0.5],
                              [0.5, 0.5]]),
        'gyration_tensor-sem': array([[0.5, 0.5],
                              [0.5, 0.5]])
    }
    """
    n_wholes_sqrt = ens_data.shape[0] ** 0.5
    ens_avg = {
        f"{ens_prop}-mean": np.mean(ens_data, axis=0),
        f"{ens_prop}-var": np.var(ens_data, axis=0, ddof=1),
        f"{ens_prop}-sem": np.std(ens_data, axis=0, ddof=1) / n_wholes_sqrt
    }
    return ens_avg


def ens_avg_from_df(
    ens_prop: str,
    ens_data: pd.DataFrame,
    exclude: List[str]
) -> pd.DataFrame:
    """
    Calculates the mean, variance, and standard error of the mean (SEM) for
    columns representing "whole" data within an "ensemble" DataFrame. The
    resulting DataFrame retains excluded columns (e.g., bin_center), while
    replacing "whole" columns with computed statistics.

    Parameters
    ----------
    ens_prop : str
        The name of the property to which this "ensemble" data belongs, used
        as a prefix for the generated statistics columns.
    ens_data : pd.DataFrame
        A DataFrame containing ensemble data. Each "whole" dataset is in its
        own column, and any columns specified in `exclude` remain unaffected.
    exclude : list of str
        A list of column names to exclude from averaging. For example,
        `exclude` might contain "bin_center" for histogram-type ensembles.

    Returns
    -------
    pd.DataFrame
        The input `ens_data` DataFrame with "whole" columns replaced by
        calculated statistics: mean, variance, and SEM of each row across the
        "whole" data columns.

    Notes
    -----
    - The function averages only the "whole" columns, excluding any specified
      in `exclude`.
    - The standard error of the mean (SEM) is computed with a
    degrees-of-freedom adjustment (ddof=1) for an unbiased estimate.

    Examples
    --------
    >>> import pandas as pd
    >>> ens_data = pd.DataFrame({
            "bin_center": [1, 2, 3],
            "whole1": [4, 5, 6],
            "whole2": [5, 6, 7]
        })
    >>> exclude = ["bin_center"]
    >>> ens_avg_from_df("density", ens_data, exclude)
       bin_center  density-mean  density-var  density-sem
    0           1           4.5        0.5       0.5
    1           2           5.5        0.5       0.5
    2           3           6.5        0.5       0.5
    """
    wholes = [col for col in ens_data.columns if col not in exclude]

    ens_data[ens_prop + "-mean"] = ens_data[wholes].mean(axis=1)
    ens_data[ens_prop + "-var"] = ens_data[wholes].var(axis=1, ddof=1)
    ens_data[ens_prop + "-sem"] = ens_data[wholes].sem(axis=1, ddof=1)

    ens_data.drop(columns=wholes, inplace=True)
    return ens_data


def ensemble_avg(
    property_: str,
    ensembles: Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]],
    geometry: Literal["cylindrical", "slit", "cubic"],
    group: Literal["bug", "nucleoid", "all"],
    ens_type: Literal["dataframe", "ndarray", "bin_edge"],
    exclude: Optional[List[str]] = None,
    save_to: Optional[str] = None,
) -> Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]]:
    """
    Generates ensemble-averaged data for a specified `property_` by averaging
    "whole" data in each "ensemble" DataFrame or array within `ensembles`.
    Columns listed in `exclude` are omitted from averaging.

    Parameters
    ----------
    property_ : str
        The physical property for which ensemble averages are calculated.
    ensembles : Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]]
        A dictionary where each key is an ensemble name and each value is
        either a DataFrame or ndarray representing "whole" data.
    geometry : {"cylindrical", "slit", "cubic"}
        The shape of the simulation box.
    group : {"bug", "nucleoid", "all"}
        Type of the particle group.
    ens_type : {"dataframe", "ndarray", "bin_edge"}
        The data format for each ensemble:

        - "dataframe": DataFrame format, with columns as "whole" simulations
          or items in `exclude`.
        - "ndarray": 3D ndarray format, with each slice along axis=0
          representing a "whole" matrix.
        - "bin_edge": 1D ndarray format for bin edges, with bin edges assumed
          to be identical across wholes.

    exclude : List[str], default None
        List of column names or items not included in ensemble averaging,
        such as "bin_center" for histograms.
    save_to : str, optional
        Path to the directory where the averaged ensemble data will be saved.

    Returns
    -------
    Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]]
        A dictionary where each key is an ensemble name, and each value is a
        DataFrame or ndarray of ensemble-averaged data.

    Raises
    ------
    ValueError
        If `geometry`, `group`, or `ens_type` contains invalid values.

    Examples
    --------
    >>> ensembles = {
            "ensemble1": pd.DataFrame({"whole1": [1, 2, 3],
                                       "whole2": [2, 3, 4]}),
            "ensemble2": pd.DataFrame({"whole1": [3, 4, 5],
                                       "whole2": [4, 5, 6]})
        }
    >>> avg_ensembles = ensemble_avg(
            property_="density",
            ensembles=ensembles,
            geometry="cubic",
            group="bug",
            ens_type="DataFrame",
            exclude=["bin_center"],
            save_to="/output/path"
        )
    """
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])

    if exclude is None:
        exclude = ['bin_center']

    # Compute ensemble averages
    ens_avgs: Union[Dict[str, np.ndarray], Dict[str, pd.DataFrame]] = {}
    if ens_type == "bin_edge":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_bin_edge(ens_data)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    elif ens_type == "ndarray":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_ndarray(property_, ens_data)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    elif ens_type == "dataframe":
        for ens_name, ens_data in ensembles.items():
            ens_avg = ens_avg_from_df(
                property_, ens_data, exclude)  # type: ignore
            ens_avgs[ens_name] = ens_avg  # type: ignore
    else:
        invalid_keyword(ens_type, ["dataframe", "ndarray", "bin_edge"])

    # Save averaged ensemble data if save_to path is provided
    if save_to is not None:
        save_property = property_ + "-ensAvg"
        for ens_name, data in ens_avgs.items():
            ens_fullname = create_fullname(ens_name, group, save_property)
            save_parent(ens_fullname, data, save_to)

    return ens_avgs


def children_stamps(
    stamps: List[Tuple[str]],
    lineage: Literal["segment", "whole"],
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """
    Combines individual "stamp" CSV files into a single DataFrame for a given
    `lineage` in a *space*. Optionally, saves the resulting
    DataFrame to a specified directory.

    Parameters
    ----------
    stamps : List[Tuple[str]]
        List of tuples where each tuple contains a single file path to a
        "stamp" CSV file.
    lineage : {"segment", "whole"}
        Lineage type, either "segment" or "whole".
    save_to : str, optional
        Directory path where the combined DataFrame is saved, if specified.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame containing data from all the "stamp" files.

    Raises
    ------
    ValueError
        If `group` or `lineage` contains invalid values.

    Notes
    -----
    - It is assumed that all the stamp files have belong to the same *group* in
      a given *space*.
    - If `lineage` is "whole", columns "segment" and "segment_id" are removed,
    if present, as they are redundant for whole-lineage data.

    Examples
    --------
    >>> stamps = [('path/to/stamp1.csv',), ('path/to/stamp2.csv',)]
    >>> df = children_stamps(stamps, lineage='whole')
    """
    invalid_keyword(lineage, ["segment", "whole"])

    # Load and concatenate all stamps
    space_stamps = pd.concat([pd.read_csv(stamp[0]) for stamp in stamps],
                             ignore_index=True)

    # Drop specific columns for "whole" lineage
    cols_to_drop = ['segment', 'segment_id']
    if lineage == "whole":
        try:
            space_stamps.drop(columns=cols_to_drop, inplace=True)
            warnings.warn(
                "'segment' and 'segment_id' columns are dropped when"
                " individual 'whole' stamps combined to create a single"
                " dataframe of 'whole' stamps by 'children_stamps'.",
                UserWarning
            )
        except KeyError:
            print(f"'{cols_to_drop}' are not among columns.")

    # Save the resulting DataFrame to a CSV if specified
    if save_to is not None:
        space_name = str(space_stamps.loc[0, "space"])
        group = str(space_stamps.loc[0, "group"])
        filename = f"{space_name}-{group}-{lineage}-stamps.csv"
        space_stamps.to_csv(os.path.join(save_to, filename), index=False)

    return space_stamps


def parents_stamps(
    stamps: pd.DataFrame,
    lineage: Literal["segment", "whole"],
    properties: Optional[Dict[str, Callable]] = None,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """
    Aggregates data from child "stamps" into a parent stamp by applying
    specified aggregation functions in a *space*. Optionally, saves the
    resulting DataFrame to a CSV file.

    Parameters
    ----------
    stamps : pd.DataFrame
        DataFrame containing child stamp data.
    lineage : {"segment", "whole"}
        Lineage type, either "segment" or "whole".
    properties : dict of str -> Callable, optional
        Dictionary specifying aggregation functions for specific properties.
    save_to : str, optional
        Directory path where the parent DataFrame is saved, if specified.

    Returns
    -------
    pd.DataFrame
        Aggregated parent DataFrame.

    Raises
    ------
    ValueError
        If `group`, or `lineage` contain invalid values.

    Notes
    -----
    - It is assumed that all the stamp files have belong to the same *group* in
      a given *space*.
    - If `lineage` is "segment", stamps correspond to individual segments,
      each with a unique `segment_id`.
    - If `lineage` is "whole", stamps correspond to whole simulations, each
      identified by an `ensemble_id`.

    Examples
    --------
    >>> df = pd.DataFrame({'segment': [1, 2], 'value': [10, 20]})
    >>> parent_df = parents_stamps(df, group="bug",
                                   lineage="segment")
    """
    invalid_keyword(lineage, ["segment", "whole"])

    # Base aggregation functions for columns
    base_columns = list(stamps.columns)
    try:
        base_columns.remove("lineage_name")
        base_columns.remove(lineage)
    except ValueError:
        print(
            f"'lineage_name' and '{lineage}'"
            " columns are not among in stamps column:"
            f"'{base_columns}', they are probably removed in"
            " a previous call of 'parents_stamps' function."
        )
    agg_funcs: Dict[str, Union[str, Callable]] = \
        {col: "last" for col in base_columns}
    if properties is not None:
        agg_funcs.update(properties)

    # Define grouping column and lineage-specific aggregations
    if lineage == "whole":
        parent_groupby = "ensemble_long"
        agg_funcs.update({"ensemble_id": "count", "n_frames": "last"})
    else:
        parent_groupby = "whole"
        agg_funcs.update({"segment_id": "count", "n_frames": "sum"})

    # Perform the groupby and aggregate
    agg_funcs.pop(parent_groupby)
    parent_stamps = stamps.groupby(parent_groupby).agg(agg_funcs).reset_index()

    # Rename and drop columns based on lineage
    if lineage == "whole":
        parent_stamps.rename(columns={"ensemble_id": "n_ensembles"},
                             inplace=True)
        # If the 'whole' stamps are generated directly in the 'probe' phase,
        # then 'segment' and 'segment_id' columns are "N/A" and are removed
        # from the list of stamps columns that are added to the parents
        # stamps.
        # There is no need to have the "n_segment" column in the parent
        # stamps, so it is removed. The "whole" stamps directly generated
        # in the "probe" phase do not have such a column, but those generated
        # from "segment" stamps have.
        # Dropping redundant columns silently:
        parent_stamps.drop(columns=["n_segments", "segment_id", "segment"],
                           inplace=True, errors="ignore")
    else:
        parent_stamps.rename(columns={"segment_id": "n_segments"},
                             inplace=True)

    if save_to is not None:
        space_name = str(parent_stamps.loc[0, "space"])
        group = str(parent_stamps.loc[0, "group"])
        file_suffix = "ensAvg" if lineage == "whole" else "whole"
        filename = f"{space_name}-{group}-stamps-{file_suffix}.csv"
        parent_stamps.to_csv(os.path.join(save_to, filename), index=False)

    return parent_stamps


def find_unique_properties(
    filepath: str,
    prop_idx: int,
    extensions: List[str],
    drop_properties: Optional[List[str]] = None,
    sep: str = "-"
) -> Tuple[List[str], List[str]]:
    """
    Extracts unique physical properties and property-measures from filenames
    matched by a glob pattern. The function identifies unique segments in
    filenames based on specified extensions and index position, then sorts and
    returns them.

    Parameters
    ----------
    filepath : str
        The glob-friendly pattern used to locate filenames, e.g.,
        `path/to/files/*`.
    prop_idx : int
        The index position in the filename where the property or
        property-measure name starts after splitting by the separator.
    extensions : List[str]
        A list of suffixes that indicate the end of a property or
        property-measure name (e.g., `"-ensAvg"`, `"-ens"`, `"-whole"`). These
        are not file extensions like `".csv"` or `".npy"`.
    drop_properties : Optional[List[str]], default None
        A list of property names to ignore when determining unique properties
        and property-measures.
    sep : str, default "-"
        The separator used to split properties and measures within a filename.

    Returns
    -------
    Tuple[List[str], List[str]]
        A tuple containing:
        - **uniq_props** (*List[str]*): A sorted list of unique physical
          properties.
        - **uniq_prop_measures** (*List[str]*): A sorted list of unique
          property-measures.

    Examples
    --------
    Given the following filenames:

    - `path/to/files/gyrT-acf-ensAvg.npy`
    - `path/to/files/gyrR-acf-ensAvg.npy`
    - `path/to/files/temp-ens.npy`

    The function can be used as follows:

    >>> find_unique_properties("path/to/files/*", prop_idx=0,
                               extensions=["-ensAvg", "-ens"])
    (['gyrT', 'fsdT'], ['gyrT-acf', 'fsdT-acf'])

    Notes
    -----
    - Ensure `prop_idx` accurately reflects the location in the filename where
      the property or measure name starts after splitting by `sep`.
    - This function assumes a consistent filename format where property-measure
      segments follow each other with defined separators and suffixes.
    """
    props_measures = glob(filepath)
    uniq_prop_measures = set()
    for ext in extensions:
        for prop in props_measures:
            property_name = \
                sep.join(
                    prop.split("/")[-1].split(ext)[0].split(sep)[prop_idx:])
            uniq_prop_measures.add(property_name)

    if drop_properties is not None:
        uniq_prop_measures.difference_update(drop_properties)

    # Extract unique properties (first segment before `sep`)
    # from property-measures
    uniq_props = {prop.split(sep)[0] for prop in uniq_prop_measures}

    # Remove any full properties that match unique property names from
    # prop-measures
    uniq_prop_measures.difference_update(uniq_props)

    # Sort and return the results as lists
    return sorted(uniq_props), sorted(uniq_prop_measures)


def space_tseries(
    input_database: str,
    property_: str,
    parser: ParserT,
    hierarchy: str,
    physical_attrs: List[str],
    group: Literal["bug", "nucleoid", "all"],
    geometry: Literal["cylindrical", "slit", "cubic"],
    topology: Literal["ring", "linear", "branched"],
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False,
) -> pd.DataFrame:
    """
    Aggregates ensemble-averaged time-series data for a specified physical
    property across multiple files, adds specified physical attributes as
    columns, and concatenates into a single DataFrame.

    Parameters
    ----------
    input_database : str
        Path to the directory containing time-series data files.
    property_ : str
        Name of the physical property of interest.
    parser : ParserT
        Parser class to infer file-specific information from filenames or
        paths.
    hierarchy : str
        Pattern prefix for the time-series filenames (e.g., `"N*"`).
    physical_attrs : List[str]
        List of physical attributes to add as new columns in the output
        DataFrame.
    group : {'bug', 'nucleoid', 'all'}
        Particle group type.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Simulation box geometry.
    topology : {}
        Polymer topology.
    divisor : float, default 0.025
        Rounding step for `phi_c_bulk` attribute.
    round_to : int, default 3
        Number of decimal places for rounding `phi_c_bulk` values.
    is_save : bool, default False
        If True, saves the concatenated DataFrame to a CSV file.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of all time-series with added physical
        attributes.

    Examples
    --------
    >>> df = space_tseries("path/to/database", "density", parser=SomeParser,
                           hierarchy="N*", physical_attrs=["dmon", "dcyl"],
                           group="all", geometry="cubic", topology="linear")
    >>> df.head()

    Notes
    -----
    - This function assumes the presence of `phi_c_bulk` attribute in the
      parser output.
    - Requires a parser class with methods to retrieve attribute information
      for each file.
    """
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(topology, ["ring", "linear", "branched"])

    property_ext = f"-{property_}-ensAvg.csv"
    ens_avg_csvs = sort_filenames(
        glob(input_database + hierarchy + property_ext), [property_ext])
    property_csvs = []

    parser_name = getattr(parser, "__name__", "unknown")
    if parser_name == "unknown":
        raise ValueError(f"'{parser}' does not have a name!")

    # Mapping of dumping frequency based on parser names
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
        property_info = \
            parser(ens_avg_csv[0], "ensemble_long", geometry, group, topology)

        ens_avg.reset_index(inplace=True)
        ens_avg.rename(columns={"index": "t_index"}, inplace=True)

        # Calculate `t_index` and `time` columns based on `dumping_freq`
        ens_avg["t_index"] *= getattr(property_info, dumping_freq[parser_name])
        ens_avg["time"] = ens_avg["t_index"] * property_info.dt

        # Add physical attributes
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)

        # Apply rounding to `phi_c_bulk`
        ens_avg["phi_c_bulk_round"] = \
            ens_avg["phi_c_bulk"].apply(
                round_up_nearest, args=[divisor, round_to])

        property_csvs.append(ens_avg)

    # Concatenate all time-series DataFrames
    property_db = pd.concat(property_csvs, axis=0)
    property_db.reset_index(drop=True, inplace=True)

    # Optionally save to file
    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group)
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = save_to_space + f"{space}-{group}-{property_}-space.csv"
        property_db.to_csv(filepath, index=False)

    return property_db


def space_hists(
    input_database: str,
    property_: str,
    parser: ParserT,
    hierarchy: str,
    physical_attrs: List[str],
    group: Literal["bug", "nucleoid", "all"],
    geometry: Literal["cylindrical", "slit", "cubic"],
    topology: Literal["ring", "linear", "branched"],
    bin_center: Optional[np.ndarray] = None,
    normalize: Optional[bool] = False,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False,
) -> pd.DataFrame:
    """
    Aggregates ensemble-averaged histogram data for a specified physical
    property across multiple files, normalizes data if specified, and
    concatenates into a single DataFrame.

    Parameters
    ----------
    input_database : str
        Path to the directory containing histogram data files.
    property_ : str
        Name of the physical property of interest.
    parser : ParserT
        Parser class to infer file-specific information from filenames or
        paths.
    hierarchy : str
        Pattern prefix for the histogram filenames (e.g., `"N*"`).
    physical_attrs : List[str]
        List of physical attributes to add as new columns in the output
        DataFrame.
    group : {'bug', 'nucleoid', 'all'}
        Particle group type.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Simulation box geometry.
    topology : str
        Polymer topology.
    bin_center : np.ndarray, optional
        Array of bin centers. If not provided, must be present in the
        DataFrames.
    normalize : bool, default False
        If True, normalizes the histogram data.
    divisor : float, default 0.025
        Rounding step for `phi_c_bulk` attribute.
    round_to : int, default 3
        Number of decimal places for rounding `phi_c_bulk` values.
    is_save : bool, default False
        If True, saves the concatenated DataFrame to a CSV file.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of all histograms with added physical
        attributes.

    Examples
    --------
    >>> df = space_hists("path/to/database", "density", parser=SomeParser,
                         hierarchy="N*", physical_attrs=["temperature"],
                         group="all", geometry="cylindrical",
                         topology="linear")
    >>> df.head()

    Notes
    -----
    - If `normalize` is True, histogram values will be scaled to sum to 1.
    - The `bin_center` should be provided if it is not available in the input
      data. When all the ensemble-averaged DataFrames have the same bin
      center, like "clustersHistTFoci" or "bondsHistTFoci" properties, the
      `bij_center` may not provided.
    """
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(topology, ["ring", "linear", "branched"])

    property_ext = f"-{property_}-ensAvg.csv"
    ens_avg_csvs_ungrouped = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs_ungrouped, fmts=[property_ext])
    property_csvs = []

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        property_info = \
            parser(ens_avg_csv[0], "ensemble_long", geometry, group, topology)

        # Handle bin_center if provided
        if bin_center is not None:
            ens_avg["bin_center"] = bin_center.tolist()
            ens_avg["bin_center-norm"] = \
                ens_avg["bin_center"] / ens_avg["bin_center"].max()

        # Normalize data if specified
        if normalize:
            normalizer = ens_avg[property_ + "-mean"].sum()
            if normalizer != 0:
                ens_avg[property_ + "-norm"] = \
                    ens_avg[property_ + "-mean"] / normalizer
            else:
                warnings.warn(
                    "All values are zero; normalized values set to zero.",
                    UserWarning)
                ens_avg[property_ + "-norm"] = 0

        # Add physical attributes
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)

        # Apply rounding to `phi_c_bulk`
        ens_avg["phi_c_bulk_round"] = ens_avg["phi_c_bulk"].apply(
            lambda x: round_up_nearest(x, divisor, round_to)
        )

        property_csvs.append(ens_avg)

    # Concatenate all histograms
    property_db = pd.concat(property_csvs, axis=0)
    property_db.reset_index(drop=True, inplace=True)

    # Optionally save to file
    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group)
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = save_to_space + f"{space}-{group}-{property_}-space.csv"
        property_db.to_csv(filepath, index=False)

    return property_db
