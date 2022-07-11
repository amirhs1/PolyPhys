"""A module that organizes and combines pandas dataframes based on patterns in
filenames.


The organizer module is for post-processing phase in which clean data files
are generated from experiments. This module uses the following terminology
to distinguish different files and directories from each other:

experiment:
    another name for a simulation or run on HPC using 'LAMMPS' or similar
    packages.

lineage: {'segment', 'whole', 'ensemble', 'space'}
    a keyword that shows whether a file or directory is a 'segment',
    'whole', 'ensemble', or 'space'.

    segment:
    The output of a long simulation or run is usually split into two or more
    segments. 'segment' refers to one of multiple chunks of a 'whole' file.

    whole:
        A complete file; a collection of 'segments'.

    ensemble:
        A collection of 'wholes' files that differs only in their initial
        conditions (e.g., random number seed).
        An ensemble is a group of themodynamically-equivalent experiments
        (more precisely states) of a system. In in-situ experimentation,
        such experiments or states or virtual copies are created by doing an
        experiment with different seudo-random generator seeds while keeping
        all the other attributes fixed. In other words, experiments in an
        ensemble differs only in their initial conditions (the random
        generator seed) in the phase space, but have the same attributes
        (macroscopic property). The 'name' of an ensmble does not have 'ens'
        'segment' attribute. An ensemble can be varied by changing ONE
        macroscopic attribute of a system at a time. A group of enembles with
        ONE varring attribute belongs to an 'space' group.

    space:
        A collection of ensembles.
        A 'space' is created by fixing all the attributes of the system under
        study, except one.

attribute:
    a macroscopic feature that characterizes a system under study; for
    instance, the length of a cubic simulation box or the toplogy of a
    polymer (see `parser` module lists of 'attribute' keywords of a system.
    For the list of valid attributes, please see `parser` module and the
    classes defined therein.

property_:
    a phyiscal feature that is measured in an experiment; for instance, the
    bin edges and frequencies of the local number denisty of a given species
    or the end-to-end size of a polymer over time. a property 'keyword' is
    camelCase or CamelCase, and it has information about the direction, name,
    and species of a physica property; for example, "rHistMon" means a
    histogram file in the radial direction for monomers.

species:
    The keyword for the type of a species; for instance, 'Mon' stands for
    monomers. A list of species is defined within a class or function that
    works differently based on the type of species.

direction:
    The short name of a direction in a given coordiante system; for instance,
    'r' for the radial direction in the spherical coordinate system.  A list
    of species is defined within a class or function that works differently
    based on the type of species.

dir_long:
    The long name of a direction; for instance, 'radial' for 'r'.  A list of
    species is defined within a class or function that works differently
    based on the type of species.

group:
    a 'group' keyword, such as 'bug', identifies a collection of atoms as
    they belong to a 'LAMMPS' group (see LAMMPS 'group' command). For the
    list of valid group keywords, see `parser` module and the classes
    defined therein.

run: {'restart', 'cont}
    a 'run' keyword shows whether a simulation trajectory file ('lammpstrj'
    or 'trj' extentions) is 'restarted' from a previously broken simulation
    or a 'continued' from a sucessfully finished simulation.

phase: {'simulationsAll', 'simulationsCont', 'logs', 'trjs', 'probe',
        'analyze', and 'viz'}
    A keyword ised the name of a directory and is related to one of the
    post-processing phases

    simulationsAll:
    The directory that contains all the Lammps-related files for running all
    the simulations in a 'space', and the bash scripts for sanity checks and
    organzing simulation files on a cluster.

    simulationCont:
    The directory that is similar to the 'simulationAll', but it is about a
    'space' group that is re-run or continued.

    logs:
    The directory that contains the LAMMPS logs in a 'space' group.

    trjs:
    The directory that contains trajectory and toplogy files of a 'space'.

    probe:
    A directory or phase that coantins the files resulted from direct probing
    of trajectories in a 'space'. Measuring a polymer's length over time (a
    time series) or counting the number of species on a spatial grid (a
    histogram) are two popular probing operations.

        observation:
            a file resulted from probing a trjectory.

    analyze:
    A directory or phase that contains the files resulted from analyzing
    observations in a 'space'. Measuring the auto-correlation function of
    the instantaneous size of a polymer or calculating the radial volume
    fraction from a histogram in the radial direction are popular examples.
    Ensemble or ensemble-averaged files are also generated in this phase.

stage: {'wholeSim', 'ens', 'ensAvg'}
    This keyword shows the stage to which a file belongs and is used in
    naming files or directories.

    wholeSim:
    A directory that contains all the 'whole' files in a 'space', where each
    'whole' is created by combining all the 'segment' files belonging to a
    'whole'.

    ens:
    A directory that contains all the 'ensemble' files of a 'space'.

    ensAvg:
    A 'esneAvg' directory contains all the ensemble-averaged files of a
    'space' while an 'ensAvg' file is created by averaging over all the
    'whole' data in an 'ensemble' file.

lineage_name:
    The unique sub-name of a filename which follows one of 'lineage' patterns
    define in a class in the `parser` module.

name:
    A file or directory 'name' is a sequence of 'lineage_name', 'property_',
    'phase', 'group', 'stage', or file 'format', separated by a splitter.

    Period spliiter ".":
    This spliiter is used in naming the files the trajectory ('lammpstrj'
    and 'trj' extentions) and topology ('data' extention) files that are
    generated by LAMMMPS package.

    Hypen spliiter "-":
    This splitter is used in all other file formats, except 'trajectory' and
    'topology' files.

    A file name has onne of the following general patterns:
        filenames:
        whole|segment.group.run.lammpstrj
        whole.group.data
        whole.log
        whole.txt
        lineage-phase-group-property_-stage.format
        allInOne-property_.format

        Directories:
        lineage-phase-group-stage

I   A  file or directory may  have all or part of the aboe period- or hypen-
    separated keywords. See below for more information about 'AllInOne' files.

Different files and directories ceated at different 'phases' to oragnize
'segment', 'whole', 'ensmeble', or 'ensAvg' files. Below, some general
categories are difined for this files and directories.

    'whole' directories:
    A 'whole' directory contains all the files that belong to a 'whole'
    simulation. The 'name' of a 'whole' directory has the same pattern as
    the 'whole' lineage pattern.

    'ens' files:
    See above for the definition of the 'ensemble' 'lineage' and 'ens' 'stage'.

    'ensAvg' files:
    A file that is created by averageing over all the 'whole' files in an
    ensemble. An ensemble-averaged file has 'ensAvg' keyword in its filename.
    An ensemble-averaged file is created by changing ONE attribute of the
    system at a time (similar to an ensemble file). A group of
    enemble-averaged files with ONE varring attribute belongs to an 'space'
    group.

    'space' directories:
    A directory contains all the files of all the 'properties' of a 'group'
    in a 'space' in a given 'stage' at a given 'phase'.
    A 'space' directory may or may not have the 'group' keyword in tis name.

    'all-in-one' files:
    A clean file (dataset) that contains all the measurements in all the
    'space' groups about a given 'property_'.

A 'space' group results in a curve in the 'viz' phase. An 'ensAvg' file gives
a data point. For a given 'property_', if there are N esmebles, each with M
'whole' files, then there are N ensemble-average groups and N*M 'whole' files
in the space group.
"""
from typing import (
    List,
    Dict,
    Tuple,
    Optional,
    Union,
    Callable
)
import pathlib
from glob import glob
import re
import numpy as np
import pandas as pd
import warnings

from .utilizer import round_up_nearest


def camel_case_split(
    word: str
) -> List[str]:
    """split a camleCase or CamelCase `str` to its constituent sub-strings.

    Parameters
    ----------
    word: str
        Word to be split.

    Return
    ------
    list: list of str

    Reference
    ---------
    https://stackoverflow.com/questions/29916065/how-to-do-camelcase-split-in-python
    """
    return re.findall(r'[A-Z]?[a-z]+|[A-Z]+(?=[A-Z]|$)', word)


def isfloat(
    string: str
) -> Union[float, str]:
    """convert `string` to float if possible, otherwise returns the
        original string.

    Parameters
    ----------
    string : str

    Return
    ------
    possibly_float : str or float
        the float-version of the string or the original string.

    """
    try:
        possibly_float = float(string)
    except ValueError:
        possibly_float = string
    return possibly_float


def sort_by_alphanumeric(
    alphanumeric: str
) -> List[Union[int, str, float]]:
    """split an `alphanumeric` into words and integers.

    Parameters
    ----------
    alphanumeric : char
        an alphanumeric string.

    Return
    ------
    mixed_alphanum : list of str, int, and float
        a mixed list of words, integers, and floats.
    """
    number_pattern = re.compile(r'(\d+\.*\d*)')  # a integer or float number
    words = number_pattern.split(alphanumeric)
    mixed_alphanum = \
        [int(word) if word.isdigit() else isfloat(word) for word in words]
    return mixed_alphanum


def sort_filenames(
    fnames: List[str],
    fmts: List[str] = ['data', 'lammpstrj'],
    report: bool = False
) -> List[str]:
    """return an alphanumerically sorted list of strings.


    Groups `fnames` with the same names into tuples where the length of
    tuples is the length of `fmts`. For instance, A LAMMPS molecular
    dynamics simulation output usually is composed of two type of files;
    a trajectroy (trj or lammpstrj extentions) constisting of several
    snapshots and a topology (data extention) containing information
    about simulation box, particle and bond types.

    Parameters
    ----------
    fnames : list of str
        a list of filenames.
    fmts : list of str or tuple, defualt=['data',('lammpstrj','trj')]
        a list of formats where each format can have one or more extensions
        (passed as a tuple).
    report : bool, default=False
        shows a report or not.

    Returns
    -------
    filenames_sorted: list of tuple
        a sorted list of tuples where each tuple has `len(formats)` filenames.

    """
    fnames_by_fmt = [None] * len(fmts)
    # a nested list where each sublist has all the files with the same
    # extension:
    for fmt_idx, fmt_exts in enumerate(fmts):
        fnames_by_fmt[fmt_idx] = \
            [fname for fname in fnames if fname.endswith(fmt_exts)]
    for fmt_idx, fnames_same_fmt in enumerate(fnames_by_fmt):
        fnames_by_fmt[fmt_idx] = \
            sorted(fnames_same_fmt, key=sort_by_alphanumeric)
    fnames_sorted = list(zip(*fnames_by_fmt))
    if report:
        print("Total number of files is ", len(fnames_sorted))
        print("Path to the first tuple of the sorted file: ", fnames_sorted[0])
    return fnames_sorted


def invalid_keyword(
    keyword: str,
    valid_keywords: List[str]
) -> None:
    """raise an error if `keyowrd` is not in `valid_keywords`.

    Parameters
    ----------
    keyword: str
        Name of the `keyword`
    valid_keywords: array of str
        Array of valid keywords
    """
    if keyword not in valid_keywords:
        raise ValueError(
            f"'{keyword}'"
            " is not a valid option. Please select one of "
            f"{valid_keywords} options.")


def save_parent(
    name: str,
    data: Union[np.ndarray, pd.DataFrame],
    property_: str,
    save_to: str,
    group: str = 'bug',
    ext: str = 'csv',
) -> None:
    """save the `data` to memory as a file with extension `ext`.

    Parameters
    ----------
    name: str
        Name of the dataframe
    df: pd.DataFrame
        Dataframe to be saved.
    property_ : str
        The physical property.
    save_to : str
        An/a absolute/relative path of a directory to which outputs are saved.
    group: {'bug', 'all'}, default 'bug'
        Type of the particle group.
    ext: {'csv', 'npy'}, defualt 'csv'
    """
    invalid_keyword(ext, ['csv', 'npy'])
    filename = "-".join([name, group, property_])
    if ext == 'csv':
        data.to_csv(save_to + filename + ".csv", index=False)
    else:
        np.save(save_to + filename + ".npy", data)


def database_path(
    input_database: str,
    phase: str,
    stage: Optional[str] = None,
    group: Optional[str] = None
) -> str:
    """create a `stage` directory for a `group` at a given `phase` in
    a `phase` directory. If the directory exists, raises error, and
    return the path of the existing directory. If the `phase`
    directory does not exist, `new_directory` creates it.

    The path of 'phase' directory is inferred from the `old_path` and
    is at the same level of the 'phase' level of the `old_path`.

    General hierarchy of the `input_path`:
        "root/parent1/.../parentN/old_phase/old_directory"
    where the old_directory and the new_directory created below have
    the similar pattern:

        old_directory: "space-old_phase-old_group-old_stage"

        old_directory: "space-old_phase-old_stage"
            if the 'group' of old database is not important.

        new_directory: "space-`phase`-`group`-`stage`"

    Parameters
    ----------
    input_database: str
        Path to directory.
    phase : {'simulationsAll', 'simulationsCont', 'logs', 'trjs', 'probe',
    'analysis', 'viz'}
        Name of the new phase
    stage: {'segment', 'wholeSim', 'ens', 'ensAvg', 'space'}, default None
        Stage of the new directory.
    group: {'bug', 'all'}, default None
        Type of the particle group.

    Return
    ------
    output_database: str
        the string equivalent of the path to a new directory.
    """
    invalid_keyword(phase,
                    ['simulationsAll', 'simulationsCont', 'logs',
                     'trjs', 'probe', 'analysis', 'viz']
                    )
    invalid_keyword(group, ['bug', 'all', None])
    invalid_keyword(stage,
                    ['segment', 'wholeSim', 'ens', 'ensAvg', 'space', None]
                    )
    old_path = pathlib.Path(input_database)  # PurePath object
    old_space_parts = old_path.parts
    # Space name of a directory
    space_name = old_space_parts[-1].split('*')[0].split('-')[0]
    # Directory full name
    dir_name = [
        part for part in [space_name, group, stage] if part is not None
    ]
    if len(dir_name) > 1:
        dir_name = '-'.join(dir_name)
    else:
        dir_name = dir_name[0]
    # Parents of old 'phase' directory, except 'root':
    output_database = list(old_space_parts[:-2])
    output_database.append(phase)  # New 'phase' path
    output_database.append(dir_name)  # New directory path
    output_database = '/'.join(output_database)
    # The following is for fixing a bug in my design:
    if output_database[0] == "/" and output_database[1] == "/":
        output_database = output_database[1:]
    output_database = pathlib.Path(output_database)
    try:
        output_database.mkdir(parents=True, exist_ok=False)
    except FileExistsError as error:
        print(error)
        print("Files are saved/overwritten in an existing directory.")
    finally:
        return str(output_database) + "/"


def whole_from_segment(
    property_: str,
    segments: List[Tuple[str]],
    parser: Callable,
    geometry: str,
    group: str,
    relation: str,
    save_to: Optional[str] = None
) -> Dict[str, np.ndarray]:
    """generates `whole` array for `property_` of the particle `group` in the
    `geometry` of interest from its `segments`.

    Parameters
    ----------
    property_ : str
        The physical property.
    children : list of tuples
        List of tuples where each tuple at least has one member (the path to
        a csv file for the `property_`).
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    group: {'bug', 'all'}
        Type of the particle group.
    relation : {'histogram', 'tseries', 'bin_edges'}
        Relation between segments and wholes:

        'hisotgram'
            Child is a histogram-like 'segment' file, so the children of
            a parent should be sum horizontally (along "axis=0" in
            numpy's lingo,).

        'tseries'
            Child is a time-series-like 'segment' file, so the children of
            a parent should be concatnated vertically (along "axis=0" in
            Pandas's lingo).

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
    of 'geomtery', and 'group', and the definitons of their keywords.
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(group, ['bug', 'all'])
    invalid_keyword(relation, ['histogram', 'tseries', 'bin_edge'])
    mapping_func = {
        'histogram': lambda whole: (whole[0], np.sum(whole[1], axis=0)),
        'tseries': lambda whole: (whole[0], np.concatenate(whole[1])),
        # 'segments' of a 'whole' have the same bin_edge, so np.unique used
        # to pick one.
        'bin_edge': lambda whole: (whole[0], np.unique(np.array(whole[1])))
    }
    wholes = {}
    for segment in segments:
        segment_info = parser(
            segment[0],
            geometry,
            group,
            'segment'
        )
        whole_name = getattr(segment_info, 'whole')
        child_arr = np.load(segment[0])
        if not bool(wholes):  # is ens_names empty or not?
            wholes[whole_name] = [child_arr]
        elif whole_name not in wholes.keys():
            wholes[whole_name] = [child_arr]
        else:
            wholes[whole_name].append(child_arr)
    wholes = dict(
        map(
            mapping_func[relation],
            wholes.items()
        )
     )
    if save_to is not None:
        _ = dict(
            map(
                lambda whole: (
                    whole[0],
                    save_parent(
                        whole[0], whole[1],
                        property_, save_to,
                        group=group, ext='npy'
                    )
                ),
                wholes.items()
            )
        )
    return wholes


def whole_from_file(
    whole_paths: List[Tuple[str]],
    parser: Callable,
    geometry: str,
    group: str,
) -> Dict[str, np.ndarray]:
    """Loads `whole` numpy arrays for a given physical property of the
    particle `group` in the `geometry` of interest from their pathes
    `wholepathes`.

    Parameters
    ----------
    observations : list of tuples
        List of tuples where each tuple at least has one member (the path to
        a csv file for the `property_`).
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    group: {'bug', 'all'}
        Type of the particle group.

    Return
    ------
    siblings: dict of np.ndarray
        Dict of siblings where keys are 'whole' names (str) and values
        are whole data (arrays).

    Notes
    -----
    Please see the 'organizer' documentation for definitions
    of 'geomtery', and 'group', and the definitons of their keywords.
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(group, ['bug', 'all'])
    wholes = {}
    for whole_path in whole_paths:
        whole_info = parser(
            whole_path[0],
            geometry,
            group,
            'whole'
        )
        whole_name = getattr(whole_info, 'lineage_name')
        wholes[whole_name] = np.load(whole_path[0])
    return wholes


def ensemble_from_vectors(
    ens: Dict[str, np.ndarray]
) -> Tuple[str, pd.DataFrame]:
    """create an "ensemble" dataframe from a dictionary of wholes where
    each "whole" is a numpy vector or 1D array.

    Parameters
    ----------
    ens: dict of np.ndarray
        A dictionary of the wholes' names and arrays in an ensembles.

    Return
    ------
    A tuple of ensemble name and its assocaited dataframe.
    """
    return (ens[0], pd.DataFrame.from_dict(ens[1], orient='columns'))


def ensemble_from_matrix(ens: Dict[str, np.ndarray]) -> Tuple[str, np.ndarray]:
    """create an "ensemble" dataframe from a dictionary of wholes where
    each "whole" is a numpy vector or 1D array.

    Parameters
    ----------
    ens: dict of np.ndarray
        A dictionary of the wholes' names and arrays in an ensembles.

    Return
    ------
    A tuple of ensemble name and its assocaited dataframe.
    """
    return (ens[0], np.stack(list(ens[1].values()), axis=0))


def ensemble(
    property_: str,
    wholes: Dict[str, np.ndarray],
    parser: Callable,
    geometry: str,
    group: str,
    whole_type: str,
    edge_wholes: Optional[Dict[str, np.ndarray]] = None,
    save_to: str = None
) -> Dict[str, Union[pd.DataFrame, np.ndarray]]:
    """Generates ensembles from `wholes` for the physical property `property_`
    of a particle `group` in a `geometry` of interest.

    Parameters
    ----------
    property_ : str
        The physical property.
    wholes : dict of np.ndarray
        A dictionary in which keys are 'whole' names and values are data
        (1D array).
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        The shape of the simulation box.
    group: {'bug', 'all'}
        The type of the particle group.
    whole_type: {'vector', 'matrix'}
        The type of "whole" files.

        'vector':
            A numpy 1D array; for example, a time series or a histogram.

        'matrix':
            A numpy 2D array; for example, the gyration matrix.

    edge_wholes:  dict of np.ndarray, default None
        A dictionary in which keys are 'whole' names and values are bin_edges.
        This option is used if `wholes` are histograms. Since the bin edges (
        and thus bin centers) are the same for all the "whole" histogram in an
        ensemble, it is written just once as a column to the final "ensemble"
        dataframe. `edge_wholes` is meaningful if "vector" `whole_type` is
        used.
    save_to : str, default None
        An/a absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    ensembles : dict of pd.DataFrame
        A dictionary for `property_` in which keys are ensemble names and
        values are dataframes. In each dataframe, the columns are wholes of
        that ensemble.
    """
    # Averging over ensembles with simailar initial paramters
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(group, ['bug', 'all'])
    invalid_keyword(whole_type, ['vector', 'matrix'])
    ensembles = {}
    bin_centers = {}
    for w_name, w_arr in wholes.items():
        w_info = parser(
            w_name,
            geometry,
            group,
            'whole',
            ispath=False
        )
        ens_name = getattr(w_info, 'ensemble_long')
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
        "vector": {
            "mapping_func": ensemble_from_vectors,
            "ext": "csv"
        },
        "matrix": {
            "mapping_func": ensemble_from_matrix,
            "ext": "npy"
        }
    }
    ensembles = dict(
        map(
            whole_types[whole_type]["mapping_func"],
            ensembles.items()
        )
    )
    if edge_wholes is not None:
        for ens_name in ensembles.keys():
            ensembles[ens_name]['bin_center'] = bin_centers[ens_name]
    if save_to is not None:
        _ = dict(
            map(
                lambda ensemble: (ensemble[0],
                                  save_parent(
                                    ensemble[0], ensemble[1],
                                    property_, save_to,
                                    group=group,
                                    ext=whole_types[whole_type]["ext"]
                                    )
                                  ),
                ensembles.items()
            )
        )
    return ensembles


def ens_avg_from_df(
    ens_property: str,
    ens_data: pd.DataFrame,
    exclude: List[str]
) -> Tuple[str, pd.DataFrame]:
    """Creates an "ensAvg" dataframe from a "ensemble" dataframe. The columns
    in the "ensemble" dataframe `ens_data` are the "whole" data and any other
    variable/data given by `exclude`; for instance, if the ensembles is a
    histogram, then there is "bin_center" column in addition to the "whole"
    columns.

    Parameters
    ----------
    ens_property: str
        The property name to which theis "ensemble" data belongs.
    ens_data: pd.DataFrame
        The dataframe of "ensemble"data.
    exclude: list of str
        The list of columns other than "whole" columns.

    Return
    ------
    ens_data: pd.DataFrame
        Update the `ens_data` with avergaing statisicts.
    """
    wholes = [col for col in ens_data.columns if col not in exclude]
    ens_data[ens_property + '-mean'] = ens_data[wholes].mean(axis=1)
    ens_data[ens_property + '-var'] = ens_data[wholes].var(axis=1, ddof=1)
    ens_data[ens_property + '-sem'] = ens_data[wholes].sem(axis=1, ddof=1)
    ens_data.drop(columns=wholes, inplace=True)
    return ens_data


def ens_avg_from_ndarray(
    ens_property: str,
    ens_data: np.ndarray,
    exclude: List[str]
) -> Dict[str, np.ndarray]:
    """Creates an "ensAvg" matrix from a "ensemble" matrix. The first axis of
    the "ensemble" matrix (axis=0 in numpy lingo) contains the whole matrices
    and statsitical measurements are performed on this axis.

    Parameters
    ----------
    ens_property: str
        The property name to which theis "ensemble" data belongs.
    ens_data: pd.DataFrame
        The matrix of "ensemble" data.
    exclude: list of str
        The list of columns other than "whole" columns. It does not have any
        application but is defined for consisitnacy with `ens_avg_from_df`.

    Return
    ------
    ens_data: dict
        A dictionary in which values of various statistical measures and
        values are the ensemble-averages.
    """
    ens_avg = {}
    ens_avg[ens_property + '-mean'] = np.mean(ens_data, axis=0)
    ens_avg[ens_property + '-var'] = np.var(ens_data, axis=0, ddof=1)
    ens_avg[ens_property + '-sem'] = \
        np.std(ens_data, axis=0, ddof=1) / ens_data.shape[0]**0.5
    return ens_avg


def ensemble_avg(
    property_: str,
    ensembles: Dict[str, pd.DataFrame],
    geometry: str,
    group: str,
    ens_type: str,
    exclude: list = ['bin_center'],
    save_to: str = None
) -> Dict[str, pd.DataFrame]:
    """perform averaging over the "whole" data in each "ensemble" data in the
    `ensembles` of the physical property `property_` of a particle `group`, if
    that columns is a valid 'whole' simulation name, not a `exclude` columns.

    The typ

    Parameters
    ----------
    property_ : str
        The physical property.
    siblings : dict of DataFrame
        Dictionary of siblings (dataframes) where keys are parent names and
        values are siblings (dataframes). In each siblings' dataframe, the
        number of columns is equal to or more than the number of
        children of that parent. Columns' names are the children' names and
        the names given by `exclude`.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        The shape of the simulation box.
    group: {'bug', 'all'}
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
    # Averging over ensembles with simailar initial paramters
    invalid_keyword(group, ['bug', 'all'])
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(ens_type, ['dataframe', 'ndarray'])
    ens_avgs = {}
    ens_types = {
        "dataframe": {
            "ens_avg_func": ens_avg_from_df,
            "ext": "csv"
        },
        "ndarray": {
            "ens_avg_func": ens_avg_from_ndarray,
            "ext": "npy"
        }
    }
    for ens, ens_data in ensembles.items():
        ens_property = ens + '-' + property_
        ens_avg = ens_types[ens_type]["ens_avg_func"](
           ens_property,
           ens_data,
           exclude
        )
        ens_avgs[ens] = ens_avg
    if save_to is not None:
        property_ = property_ + '-ensAvg'
        _ = dict(
            map(
                lambda ens: (ens[0],
                             save_parent(
                                 ens[0],
                                 ens[1],
                                 property_, save_to,
                                 group=group,
                                 ext=ens_types[ens_type]["ext"]
                             )
                             ),
                ens_avgs.items()
            )
        )
    return ens_avgs


def children_stamps(
    stamps: List[Tuple[str]],
    group: str,
    lineage: str,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """generates a dataset of the phyiscal attributes and equilibrium
    (time-averaged) physical properties of all the `lineage` simulations of a
    particle `group` in a 'space' in a `geometry`.

    The name of space is created from the first value of 'filename' in
    the generated dataset.

    Parameters
    ----------
    stamps: list of tuple
        List of tuples where each tumple has one member and that member is a
        filepath to the stamp of a 'segment' simulation in a space.
    group: {'bug', 'all'}
        Type of the particle group.
    lineage: {'segment', 'whole'}
        Lineage type of children' stamps
    save_to : str, default None
        Absolute or relative path of a directory to which outputs are saved.

    Return
    ------
    space_stamps: pd.DataFrame
        Dataframe of all the children stamps in the `group` in a space.
    """
    invalid_keyword(group, ['bug', 'all'])
    invalid_keyword(lineage, ['segment', 'whole'])
    space_stamps = [pd.read_csv(stamp[0]) for stamp in stamps]
    space_stamps = pd.concat(space_stamps)
    space_stamps.reset_index(inplace=True, drop=True)
    if lineage == 'whole':
        # Some older version of parsers use 'segment'= "N/A" and
        # 'segment_id'="N/A" in a "whole" stamp when "whole" linage
        # is used.
        try:
            cols_to_drop = ['segment', 'segment_id']
            space_stamps.drop(columns=cols_to_drop, inplace=True)
            warnings.warn(
                "'segment' and 'segment_id' columns are dropped when"
                " individual 'whole' stamps combined to create a single"
                " dataframe of 'whole' stamps by 'children_stamps'.",
                UserWarning
            )
        except KeyError:
            print(f"'{cols_to_drop}' are not among columns.")
    if save_to is not None:
        space_name = space_stamps.loc[0, 'space']
        filename = '-'.join([space_name, group, lineage, 'stamps.csv'])
        space_stamps.to_csv(save_to + filename, index=False)
    return space_stamps


def parents_stamps(
    stamps: pd.DataFrame,
    geometry: str,
    group: str,
    lineage: str,
    properties: Optional[Dict[str, Callable]] = None,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """perform merging/ensemble-averaging over all the 'segment/'whole'
    simulation stamps in a 'space' in a given `geometry` for a given
    `group` basedon the given `lineage`.

    Parameters
    ----------
    stamps: DataFrame
        Dataframe of all the simulation stamps in the `group` in a space.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box
    group: {'bug', 'all'}
        Type of the particle group.
    lineage: {('segment', 'whole'}
        Lineage type of children's stamps.
    properties: dict of str
        A dictionary in which the keys are properties such as the time-
        averaged radius of gyration which are measured during the 'probe'
        phase and the values are user-defined or numpy functions which are
        used as the aggregation function bu pandas.
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
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    invalid_keyword(group, ['bug', 'all'])
    invalid_keyword(lineage, ['segment', 'whole'])
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
            "a previous call of 'parents_stamps' function."
        )
    # aggregation dictionary: See Note above.
    agg_funcs = dict()
    attr_agg_funcs = ['last'] * len(stamps_cols)
    agg_funcs.update(zip(stamps_cols, attr_agg_funcs))
    if properties is not None:  # add/update agg funcs for properties.
        agg_funcs.update(properties)
    # Handing 'lineage'-specific details:
    if lineage == 'whole':
        parent_groupby = 'ensemble_long'
        # aggregating functions for properties
        agg_funcs['ensemble_id'] = 'count'
        agg_funcs['n_frames'] = 'last'
        file_lastname = 'ensAvg'
    else:
        parent_groupby = 'whole'
        # aggregating functions for properties
        agg_funcs['segment_id'] = 'count'
        agg_funcs['n_frames'] = 'sum'
        file_lastname = 'whole'
    agg_funcs.pop(parent_groupby)
    parents_stamps = stamps.groupby([parent_groupby]).agg(agg_funcs)
    parents_stamps.reset_index(inplace=True)
    if lineage == 'whole':
        parents_stamps.rename(
            columns={'ensemble_id': 'n_ensembles'},
            inplace=True
        )
        # If the 'whole' stamps are generated directly in the 'probe' phase,
        # then 'segment' and 'segment_id' columns are "N/A" and are removed
        # from the list of stamps columns that are added to the parents
        # stamps.
        # There is no need to have the "n_segment" column in the parent
        # stamps, so it is removed. The "whole" stamps directly generated
        # in the "probe" phase do not have such a column, but those generated
        # from "segment" stamps have.
        # Droping redundant columns silently:
        parents_stamps.drop(
            columns=['n_segments', 'segment_id', 'segment'],
            inplace=True,
            errors='ignore'
            )
    else:
        parents_stamps.rename(
            columns={'segment_id': 'n_segments'},
            inplace=True
        )
    if save_to is not None:
        space_name = parents_stamps.loc[0, 'space']
        filename = '-'.join(
            [space_name, group, 'stamps', file_lastname + '.csv']
        )
        parents_stamps.to_csv(save_to + filename, index=False)
    return parents_stamps


def unique_property(
    filepath: str,
    prop_idx: int,
    extensions: List[str],
    drop_properties: Optional[List[str]] = None,
    sep: Optional[str] = "-"
) -> Tuple[List[str], List[str]]:
    """Finds unique physical properties and physical property-measures by
    spliting filenames given by the glob-friendly 'filepath'.

    A measure refers to some measurement done on a physical property.

    A physcial property-measure is defined as a measurement done on a physical;
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
        The proeprties that should be ignored.
    sep: str, default "-"
        The seperator between a "property" and its "measure".

    Return
    ------
    uniq_props: list of str
        A sorted liste of unique physcial properties.
    uniq_prop_measures: list of str
        A sorted list of unique property-measures.
    """
    props_measures = glob(filepath)
    uniq_prop_measures = []
    for ext in extensions:
        prop_measure_per_ext = list(
            set(
                [sep.join(
                    prop.split("/")[-1].split(ext)[0].split(sep)[prop_idx:]
                    ) for prop in props_measures]
                )
            )
        uniq_prop_measures.extend(prop_measure_per_ext)
    for drop_property in drop_properties:
        try:
            uniq_prop_measures.remove(drop_property)
        except ValueError:
            print(
                f"'{drop_property}' is not among unique properties."
                )
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
    parser: Callable,
    hierarchy: str,
    physical_attrs: List[str],
    species: str,
    group: str,
    geometry: str,
    is_save: Optional[bool] = False
) -> pd.DataFrame:
    """Takes the `property_path` to 'ensAvg' time series of a given `property_`
    in a given space `input_database`,  adds the `physical_attrs` of interest
    as the new columns to each 'ensAvg' dataframe, and merges all the 'ensAvg'
    dataframes into one 'space' dataframe along the 0 (or 'row' or 'index')
    in pandas's lingo,

    In each 'ensemble-averaged' dataframe, there are 3 columns with this name
    pattern:

    column name = '[long_ensemble]-[porperty_][-measure]-[stat]'

    where '[-measure]' is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. the 'stat' keyword is either 'mean',
    'ver', or 'sem'.

    Parameters
    ----------
    property_path: str
        Path to the the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    property_pattern: str
        The pattern by which the filenames of timeseries are started with; for
        instance, "N*" means files start with "N"
    attributes: list of str
        The physical attributes that will added as new columns to the
        concatenated timeseries.
    species: {'Mon', 'Crd', 'Foci'}
        The species of particles.
    group: {'bug', 'all'}
        The type of the particle group.
    geometry : {'biaxial', 'slit', 'box'}
        The shape of the simulation box.
    single_space:
        Whether the all-in-one file is for all the timeseries properties of a
        single space or all the space in a project.
    is_save : bool, default False
        whether to save output to file or not.

    Return
    ------
    all_in_one: pandas.DataFrame
        a dataframe in which all the timeseries are concatenated along `orient`
        of interest, and "properties and attributes" of interest are added to
        it as the new columns.
    """
    invalid_keyword(species, ['Mon', 'Crd', 'Foci'])
    property_ext = "-" + property_ + "-ensAvg.csv"
    ens_avg_csvs = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs, fmts=[property_ext])
    property_db = []
    for ens_avg_csv in ens_avg_csvs:
        ens_avg_df = pd.read_csv(ens_avg_csv[0], header=0)
        # the first column of porperty_df is used to extract
        # the information about the property and the space it
        # belongs to.
        # columns in ens_avg_df are ensemble long names
        ens_long_name = ens_avg_df.columns[0].split('-')[0]
        property_info = parser(
            ens_long_name,
            geometry,
            group,
            'ensemble_long',
            ispath=False
        )
        # Column names wihtout 'ensemble_long' name
        # See the explanantion in doc abbut column names.
        col_names = {
            col: "-".join(col.split('-')[1:]) for col in ens_avg_df.columns
        }
        ens_avg_df.rename(columns=col_names, inplace=True)
        ens_avg_df.reset_index(inplace=True)
        ens_avg_df.rename(columns={'index': 'time'}, inplace=True)
        ens_avg_df['time'] = ens_avg_df['time'] * property_info.dt
        for attr_name in physical_attrs:
            ens_avg_df[attr_name] = getattr(property_info, attr_name)
        property_db.append(ens_avg_df)
    property_db = pd.concat(property_db, axis=0)
    property_db.reset_index(inplace=True, drop=True)
    if is_save is not False:
        save_to_space = database_path(
            input_database,
            'analysis',
            stage='space',
            group='bug'
            )
        space = save_to_space.split("/")[-2].split("-")[0]
        output = "-".join([space, group, property_]) + "-space.csv"
        property_db.to_csv(save_to_space + output, index=False)
    return property_db


def space_hists(
    input_database: str,
    property_: str,
    parser: Callable,
    hierarchy: str,
    physical_attrs: List[str],
    species: str,
    group: str,
    geometry: str,
    bin_center: Optional[np.ndarray] = None,
    normalize: Optional[bool] = False,
    divisor: Optional[float] = 0.025,
    round_to: Optional[int] = 3,
    is_save: Optional[bool] = False
) -> pd.DataFrame:
    """Takes the `property_path` to 'ensAvg' time series of a given `property_`
    in a given space `input_database`,  adds the `physical_attrs` of interest
    as the new columns to each 'ensAvg' dataframe, and merges all the 'ensAvg'
    dataframes into one 'space' dataframe along the 0 (or 'row' or 'index')
    in pandas's lingo,

    In each 'ensemble-averaged' dataframe, there are 4 columns with this name
    pattern:

    column name = '[long_ensemble]-[porperty_][-measure]-[stat]'

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
    property_path: str
        Path to the the timeseries of the physical property of interest.
    property_: str
        Name of the physical property of interest.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    property_pattern: str
        The pattern by which the filenames of timeseries are started with; for
        instance, "N*" means files start with "N"
    attributes: list of str
        The physical attributes that will added as new columns to the
        concatenated timeseries.
    species: {'Mon', 'Crd', 'Foci'}
        The species of particles.
    group: {'bug', 'all'}
        The type of the particle group.
    geometry : {'biaxial', 'slit', 'box'}
        The shape of the simulation box.
    single_space:
        Whether the all-in-one file is for all the timeseries properties of a
        single space or all the space in a project.
    bin_center: numpy array, default None
        The bin centers. The argument should be given if the 'bin_center' is
        not in the ensemble-averaged dataframe and the same array of bin
        centers is used in all different ensemble-averaged dataframes. This is
        the case for "clustersHistTFoci" or "bondsHistTFoci" properties, but
        not for "zHistMon" or "rHistCrd".
    normalize: bool, default False
        Whether normalize hists or not.
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

    Requirenents
    ------------
    PolyPhys, Pandas
    """
    invalid_keyword(species, ['Mon', 'Crd', 'Foci'])
    property_ext = "-" + property_ + "-ensAvg.csv"
    ens_avg_csvs = glob(input_database + hierarchy + property_ext)
    ens_avg_csvs = sort_filenames(ens_avg_csvs, fmts=[property_ext])
    property_db = []
    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        # the first column of porperty_df is used to extract
        # the information about the property and the space it
        # belongs to.
        # columns in ens_avg_df are ensemble long names
        ens_long_name = ens_avg.columns[0].split('-')[0]
        property_info = parser(
            ens_long_name,
            geometry,
            group,
            'ensemble_long',
            ispath=False
        )
        # Column names wihtout 'ensemble_long' name
        # See the explanantion in doc about column names.
        col_names = [col for col in ens_avg.columns if col != 'bin_center']
        col_names = {col: "-".join(col.split('-')[1:]) for col in col_names}
        col_names['bin_center'] = property_ + '-bin_center'
        if bin_center is not None:
            ens_avg['bin_center'] = bin_center.tolist()
        ens_avg.rename(columns=col_names, inplace=True)
        if normalize is True:
            ens_avg[property_+'-norm'] = \
                ens_avg[property_+'-mean'] / ens_avg[property_+'-mean'].sum()
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)
        ens_avg['phi_c_bulk_round'] = ens_avg['phi_c_bulk'].apply(
            round_up_nearest, args=[divisor, round_to])
        property_db.append(ens_avg)
    property_db = pd.concat(property_db, axis=0)
    property_db.reset_index(inplace=True, drop=True)
    if is_save is not False:
        save_to_space = database_path(
            input_database,
            'analysis',
            stage='space',
            group='bug'
            )
        space = save_to_space.split("/")[-2].split("-")[0]
        output = "-".join([space, group, property_]) + "-space.csv"
        property_db.to_csv(save_to_space + output, index=False)
    return property_db
