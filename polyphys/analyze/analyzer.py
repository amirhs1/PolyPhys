from typing import (
    Callable,
    List,
    NewType,
    Tuple,
    Optional,
    Union
)
from glob import glob
from polyphys.manage.organizer import (
    invalid_keyword,
    sort_filenames,
    database_path,
    whole_from_file,
    whole_from_segment,
    ensemble,
    ensemble_avg,
    children_stamps,
    parents_stamps
)

from polyphys.analyze.distributions import distributions_generator
from polyphys.analyze.correlations import acf_generator
import numpy as np
import pandas as pd

PropertyT = NewType('PropertyT', str)
SpeciesT = NewType('SpeciesT', str)
GroupT = NewType('GroupT', str)
DirectionT = NewType('DirectionT', str)
AxisT = NewType('AxisT', int)

TimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT]
NonScalarTimeSeriesT = Tuple[PropertyT, SpeciesT, GroupT, AxisT]
HistogramT = Tuple[DirectionT, SpeciesT, GroupT]


def time_series(
    observations: List[str],
    parser: Callable,
    geometry: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    tseries_properties: Optional[List[TimeSeriesT]] = None,
    acf_tseries_properties: Optional[List[TimeSeriesT]] = None,
    nlags: int = 7000,
    alpha: float = 0.05
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `TimeSeriesT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `is_segment` is `True`, `observations` are "segements" and
    are "vertically" merged to create "wholes".

    Issue
    -----
    If the order of the acf_generator and the first ensemble_avg is
    change, the column names in each ensemble dataframe in 'ensembles'
    changes.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    is_segment: bool
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-aveages are saved.
    tseries_properties: list of HistogramT
        A list of tumple where each tumple has three members: the property
        name, species, and group of a 'time-series' property.
    acf_tseries_properties: list of HistogramT
        A list of tumple where each tumple has three members: the property
        name, species, and group of a 'time-series' property. For
        `cls_tseries_properties`, the the auto correlation function (AFC) is
        also computed.
    nlags: int, default 7000
        Maximum lag in the auto correlation function (AFC).
    alpha: float, default 0.05
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett”s formula.
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    if tseries_properties is not None:
        for property_, species, group in tseries_properties:
            tseries = sort_filenames(
                observations,
                fmts=['-' + property_ + species + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    relation='tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group
                )
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
            )
    if acf_tseries_properties is not None:
        for property_, species, group in acf_tseries_properties:
            tseries = sort_filenames(
                observations,
                fmts=['-' + property_ + species + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    'tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group
                )
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens
                )
            acfs, lower_cls, upper_cls = acf_generator(
                property_ + species,
                ensembles,
                nlags,
                alpha,
                group,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acf',
                    acfs,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acfLowerCi',
                    lower_cls,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
                )
            _ = ensemble_avg(
                    property_ + species + '-acfUpperCi',
                    upper_cls,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
                )


def histograms(
    observations: List[str],
    parser: Callable,
    geometry: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    hist_properties: Optional[List[HistogramT]] = None,
    rho_phi_hist_properties: Optional[List[HistogramT]] = None
) -> None:
    """Runs various statistical analyses on `observations` of
    each of `HistogramT` types in a given `geometry` and then
    writes the ensembles and ensemble-averages of time series and its
    associated analyses to file.

    If the `segment` is `True`, `observations` are "segements" and
    are "horizontally" merged to create "wholes".

    Issue
    -----
    HThis function only work for spatial distributions created by
    SpatialDistribution class.

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}, default biaxial
        Shape of the simulation box.
    is_segment: bool, default False
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-aveages are saved.
    rho_phi_hist_properties: list of HistogramT, default None
        A list of tumple where each tumple has four members: the direction,
        direction long name, species, and group of a 'histogram' property.
        These histogram properties are then used to calculate the local
        number density and volume fraction.
    hist_properties: list of HistogramT default None
        A list of tumple where each tumple has three members: the direction,
        species, and group of a 'histogram' property.
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    # Histograms:
    # Two types of histograms with and without rhos and phis:
    # rho: local number density
    # phi: locla volume fraction
    if rho_phi_hist_properties is not None:
        for direction, species, group in rho_phi_hist_properties:
            hists = sort_filenames(
                observations,
                fmts=['-' + direction + 'Hist' + species + '.npy']
            )
            edges = sort_filenames(
                observations,
                fmts=['-' + direction + 'Edge' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    'histogram',
                    save_to=save_to_whole
                )
                edge_wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    edges,
                    parser,
                    geometry,
                    group,
                    'histogram',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group
                )
                edge_wholes = whole_from_file(
                    edges,
                    parser,
                    geometry,
                    group
                )
            # 'whole' dataframes, each with a 'whole' columns.
            rho_wholes, phi_wholes = distributions_generator(
                wholes,
                edge_wholes,
                group,
                species,
                geometry,
                direction,
                parser,
                save_to=save_to_whole
            )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                parser,
                geometry,
                group,
                save_to=save_to_ens_avg
            )
            ensembles = ensemble(
                direction + 'Rho' + species,
                rho_wholes,
                parser,
                geometry,
                group,
                edge_wholes=edge_wholes,
                save_to=None
            )
            _ = ensemble_avg(
                direction + 'Rho' + species,
                ensembles,
                parser,
                geometry,
                group,
                save_to=save_to_ens_avg
            )
            ensembles = ensemble(
                direction + 'Phi' + species,
                phi_wholes,
                parser,
                geometry,
                group,
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Phi' + species,
                ensembles,
                parser,
                geometry,
                group,
                save_to=save_to_ens_avg
            )
        del rho_wholes, phi_wholes, ensembles
    if hist_properties is not None:
        for species, direction, group in hist_properties:
            hists = sort_filenames(
                observations,
                fmts=['-' + direction + 'Hist' + species + '.npy']
            )
            edges = sort_filenames(
                observations,
                fmts=['-' + direction + 'Edge' + species + '.npy']
            )
            if is_segment is True:
                wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    hists,
                    parser,
                    geometry,
                    group,
                    'histogram',
                    save_to=save_to_whole
                )
                edge_wholes = whole_from_segment(
                    direction + 'Hist' + species,
                    edges,
                    parser,
                    geometry,
                    group,
                    'histogram',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    hists,
                    parser,
                    geometry,
                    group
                )
                edge_wholes = whole_from_file(
                    edges,
                    parser,
                    geometry,
                    group
                )
            ensembles = ensemble(
                direction + 'Hist' + species,
                wholes,
                parser,
                geometry,
                group,
                edge_wholes=edge_wholes,
                save_to=save_to_ens
            )
            _ = ensemble_avg(
                direction + 'Hist' + species,
                ensembles,
                parser,
                geometry,
                group,
                save_to=save_to_ens_avg
            )


def nonscalar_time_series(
    observations: List[str],
    parser: Callable,
    geometry: str,
    is_segment: bool,
    save_to: Tuple[Union[str, None], Union[str, None], Union[str, None]],
    nonscalar_hist_properties: Optional[List[NonScalarTimeSeriesT]] = None,
    nonscalar_matrix_properties: Optional[List[NonScalarTimeSeriesT]] = None,
) -> None:
    """Runs overs all 'segment' `observations` of each of
    "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
    a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
    and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
    associated analyses to file. `nonscalar_hist_properties` are time-varying
    histograms while `nonscalar_matrix_properties` are time-varying matrixes.

    A non_scalar property is either a time-varing 1D array (a vector) or a
    time-varing 2D one (a matrix). See the "Notes" and "Issues" below.

    Notes
    -----
    Currently, the vector-like properties are histograms collected over time so
    thery are passed by `nonscalar_hist_properties` argument.

    If the `is_segment` is `True`, `observations` are "segements" and
    are "vertically" merged to create "wholes".

    Issues
    ------
    1. Based on the Notes above, different type of nonscalar properties should
    be passed to this function and new functions should be defined for
    matrix-like nonscalar properties, so it can be possible to define the
    ensemble and ensemble-averaged versions.

    2. `nonscalar_hist_properties` is currently list of taime-varing histogram
    that do not need any more process as those passed to "histograms".

    Parameters
    ----------
    observations: list of str
        List of path to different observations generated by 'probed'
        trajectory files.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    is_segment: bool
        Whether `observations` are 'segment' or 'whole'
    save_to : tuple of three str
        Absolute or relative path of the directories to which wholes,
        ensembles, and ensemble-aveages are saved.
    nonscalar_hist_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has four string members. The
        first string is the name of a physical property, the second one is
        the particletype, the third one is `group` type, and the last one
        is the axis over which the. These physical
        properties are all of nonscalar form.
    nonscalar_matrix_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    """
    save_to_whole, save_to_ens, save_to_ens_avg = save_to
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    if nonscalar_hist_properties is not None:
        for property_, species, group, avg_axis in nonscalar_hist_properties:
            tseries = sort_filenames(
                    observations,
                    fmts=['-' + property_ + species + '.npy']
                )
            if is_segment is True:
                wholes = whole_from_segment(
                    property_ + species,
                    tseries,
                    parser,
                    geometry,
                    group,
                    'tseries',
                    save_to=save_to_whole
                )
            else:
                wholes = whole_from_file(
                    tseries,
                    parser,
                    geometry,
                    group
                )
            wholes = {whole_name: np.mean(whole_array, axis=avg_axis)
                      for whole_name, whole_array in wholes.items()
                      }
            ensembles = ensemble(
                    property_ + species,
                    wholes,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens
            )
            _ = ensemble_avg(
                    property_ + species,
                    ensembles,
                    parser,
                    geometry,
                    group,
                    save_to=save_to_ens_avg
            )
    if nonscalar_matrix_properties is not None:
        raise NotImplementedError("This part of function is not defined yet!")


def analyze_bug(
    input_database: str,
    hierarchy: str,
    parser: Callable,
    geometry: str,
    is_segment: bool, 
    tseries_properties: Optional[List[TimeSeriesT]] = None,
    acf_tseries_properties: Optional[List[TimeSeriesT]] = None,
    hist_properties: Optional[List[HistogramT]] = None,
    rho_phi_hist_properties: Optional[List[HistogramT]] = None,
    nonscalar_hist_properties: Optional[List[NonScalarTimeSeriesT]] = None,
    nonscalar_matrix_properties: Optional[List[NonScalarTimeSeriesT]] = None,
    nlags: int = 7000,
    alpha: float = 0.05
) -> None:
    """read in the 'probe' observations of the 'group' particles based on the
    `hierarchy` of directories and files from the `input_database` path to the
    'probe' phase of a 'space' and creates the 'analysis' phase at that parent
    directory of the 'probe' of that 'space', infers 'space' names from
    `input_database` path and creates a 'space' directories at various stages
    in the 'analysis' directory for both 'bug' and 'all' groups.

    `tseries_properties`, `hists_properties`, `rho_phi_hists_properties` are
    list of tuples in which each tuple has three string members. The first
    string is the name of a physical property, the second one is the particle
    type, and the last one is `group` type.

    Parameters
    ----------
    input_database: str
        Path to the input_database; a 'space' directory at a given 'phase'.
    hierarchy: str
        Hierarchy of the directories and files within the `input_database`;
        for instance, "/N*/N*" means files that starts with "N" and are
        located in directories starting with "N".
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    is_segment: bool
        Whether `observations` are 'segment' (True) or 'whole' (False)
    nonscalar_properties: list of TimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    tseries_properties: list of TimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of the time-series form.
    acf_tseries_properties: list of TimeSeriesT, default None
        A list of tumple where each tumple has three members: the property
        name, species, and group of a 'time-series' property. For
        `cls_tseries_properties`, the the auto correlation function (AFC) is
        also computed.
    hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of the histogram form.
    rho_phi_hist_properties: list of HistogramT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of the histogram form; however, in contrast to
        `hists_properties`, the local number denisty and volume fraction of
        `rho_phi_hists_properties` are also calculated.
    nonscalar_hist_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has four string members. The
        first string is the name of a physical property, the second one is
        the particletype, the third one is `group` type, and the last one
        is the axis over which the. These physical
        properties are all of nonscalar form.
    nonscalar_matrix_properties: list of NonScalarTimeSeriesT, default None
        A list of tuples in which each tuple has three string members. The
        first string is the name of a physical property, the second one is
        the particletype, and the last one is `group` type. These physical
        properties are all of nonscalar form.
    nlags: int, default 7000
        Maximum lag in the auto correlation function (AFC).
    alpha: float, default 0.05
        If a number is given, the confidence intervals for the given level
        are returned. For instance if alpha=.05, 95 % confidence intervals
        are returned where the standard deviation is computed according to
        Bartlett”s formula.
    """
    invalid_keyword(geometry, ['biaxial', 'slit', 'box'])
    observations = glob(input_database + hierarchy)
    if observations == []:
        raise ValueError(
            "File not found in "
            f"'{input_database + hierarchy}'"
            )
    # 'bug' save_to paths:
    save_to_ens = database_path(
        input_database, 'analysis', stage='ens', group='bug'
    )
    save_to_ens_avg = database_path(
        input_database, 'analysis', stage='ensAvg', group='bug'
    )
    # stamps:
    stamp_files = sort_filenames(observations, fmts=['-stamps.csv'])

    if is_segment is True:
        save_to_whole = database_path(
            input_database, 'analysis', stage='wholeSim', group='bug'
        )
        segments_stamps = children_stamps(
            stamp_files,
            'bug',
            'segment',  # lineage of the children stamps
            save_to=save_to_whole  # save all the segment stamps in one file
        )
        whole_stamps = parents_stamps(
            segments_stamps,
            geometry,
            'bug',
            'segment',  # lineage of the children stamps
            save_to=save_to_ens  # save all the whole stamps
        )
    else:
        save_to_whole = None
        whole_stamps = children_stamps(
            stamp_files,
            geometry,
            'bug',
            'whole',  # lineage of the children stamps
            save_to=save_to_ens  # save all the whole stamps
        )
    _ = parents_stamps(
        whole_stamps,
        geometry,
        'bug',
        'whole',  # lineage of the children stamps
        save_to=save_to_ens_avg  # save all the ensemble-averaged stamps
    )
    # physical properties
    if tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            tseries_properties=tseries_properties,
        )
    if acf_tseries_properties is not None:
        time_series(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            acf_tseries_properties=acf_tseries_properties,
            nlags=nlags,
            alpha=alpha
        )
    if hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            hist_properties=hist_properties,
        )
    if rho_phi_hist_properties is not None:
        histograms(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            rho_phi_hist_properties=rho_phi_hist_properties,
        )
    if nonscalar_hist_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_hist_properties=nonscalar_hist_properties,
        )
    if nonscalar_matrix_properties is not None:
        nonscalar_time_series(
            observations,
            parser,
            geometry,
            is_segment,
            (save_to_whole, save_to_ens, save_to_ens_avg),
            nonscalar_matrix_properties=nonscalar_matrix_properties,
        )


def error_calc_block(
    data: np.ndarray,
    save_to: Optional[str] = None
) -> pd.DataFrame:
    """
    computes the statistical inefficiency (si) and uncertainty associated with
     a physical quantity calculated in a simulation.

    Using Flybbjerg-Peterson block average method, the plateau should be
     evident after 6-8 transformations.

    Parameters
    ----------
    data: np.array
        Inpute data
    filename: str
        Name of output file
    to_file: bool
        Whether save results to file or not.

    Return:
    block_analysis: pd.dataframe
        A pandas dataframe in which there are sevweral columns of data about
         block-averaging error analysis.

    References
    ----------
    Original python code snippets are from "Computer simulation of liquids",
    Allen MP Tildesley DJ (2017):
    https://github.com/Allen-Tildesley/examples/blob/master/python_examples/error_calc.py

    https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py

    "Error estimates on averages of correlated data",HG Flyvbjerg and
    H Petersen. J Chem Phys, 91, 461 (1989): https://doi.org/10.1063/1.457480
    """
    nframes = len(data)  # size of the data
    var = np.var(data, ddof=1)  # Bias-corrected sample variance
    var_err = np.sqrt(2 / (nframes-1)) * var  # error in bias-corrected var
    sem = np.sqrt(var / nframes)  # correlations are neglected
    sem_err = np.sqrt(1 / (2*(nframes-1))) * sem  # error in SEM
    blocks = data.copy()
    ntransfroms = np.zeros(0, dtype=np.int8)  # number of block tranformartion
    ntransfroms = np.append(ntransfroms, 0)
    nblocks = np.zeros(0, dtype=np.int8)  # number of blocks
    nblocks = np.append(nblocks, nframes)
    bsize = np.zeros(0, dtype=np.int8)  # size of each block
    bsize = np.append(bsize, 1)
    bvar = np.zeros(0)  # block variances
    bvar = np.append(bvar, var)
    bvar_err = np.zeros(0)  # error in block variances
    bvar_err = np.append(bvar_err, var_err)
    bsem = np.zeros(0)  # block sems
    bsem = np.append(bsem, sem)
    bsem_err = np.zeros(0)  # error in sem
    bsem_err = np.append(bsem_err, sem_err)
    si = np.zeros(0)  # statistical inefficiency (si)
    si_initial = bsize[-1] * bvar[-1] / var  # inintial si
    si = np.append(si, si_initial)
    si_err = np.zeros(0)  # error in si
    si_initial_err = np.sqrt(2 / (nframes-1)) * si_initial  # initial si error
    si_err = np.append(si_err, si_initial_err)

    while True:  # Loop over number, and hence length, of blocks
        # halving nblocks and rounding if it is odd:
        if nblocks[-1] <= 3:  # loop counter
            break
        nblocks = np.append(nblocks, nblocks[-1] // 2)
        bsize = np.append(bsize, bsize[-1] * 2)
        blocks[0:nblocks[-1]] = (
            blocks[0:2*nblocks[-1]-1:2] + blocks[1:2*nblocks[-1]:2]
            ) / 2.0  # Blocking transformation, halving the data set
        bvar = np.append(bvar, np.var(blocks[0:nblocks[-1]], ddof=1))
        bvar_err = np.append(
            bvar_err, np.sqrt(2 / (nblocks[-1]-1)) * bvar[-1]
        )
        bsem = np.append(bsem, np.sqrt(bvar[-1] / nblocks[-1]))
        bsem_err = np.append(
            bsem_err, np.sqrt((1 / (2 * (nblocks[-1]-1)))) * bsem[-1]
            )
        si = np.append(si, bsize[-1] * bvar[-1] / var)
        si_err = np.append(
            si_err, np.sqrt((1 / (2 * (nblocks[-1]-1)))) * si[-1]
            )
        ntransfroms = np.append(ntransfroms, ntransfroms[-1]+1)

    cols = [
        "ntransfroms", "bsize", "nblocks", "var", "var_err", "sem", "sem_err",
        "si", "si_err"
        ]
    block_analysis = pd.DataFrame(
        data=np.stack(
            (ntransfroms, bsize, nblocks, bvar, bvar_err, bsem, bsem_err, si,
                si_err),
            axis=1
        ),
        columns=cols
    )
    if save_to is not None:
        block_analysis.to_csv(save_to + '-block_average.csv', index=False)
    return block_analysis
