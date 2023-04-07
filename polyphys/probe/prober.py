from typing import Optional, Dict, Any, Union
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mda_dist
from MDAnalysis import transformations as mda_trans
from MDAnalysis.analysis.base import AnalysisFromFunction
import numpy as np

from polyphys.manage.parser import (
    SumRuleCyl, TransFociCyl, TransFociCub, HnsCub, HnsCyl
    )
import polyphys.manage.typer as ptyper
from polyphys.manage.organizer import invalid_keyword
from polyphys.analyze import clusters, correlations
from polyphys.analyze.measurer import transverse_size, fsd, end_to_end
import warnings


def stamps_report_with_measures(
    report_name: str,
    sim_info: Union[TransFociCub, TransFociCyl],
    n_frames: int,
    measures: Dict[str, float],
) -> None:
    """
    Writes a summary of stamps (properties and attributes) of a simulation
    to file.

    `stamps_report` generates a dataset called `report_name` of the
    values of some attributes of `sim_info`, the number of frames
    `n_frames`, all the key and value pairs in all the given dictionaries
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: ParserT
        A ParserT instant object that contains information about the name,
        parents,and physical attributes of a simulation.
    n_frames: int
        Number of frames/snapshots/configurations in a simulation.
    measures: Dict
        A dictionary of measures where a key and value pair is the name and
        value of a physical property.
    """
    with open(report_name, mode='w') as report:
        # write header
        for lineage_name in sim_info._genealogy:
            report.write(f"{lineage_name},")
        for attr_name in sim_info.attributes:
            report.write(f"{attr_name},")
        for measure_name in measures.keys():  # each measure is a dict
            report.write(f"{measure_name},")
        report.write("n_frames\n")
        # write values
        for lineage_name in sim_info.genealogy:
            attr_value = getattr(sim_info, lineage_name)
            report.write(f"{attr_value},")
        for attr_name in sim_info.attributes:
            attr_value = getattr(sim_info, attr_name)
            report.write(f"{attr_value},")
        for measure_value in measures.values():
            report.write(f"{measure_value},")
        report.write(f"{n_frames}")


def stamps_report(
    report_name: str,
    sim_info: ptyper.ParserT,
    n_frames: int
) -> None:
    """
    Writes a summary of stamps (properties and attributes) of a simulation
    to file.

    `stamps_report` generates a dataset called `report_name` of the
    values of some attributes of `sim_info`, the number of frames
    `n_frames`, all the key and value pairs in all the given dictionaries
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: ParserT
        A SumRule object that contains information about the name, parents,
        and physical attributes of a simulation.
    n_frames: int
        Number of frames/snapshots/configurations in a simulation.
    """
    with open(report_name, mode='w') as report:
        # write header
        for lineage_name in sim_info.genealogy:
            report.write(f"{lineage_name},")
        for attr_name in sim_info.attributes:
            report.write(f"{attr_name},")
        report.write("n_frames\n")
        # write values
        for lineage_name in sim_info.genealogy:
            attr_value = getattr(sim_info, lineage_name)
            report.write(f"{attr_value},")
        for attr_name in sim_info.attributes:
            attr_value = getattr(sim_info, attr_name)
            report.write(f"{attr_value},")
        report.write(f"{n_frames}")


def bin_create(
    sim_name: str,
    edge_name: str,
    bin_size: float,
    lmin: float,
    lmax: float,
    save_to: str
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generates arrays of bins and histograms

    Parameters
    ----------
    sim_name: str
        Name of the simulation.
    edge_name: str
        Name of the variable for which the histogram is computed.
    bin_size : float
        Size of each bin.
    lmin : float
        Lower bound of the system in the direction of interest.
    lmax : float
        Upper bound of the system in the direction of interest.
    save_to : str
        Whether save outputs to memory as csv files or not.

    Return
    ------
    bin_edges : numpy array of float
        The edges to pass into a histogram. Save `bin_edges` to file if
        `save_to` is not None.
    hist: array of int
        An empty histogram
    """
    bin_edges = np.arange(lmin, lmax + bin_size, bin_size)
    hist = np.zeros(len(bin_edges) - 1, dtype=np.int16)
    np.save(save_to + sim_name + '-' + edge_name + '.npy', bin_edges)
    return bin_edges, hist


def fixedsize_bins(
    sim_name: str,
    edge_name: str,
    bin_size: float,
    lmin: float,
    lmax: float,
    bin_type: str = 'ordinary',
    save_to: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Generates arrays of bins and histograms, ensuring that the `bin_size`
    guaranteed. To achieve this, it extends the `lmin` and `lmax` limits.

    To-do List
    ----------
    1. Following the idea used:
    https://docs.mdanalysis.org/1.1.1/_modules/MDAnalysis/lib/util.html#fixedwidth_bins
    Makes input array-like so bins can be calculated for 1D data (then all
    parameters are simple floats) or nD data (then parameters are supplied
    as arrays, with each entry corresponding to one dimension).
    2. Eliminate the if-statement for the periodic_bin_edges.

    Parameters
    ----------
    sim_name: str
        Name of the simulation.
    edge_name: str
        Name of the variable for which the histogram is computed.
    bin_size : float
        Size of each bin.
    lmin : float
        Lower bound of the system in the direction of interest.
    lmax : float
        Upper bound of the system in the direction of interest.
    bin_type: {'ordinary', 'nonnegative', 'periodic'}, default 'ordinary'
        The type of bin in a given direction in a given coordinate system:

        'ordinary'
            A bounded or unbounded coordinate such as any of the cartesian
            coordinates or the polar coordinate in the spherical coordinate
            system. For such coordinates, the `lmin` and `lmax` limits are
            equally extended to ensure the `bin_size`.

        'nonnegative'
            A nonnegative coordinate such as the r direction in the polar
            or spherical coordinate system. For such coordinates, ONLY `lmax`
            limit is extended to ensure the `bin_size`. `lmin` is either 0.0
            or a positive number smaller than `lmax`.

        'periodic'
            A periodic coordinate such as the azimuthal direction in the
            spherical coordinate. It is assumed that 'period'=`lmax`-`lmin`;
            therefore, if 'period' is not a multiple of `bin_size`, then an
            array of bin_edges is used; otherwise, n_bins is used.

    save_to : str, default None
        Whether save outputs to memory as npy files or not.

    Return
    ------
    bin_edges : numpy array of  float
        The edges to pass into a histogram. Save `bin_edges` to file if
        `save_to` is not None.
    hist: array of  int
        An empty histogram

    Reference:
    https://docs.mdanalysis.org/1.1.1/documentation_pages/lib/util.html#MDAnalysis.analysis.density.fixedwidth_bins
    """
    hist_collectors = 0
    bin_edges = 0
    bin_types = ['ordinary', 'nonnagative', 'periodic']
    if not np.all(lmin < lmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = bin_size
    _lmin = lmin
    _lmax = lmax
    _length = _lmax - _lmin
    n_bins: int = 0
    if bin_type == 'ordinary':
        n_bins = int(np.ceil(_length / _delta))
        dl = 0.5 * (n_bins * _delta - _length)  # excess length
        # add half of the excess to each end:
        _lmin = _lmin - dl
        _lmax = _lmax + dl
        # create empty grid with the right dimensions (and get the edges)
        hist_collectors, bin_edges = np.histogram(
            np.zeros(1),
            bins=n_bins,
            range=(_lmin, _lmax)
        )
    elif bin_type == 'nonnegative':
        n_bins = int(np.ceil(_length / _delta))
        dl = 0.5 * (n_bins * _delta - _length)
        _lmin = _lmin - dl
        _lmax = _lmax + dl
        if _lmin <= 0.0:
            _lmin = 0.0
            # add full of the excess to upper end:
            _lmax = _lmax + 2 * dl
        hist_collectors, bin_edges = np.histogram(
            np.zeros(1),
            bins=n_bins,
            range=(_lmin, _lmax)
        )
    elif bin_type == 'periodic':  # Assuming that the _length=period:
        n_bins = int(np.ceil(_length / _delta))
        print(
            f"Number of bins '{n_bins}'"
            " is more than or equal to the actual number of bins in "
            f"'periodic' bin type since the 'period=lmax-min={_length}'"
            f"and '{_delta}'"
            "(not 'n_bins') are used to created 'bin_edges'."
            )
        bin_edges = np.arange(_lmin, _lmax + _delta, _delta)
        hist_collectors, bin_edges = np.histogram(
            np.zeros(1),
            bins=bin_edges,
            range=(_lmin, _lmax)
        )
    else:
        invalid_keyword(bin_type, bin_types)
    hist_collectors = hist_collectors * 0
    hist_collectors_std = hist_collectors * 0
    if save_to is not None:
        np.save(save_to + sim_name + '-' + edge_name + '.npy', bin_edges)
    results = {
        'n_bins': n_bins,
        'bin_edges': bin_edges,
        'collector': hist_collectors,
        'collector_std': hist_collectors_std,
        'range': (_lmin, _lmax)
    }
    return results


def frame_dist_hist(
    distances: np.ndarray, hist_info: Dict[str, Any]
) -> np.ndarray:
    """
    Calculate the histogram of `distances`.

    Parameters
    ----------
    distances: np.ndarray
        The distances of a group of atoms form a given origin in a given frame.
    hist_info: Dict[str, any]
        The information about the histogram.

    Return
    ------
    frame_hist: np.ndarray
        The counts of atoms in each bin.
    """
    frame_hist, _ = np.histogram(
        distances,
        bins=hist_info['bin_edges'],
        range=hist_info['range']
    )
    return frame_hist


def write_hists(
    hist_infos: Dict[str, Any],
    sim_name: str,
    save_to: str,
    std: bool = False
) -> None:
    """
    Writes histogram per species per direction to file.

    Parameters
    ----------
    hist_infos : Dict[str, Any]
        A dict of dicts that contains the information about direction, species,
         and histograms.
    sim_name: str
        The name of simulation file to which the `hist_infos` belongs.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    std : bool, default False
        _description_, by default False
    """
    for dir_ in hist_infos.keys():
        for species, hist in hist_infos[dir_].items():
            np.save(
                save_to + sim_name + '-' + dir_ + species + '.npy',
                hist['collector']
            )
            if std is True:
                np.save(
                    save_to + sim_name + '-' + dir_ + 'Std' + species + '.npy',
                    hist['collector_std']
                )
        # end of loop
    # end of loop


def sum_rule_bug_cyl(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = './',
    continuous: bool = False
) -> None:
    """
    Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        print(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.")
    print("Setting the name of analyze file...")
    sim_info = SumRuleCyl(
        trajectory,
        lineage,
        'cylindrical',
        'bug',
        'linear',
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(sim_info.mmon * sim_info.eps_others) \
        # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt)
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups
    bug = cell.select_atoms('resid 1')  # the polymer
    # collectors
    trans_size_t = []
    fsd_t = []
    rflory_t = []
    gyr_t = []
    asphericity_t = []
    shape_parameter_t = []
    principal_axes_t = np.empty([0, 3, 3])
    # dict of bin edges:
    # bin_edges = {
    #    'rfloryEdge': {
    #        'bin_size': 0.5 * sim_info.dmon,
    #        'lmin': 0,
    #        'lmax': sim_info.nmon * sim_info.dmon
    #    }
    # }
    # distribution of the size of the end-to-end
    # rflory_hist_info: Dict[str, Any] = fixedsize_bins(
    #    sim_name,
    #    'rfloryEdgeMon',
    #    bin_edges['rfloryEdge']['bin_size'],
    #    bin_edges['rfloryEdge']['lmin'],
    #    bin_edges['rfloryEdge']['lmax'],
    #    bin_type='nonnegative',
    #    save_to=save_to
    # )
    # if any([
    #        rflory_hist_info['collector'].any() != 0,
    #        rflory_hist_info['collector_std'].any() != 0,
    #        ]):
    #    raise ValueError("One of the histogram collectors are not empty!")
    for _ in sliced_trj:
        # various measures of chain size
        trans_size_t.append(transverse_size(bug))
        fsd_t.append(fsd(bug.positions))
        gyr_t.append(bug.radius_of_gyration())
        rms = end_to_end(bug.positions)
        rflory_t.append(rms)
        # frame_hist, _ = np.histogram(
        #    rms,
        #    bins=rflory_hist_info['bin_edges'],
        #    range=rflory_hist_info['range']
        # )
        # # RDF of the end-to-end distance
        # rflory_hist_info['collector'] += frame_hist
        # rflory_hist_info['collector_std'] += np.square(frame_hist)
        # shape parameters:
        asphericity_t.append(bug.asphericity(wrap=False, unwrap=False))
        shape_parameter_t.append(bug.shape_parameter(wrap=False))
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(wrap=False)]),
            axis=0
        )

    np.save(save_to + sim_name + '-transSizeTMon.npy', np.array(trans_size_t))
    np.save(save_to + sim_name + '-fsdTMon.npy', np.array(fsd_t))
    np.save(save_to + sim_name + '-gyrTMon.npy', np.array(gyr_t))
    np.save(save_to + sim_name + '-rfloryTMon.npy', np.array(rflory_t))
    # np.save(
    #    save_to + sim_name + '-rfloryHistMon.npy',
    #    rflory_hist_info['collector']
    # )
    # np.save(
    #    save_to + sim_name + '-rfloryHistStdMon.npy',
    #    rflory_hist_info['collector_std']
    # )
    np.save(save_to + sim_name + '-asphericityTMon.npy', np.array(
        asphericity_t))
    np.save(save_to + sim_name + '-shapeTMon.npy', np.array(shape_parameter_t))
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def sum_rule_all_cyl(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """
    Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest, and saves a variety of
    outputs (mostly in the csv format) to the `save_to` directory.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        print(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.")
    print("Setting the name of analyze file...\n")
    sim_info = SumRuleCyl(
        trajectory,
        lineage,
        'cylindrical',
        'all',
        'linear'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {
            'bin_size':  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            'lmin': 0,
            'lmax': 0.5 * sim_info.dcyl
            },
        'zEdge': {
            'bin_size':  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            'lmin': -0.5 * sim_info.lcyl,
            'lmax': 0.5 * sim_info.lcyl
            },
        'thetaEdge': {
            'bin_size':  np.pi / 36,
            'lmin': -1 * np.pi,
            'lmax': np.pi
            }
        }
    # LJ time difference between two consecutive frames:
    time_unit = \
        sim_info.dmon * np.sqrt(sim_info.mmon * sim_info.eps_others)
    lj_nstep = sim_info.adump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cell = mda.Universe(
        topology,
        trajectory,
        topology_format='DATA',
        format='LAMMPSDUMP',
        lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z",
        dt=sim_real_dt
    )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds = cell.select_atoms('resid 0')  # crowders
    bug = cell.select_atoms('resid 1')  # the chain or monomers
    # bin edges and histograms in different directions:
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        'xEdge',
        bin_edges['xEdge']['bin_size'],
        bin_edges['xEdge']['lmin'],
        bin_edges['xEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        'yEdge',
        bin_edges['yEdge']['bin_size'],
        bin_edges['yEdge']['lmin'],
        bin_edges['yEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        'zEdge',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            x_hist_info['collector'].any() != 0,
            x_hist_info['collector_std'].any() != 0,
            y_hist_info['collector'].any() != 0,
            y_hist_info['collector_std'].any() != 0,
            z_hist_info['collector'].any() != 0,
            z_hist_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    # z direction of the cylindrical coordinate system
    z_hist_crd_info = fixedsize_bins(
        sim_name,
        'zEdgeCrd',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    z_hist_mon_info = fixedsize_bins(
        sim_name,
        'zEdgeMon',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # theta of the cylindrical coordinate system
    theta_hist_crd_info = fixedsize_bins(
        sim_name,
        'thetaEdgeCrd',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
    )  # in radians
    theta_hist_mon_info = fixedsize_bins(
        sim_name,
        'thetaEdgeMon',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
    )  # in radians
    # check if any of the histograms are empty or not.
    if any([
            r_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector'].any() != 0,
            z_hist_mon_info['collector'].any() != 0,
            z_hist_crd_info['collector'].any() != 0,
            theta_hist_mon_info['collector'].any() != 0,
            theta_hist_crd_info['collector'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            z_hist_mon_info['collector_std'].any() != 0,
            z_hist_crd_info['collector_std'].any() != 0,
            theta_hist_mon_info['collector_std'].any() != 0,
            theta_hist_crd_info['collector_std'].any() != 0
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # 2D hists
    # crd
    # # xy
    hist_crd_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_crd_info_xy['collector'] = np.zeros(hist_crd_info_xy['n_bins'])
    hist_crd_info_xy['collector'] *= 0
    # # xz
    hist_crd_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_xz['collector'] = np.zeros(hist_crd_info_xz['n_bins'])
    hist_crd_info_xz['collector'] *= 0
    # # yz
    hist_crd_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_yz['collector'] = np.zeros(hist_crd_info_yz['n_bins'])
    hist_crd_info_yz['collector'] *= 0
    # mon
    # # xy
    hist_mon_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_mon_info_xy['collector'] = np.zeros(hist_mon_info_xy['n_bins'])
    hist_mon_info_xy['collector'] *= 0
    # # xz
    hist_mon_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_xz['collector'] = np.zeros(hist_mon_info_xz['n_bins'])
    hist_mon_info_xz['collector'] *= 0
    # # yz
    hist_mon_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_yz['collector'] = np.zeros(hist_mon_info_yz['n_bins'])
    hist_mon_info_yz['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        pos_r = np.linalg.norm(crds.positions[:, :2], axis=1)
        frame_hist, _ = np.histogram(
            pos_r,
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
        )
        r_hist_crd_info['collector'] += frame_hist
        r_hist_crd_info['collector_std'] += np.square(frame_hist)
        # bug
        pos_r = np.linalg.norm(bug.positions[:, :2], axis=1)
        frame_hist, _ = np.histogram(
            pos_r,
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
        )
        r_hist_mon_info['collector'] += frame_hist
        r_hist_mon_info['collector_std'] += np.square(frame_hist)
        # histogram in z direction
        # crds
        pos_z = crds.positions[:, 2]
        frame_hist, _ = np.histogram(
            pos_z,
            bins=z_hist_crd_info['bin_edges'],
            range=z_hist_crd_info['range']
        )
        z_hist_crd_info['collector'] += frame_hist
        z_hist_crd_info['collector_std'] += np.square(frame_hist)
        # bug
        pos_z = bug.positions[:, 2]
        frame_hist, _ = np.histogram(
            pos_z,
            bins=z_hist_mon_info['bin_edges'],
            range=z_hist_mon_info['range']
        )
        z_hist_mon_info['collector'] += frame_hist
        z_hist_mon_info['collector_std'] += np.square(frame_hist)
        # histogram in theta
        # crds
        theta = np.arctan2(
            crds.positions[:, 1],
            crds.positions[:, 0]
        )  # in radians between [-np.pi, np.pi]
        frame_hist, _ = np.histogram(
            theta,
            bins=theta_hist_crd_info['bin_edges'],
            range=theta_hist_crd_info['range']
        )
        theta_hist_crd_info['collector'] += frame_hist
        theta_hist_crd_info['collector_std'] += np.square(frame_hist)
        # bug
        theta = np.arctan2(
            bug.positions[:, 1],
            bug.positions[:, 0]
        )  # in radian
        frame_hist, _ = np.histogram(
            theta,
            bins=theta_hist_mon_info['bin_edges'],
            range=theta_hist_mon_info['range']
        )
        theta_hist_mon_info['collector'] += frame_hist
        theta_hist_mon_info['collector_std'] += np.square(frame_hist)
        # histogram in 2 dimensional space
        # crds
        # # xy
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=hist_crd_info_xy['bin_edges'],
            range=hist_crd_info_xy['range'],
        )
        hist_crd_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=hist_crd_info_xz['bin_edges'],
            range=hist_crd_info_xz['range'],
        )
        hist_crd_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=hist_crd_info_yz['bin_edges'],
            range=hist_crd_info_yz['range'],
        )
        hist_crd_info_yz['collector'] += frame_hist
        # bug
        # # xy
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=hist_mon_info_xy['bin_edges'],
            range=hist_mon_info_xy['range'],
        )
        hist_mon_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=hist_mon_info_xz['bin_edges'],
            range=hist_mon_info_xz['range'],
        )
        hist_mon_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=hist_mon_info_yz['bin_edges'],
            range=hist_mon_info_yz['range'],
        )
        hist_mon_info_yz['collector'] += frame_hist
    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info
        },
        'zHist': {
            'Crd': z_hist_crd_info,
            'Mon': z_hist_mon_info
        },
        'thetaHist': {
            'Crd': theta_hist_crd_info,
            'Mon': theta_hist_mon_info
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': hist_crd_info_xy,
            'Mon': hist_mon_info_xy
        },
        'xzHist': {
            'Crd': hist_crd_info_xz,
            'Mon': hist_mon_info_xz
        },
        'yzHist': {
            'Crd': hist_crd_info_yz,
            'Mon': hist_mon_info_yz
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def trans_foci_bug_cyl(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = TransFociCyl(
        trajectory,
        lineage,
        'cylindrical',
        'bug',
        'ring'
        )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cluster_cutoff = sim_info.dmon_large + sim_info.dcrowd
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups
    foci = cell.select_atoms('type 2')  # the foci
    bug = cell.select_atoms('resid 1')  # the bug/polymer: small and large mons
    # defining collectors
    # -bug:
    trans_size_t = []
    fsd_t = []
    rflory_t = []
    gyr_t = []
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = []
    shape_parameter_t = []
    # -foci:
    foci_t = np.empty([0, sim_info.nmon_large, sim_info.nmon_large])
    dir_contacts_t = np.empty(
        [0, sim_info.nmon_large, sim_info.nmon_large], dtype=int
    )
    bonds_t = np.empty([0, sim_info.nmon_large], dtype=int)
    clusters_t = np.empty([0, sim_info.nmon_large+1], dtype=int)
    for _ in sliced_trj:
        # bug:
        # various measures of chain size
        trans_size_t.append(transverse_size(bug))
        fsd_t.append(fsd(bug.positions))
        gyr_t.append(bug.radius_of_gyration())
        rflory_t.append(end_to_end(bug.positions))
        # shape parameters:
        asphericity_t.append(bug.asphericity(wrap=False, unwrap=False))
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(wrap=False)]),
            axis=0
        )
        shape_parameter_t.append(bug.shape_parameter(wrap=False))
        # foci:
        dist_sq_mat = clusters.dist_sq_matrix(foci.positions)
        foci_pair_dist = np.triu(dist_sq_mat, 1)
        foci_pair_dist = np.sqrt(foci_pair_dist)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t = np.append(foci_t, np.array([foci_pair_dist]), axis=0)
        dir_contacts = clusters.direct_contact(dist_sq_mat, cluster_cutoff)
        dir_contacts_t = np.append(
            dir_contacts_t,
            np.array([dir_contacts]),
            axis=0
        )
        bonds_stat = clusters.count_bonds(dir_contacts)
        bonds_t = np.append(bonds_t, np.array([bonds_stat]), axis=0)
        contacts = clusters.contact_matrix(dir_contacts)
        clusters_stat = clusters.count_clusters(contacts)
        clusters_t = np.append(clusters_t, np.array([clusters_stat]), axis=0)
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-transSizeTMon.npy', np.array(trans_size_t))
    np.save(save_to + sim_name + '-fsdTMon.npy', np.array(fsd_t))
    np.save(save_to + sim_name + '-rfloryTMon.npy', np.array(rflory_t))
    np.save(save_to + sim_name + '-gyrTMon.npy', np.array(gyr_t))
    np.save(save_to + sim_name + '-asphericityTMon.npy',
            np.array(asphericity_t)
            )
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', np.array(shape_parameter_t))
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def trans_foci_all_cyl(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest, and saves a variety of
    outputs (mostly in the csv format) to the `save_to` directory.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = TransFociCyl(
        trajectory,
        lineage,
        'cylindrical',
        'all',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {
            'bin_size':  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            'lmin': 0,
            'lmax': 0.5 * sim_info.dcyl
            },
        'zEdge': {
            'bin_size':  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            'lmin': -0.5 * sim_info.lcyl,
            'lmax': 0.5 * sim_info.lcyl
            },
        'thetaEdge': {
            'bin_size':  np.pi / 36,
            'lmin': -1 * np.pi,
            'lmax': np.pi
            }
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    lj_nstep = sim_info.adump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cell = mda.Universe(
        topology, trajectory,
        topology_format='DATA',
        format='LAMMPSDUMP',
        lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z",
        dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    dna = cell.select_atoms('type 1')  # small monomers
    foci = cell.select_atoms('type 2')  # large monomers
    bug = cell.select_atoms('resid 1')  # the polymer
    crds = cell.select_atoms('resid 0')  # crowders
    # bin edges and histograms in different directions:
    # radial direction of the cylindrical coordinate system
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    # z direction of the cylindrical coordinate system
    z_hist_crd_info = fixedsize_bins(
        sim_name,
        'zEdgeCrd',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    z_hist_foci_info = fixedsize_bins(
        sim_name,
        'zEdgeFoci',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    z_hist_dna_info = fixedsize_bins(
        sim_name,
        'zEdgeDna',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    z_hist_mon_info = fixedsize_bins(
        sim_name,
        'zEdgeMon',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # theta of the cylindrical coordinate system
    theta_hist_crd_info = fixedsize_bins(
        sim_name,
        'thetaEdgeCrd',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
        )  # in radians
    theta_hist_foci_info = fixedsize_bins(
        sim_name,
        'thetaEdgeFoci',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
        )  # in radians
    theta_hist_dna_info = fixedsize_bins(
        sim_name,
        'thetaEdgeDna',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
        )  # in radians
    theta_hist_mon_info = fixedsize_bins(
        sim_name,
        'thetaEdgeMon',
        bin_edges['thetaEdge']['bin_size'],
        bin_edges['thetaEdge']['lmin'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
        )  # in radians
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd_info['collector'].any() != 0,
            r_hist_foci_info['collector'].any() != 0,
            r_hist_dna_info['collector'].any() != 0,
            r_hist_mon_info['collector'].any() != 0,
            z_hist_crd_info['collector'].any() != 0,
            z_hist_foci_info['collector'].any() != 0,
            z_hist_dna_info['collector'].any() != 0,
            z_hist_mon_info['collector'].any() != 0,
            theta_hist_crd_info['collector'].any() != 0,
            theta_hist_foci_info['collector'].any() != 0,
            theta_hist_dna_info['collector'].any() != 0,
            theta_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            r_hist_foci_info['collector_std'].any() != 0,
            r_hist_dna_info['collector_std'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
            z_hist_crd_info['collector_std'].any() != 0,
            z_hist_foci_info['collector_std'].any() != 0,
            z_hist_dna_info['collector_std'].any() != 0,
            z_hist_mon_info['collector_std'].any() != 0,
            theta_hist_crd_info['collector_std'].any() != 0,
            theta_hist_foci_info['collector_std'].any() != 0,
            theta_hist_dna_info['collector_std'].any() != 0,
            theta_hist_mon_info['collector_std'].any() != 0
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # bin edges and histograms in different directions:
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        'xEdge',
        bin_edges['xEdge']['bin_size'],
        bin_edges['xEdge']['lmin'],
        bin_edges['xEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        'yEdge',
        bin_edges['yEdge']['bin_size'],
        bin_edges['yEdge']['lmin'],
        bin_edges['yEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        'zEdge',
        bin_edges['zEdge']['bin_size'],
        bin_edges['zEdge']['lmin'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            x_hist_info['collector'].any() != 0,
            x_hist_info['collector_std'].any() != 0,
            y_hist_info['collector'].any() != 0,
            y_hist_info['collector_std'].any() != 0,
            z_hist_info['collector'].any() != 0,
            z_hist_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # 2D hists
    # crd
    # # xy
    hist_crd_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_crd_info_xy['collector'] = np.zeros(hist_crd_info_xy['n_bins'])
    hist_crd_info_xy['collector'] *= 0
    # # xz
    hist_crd_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_xz['collector'] = np.zeros(hist_crd_info_xz['n_bins'])
    hist_crd_info_xz['collector'] *= 0
    # # yz
    hist_crd_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_yz['collector'] = np.zeros(hist_crd_info_yz['n_bins'])
    hist_crd_info_yz['collector'] *= 0
    # mon
    # # xy
    hist_mon_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_mon_info_xy['collector'] = np.zeros(hist_mon_info_xy['n_bins'])
    hist_mon_info_xy['collector'] *= 0
    # # xz
    hist_mon_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_xz['collector'] = np.zeros(hist_mon_info_xz['n_bins'])
    hist_mon_info_xz['collector'] *= 0
    # # yz
    hist_mon_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_yz['collector'] = np.zeros(hist_mon_info_yz['n_bins'])
    hist_mon_info_yz['collector'] *= 0
    # foci
    # # xy
    hist_foci_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_foci_info_xy['collector'] = np.zeros(hist_foci_info_xy['n_bins'])
    hist_foci_info_xy['collector'] *= 0
    # # xz
    hist_foci_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_foci_info_xz['collector'] = np.zeros(hist_foci_info_xz['n_bins'])
    hist_foci_info_xz['collector'] *= 0
    # # yz
    hist_foci_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_foci_info_yz['collector'] = np.zeros(hist_foci_info_yz['n_bins'])
    hist_foci_info_yz['collector'] *= 0
    # dna
    # # xy
    hist_dna_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_dna_info_xy['collector'] = np.zeros(hist_dna_info_xy['n_bins'])
    hist_dna_info_xy['collector'] *= 0
    # # xz
    hist_dna_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_dna_info_xz['collector'] = np.zeros(hist_dna_info_xz['n_bins'])
    hist_dna_info_xz['collector'] *= 0
    # # yz
    hist_dna_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_dna_info_yz['collector'] = np.zeros(hist_dna_info_yz['n_bins'])
    hist_dna_info_yz['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        pos_r = np.linalg.norm(crds.positions[:, :2], axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_crd_info)
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # foci
        pos_r = np.linalg.norm(foci.positions[:, :2], axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_foci_info)
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # dna
        pos_r = np.linalg.norm(dna.positions[:, :2], axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_dna_info)
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # bug
        pos_r = np.linalg.norm(bug.positions[:, :2], axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_mon_info)
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # histogram in z direction
        # crds
        pos_z = crds.positions[:, 2]
        pos_hist = frame_dist_hist(pos_z, z_hist_crd_info)
        z_hist_crd_info['collector'] += pos_hist
        z_hist_crd_info['collector_std'] += np.square(pos_hist)
        # foci
        pos_z = foci.positions[:, 2]
        pos_hist = frame_dist_hist(pos_z, z_hist_foci_info)
        z_hist_foci_info['collector'] += pos_hist
        z_hist_foci_info['collector_std'] += np.square(pos_hist)
        # dna
        pos_z = dna.positions[:, 2]
        pos_hist = frame_dist_hist(pos_z, z_hist_dna_info)
        z_hist_dna_info['collector'] += pos_hist
        z_hist_dna_info['collector_std'] += np.square(pos_hist)
        # bug
        pos_z = bug.positions[:, 2]
        pos_hist = frame_dist_hist(pos_z, z_hist_mon_info)
        z_hist_mon_info['collector'] += pos_hist
        z_hist_mon_info['collector_std'] += np.square(pos_hist)
        # histogram in theta
        # crds
        theta = np.arctan2(
            crds.positions[:, 1],
            crds.positions[:, 0]
        )  # in radians between [-np.pi, np.pi]
        pos_hist = frame_dist_hist(theta, theta_hist_crd_info)
        theta_hist_crd_info['collector'] += pos_hist
        theta_hist_crd_info['collector_std'] += np.square(pos_hist)
        # foci
        theta = np.arctan2(
            foci.positions[:, 1],
            foci.positions[:, 0]
        )  # in radians between [-np.pi, np.pi]
        pos_hist = frame_dist_hist(theta, theta_hist_foci_info)
        theta_hist_foci_info['collector'] += pos_hist
        theta_hist_foci_info['collector_std'] += np.square(pos_hist)
        # dna
        theta = np.arctan2(
            dna.positions[:, 1],
            dna.positions[:, 0]
        )  # in radians between [-np.pi, np.pi]
        pos_hist = frame_dist_hist(theta, theta_hist_dna_info)
        theta_hist_dna_info['collector'] += pos_hist
        theta_hist_dna_info['collector_std'] += np.square(pos_hist)
        # bug
        theta = np.arctan2(
            bug.positions[:, 1],
            bug.positions[:, 0]
            )  # in radians between [-np.pi, np.pi]
        pos_hist = frame_dist_hist(theta, theta_hist_mon_info)
        theta_hist_mon_info['collector'] += pos_hist
        theta_hist_mon_info['collector_std'] += np.square(pos_hist)
        # histogram in 2 dimensional space
        # crds
        # # xy
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=hist_crd_info_xy['bin_edges'],
            range=hist_crd_info_xy['range'],
        )
        hist_crd_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=hist_crd_info_xz['bin_edges'],
            range=hist_crd_info_xz['range'],
        )
        hist_crd_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=hist_crd_info_yz['bin_edges'],
            range=hist_crd_info_yz['range'],
        )
        hist_crd_info_yz['collector'] += frame_hist
        # bug
        # # xy
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=hist_mon_info_xy['bin_edges'],
            range=hist_mon_info_xy['range'],
        )
        hist_mon_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=hist_mon_info_xz['bin_edges'],
            range=hist_mon_info_xz['range'],
        )
        hist_mon_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=hist_mon_info_yz['bin_edges'],
            range=hist_mon_info_yz['range'],
        )
        hist_mon_info_yz['collector'] += frame_hist
        # foci
        # # xy
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=hist_foci_info_xy['bin_edges'],
            range=hist_foci_info_xy['range'],
        )
        hist_foci_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=hist_foci_info_xz['bin_edges'],
            range=hist_foci_info_xz['range'],
        )
        hist_foci_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=hist_foci_info_yz['bin_edges'],
            range=hist_foci_info_yz['range'],
        )
        hist_foci_info_yz['collector'] += frame_hist
        # dna
        # # xy
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=hist_dna_info_xy['bin_edges'],
            range=hist_dna_info_xy['range'],
        )
        hist_dna_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=hist_dna_info_xz['bin_edges'],
            range=hist_dna_info_xz['range'],
        )
        hist_dna_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=hist_dna_info_yz['bin_edges'],
            range=hist_dna_info_yz['range'],
        )
        hist_dna_info_yz['collector'] += frame_hist
    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info,
            'Foci': r_hist_foci_info,
            'Dna': r_hist_dna_info
        },
        'zHist': {
            'Crd': z_hist_crd_info,
            'Mon': z_hist_mon_info,
            'Foci': z_hist_foci_info,
            'Dna': z_hist_dna_info
        },
        'thetaHist': {
            'Crd': theta_hist_crd_info,
            'Mon': theta_hist_mon_info,
            'Foci': theta_hist_foci_info,
            'Dna': theta_hist_dna_info
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': hist_crd_info_xy,
            'Mon': hist_mon_info_xy,
            'Foci': hist_foci_info_xy,
            'Dna': hist_dna_info_xy
        },
        'xzHist': {
            'Crd': hist_crd_info_xz,
            'Mon': hist_mon_info_xz,
            'Foci': hist_foci_info_xz,
            'Dna': hist_dna_info_xz
        },
        'yzHist': {
            'Crd': hist_crd_info_yz,
            'Mon': hist_mon_info_yz,
            'Foci': hist_foci_info_yz,
            'Dna': hist_dna_info_yz
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def trans_foci_bug_cub(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = TransFociCub(
        trajectory,
        lineage,
        'cubic',
        'bug',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cluster_cutoff = sim_info.dmon_large + sim_info.dcrowd
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups
    foci = cell.select_atoms('type 2')  # the foci: large monomers
    bug = cell.select_atoms('resid 1')  # the bug/polymer: small and large mons
    # defining collectors
    # -bug:
    gyr_t = []
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = []
    shape_parameter_t = []
    # -foci:
    foci_t = np.empty([0, sim_info.nmon_large, sim_info.nmon_large])
    dir_contacts_t = np.empty(
        [0, sim_info.nmon_large, sim_info.nmon_large], dtype=int
    )
    bonds_t = np.empty([0, sim_info.nmon_large], dtype=int)
    clusters_t = np.empty([0, sim_info.nmon_large+1], dtype=int)
    for _ in sliced_trj:
        # bug:
        # various measures of chain size
        gyr_t.append(bug.radius_of_gyration())
        # shape parameters:
        asphericity_t.append(bug.asphericity(wrap=False, unwrap=True))
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(wrap=False)]),
            axis=0
        )
        shape_parameter_t.append(bug.shape_parameter(wrap=False))
        # foci:
        dist_sq_mat = clusters.dist_sq_matrix(foci.positions)
        foci_pair_dist = np.triu(dist_sq_mat, 1)
        foci_pair_dist = np.sqrt(foci_pair_dist)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t = np.append(foci_t, np.array([foci_pair_dist]), axis=0)
        dir_contacts = clusters.direct_contact(dist_sq_mat, cluster_cutoff)
        dir_contacts_t = np.append(
            dir_contacts_t,
            np.array([dir_contacts]),
            axis=0
        )
        bonds_stat = clusters.count_bonds(dir_contacts)
        bonds_t = np.append(bonds_t, np.array([bonds_stat]), axis=0)
        contacts = clusters.contact_matrix(dir_contacts)
        clusters_stat = clusters.count_clusters(contacts)
        clusters_t = np.append(clusters_t, np.array([clusters_stat]), axis=0)
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-gyrTMon.npy', np.array(gyr_t))
    np.save(save_to + sim_name + '-asphericityTMon.npy',
            np.array(asphericity_t)
            )
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', np.array(shape_parameter_t))
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def trans_foci_all_cub(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest, and saves a variety of
    outputs (mostly in the csv format) to the `save_to` directory.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a
    trajectory or topology file; moreover, LAMMPS recenter is used to
    restrict the center of mass of "bug" (monomers) to the center of
    simulation box; and consequently, coordinates of all the particles in a
    trajectory or topology file is recentered to fulfill this constraint.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = TransFociCub(
        trajectory,
        lineage,
        'cubic',
        'all',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {  # edges for distance r from the box center
            'bin_size':  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            'lmin': 0,
            'lmax': 0.5 * sim_info.lcube
            },
        'lEdge': {  # edges in cartesian coordinates
            'bin_size':  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            'lmin': -0.5 * sim_info.lcube,
            'lmax': 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    lj_nstep = sim_info.adump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
    cell = mda.Universe(
        topology, trajectory,
        topology_format='DATA',
        format='LAMMPSDUMP',
        lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z",
        dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    dna = cell.select_atoms('type 1')  # small monomers
    foci = cell.select_atoms('type 2')  # large monomers
    bug = cell.select_atoms('resid 1')  # the polymer
    crds = cell.select_atoms('resid 0')  # crowders
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd_info['collector'].any() != 0,
            r_hist_foci_info['collector'].any() != 0,
            r_hist_dna_info['collector'].any() != 0,
            r_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            r_hist_foci_info['collector_std'].any() != 0,
            r_hist_dna_info['collector_std'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # bin edges and histograms in different directions:
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        'xEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        'yEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        'zEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            x_hist_info['collector'].any() != 0,
            x_hist_info['collector_std'].any() != 0,
            y_hist_info['collector'].any() != 0,
            y_hist_info['collector_std'].any() != 0,
            z_hist_info['collector'].any() != 0,
            z_hist_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # 2D hists
    # crd
    # # xy
    hist_crd_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_crd_info_xy['collector'] = np.zeros(hist_crd_info_xy['n_bins'])
    hist_crd_info_xy['collector'] *= 0
    # # xz
    hist_crd_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_xz['collector'] = np.zeros(hist_crd_info_xz['n_bins'])
    hist_crd_info_xz['collector'] *= 0
    # # yz
    hist_crd_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_yz['collector'] = np.zeros(hist_crd_info_yz['n_bins'])
    hist_crd_info_yz['collector'] *= 0
    # mon
    # # xy
    hist_mon_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_mon_info_xy['collector'] = np.zeros(hist_mon_info_xy['n_bins'])
    hist_mon_info_xy['collector'] *= 0
    # # xz
    hist_mon_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_xz['collector'] = np.zeros(hist_mon_info_xz['n_bins'])
    hist_mon_info_xz['collector'] *= 0
    # # yz
    hist_mon_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_yz['collector'] = np.zeros(hist_mon_info_yz['n_bins'])
    hist_mon_info_yz['collector'] *= 0
    # foci
    # # xy
    hist_foci_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_foci_info_xy['collector'] = np.zeros(hist_foci_info_xy['n_bins'])
    hist_foci_info_xy['collector'] *= 0
    # # xz
    hist_foci_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_foci_info_xz['collector'] = np.zeros(hist_foci_info_xz['n_bins'])
    hist_foci_info_xz['collector'] *= 0
    # # yz
    hist_foci_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_foci_info_yz['collector'] = np.zeros(hist_foci_info_yz['n_bins'])
    hist_foci_info_yz['collector'] *= 0
    # dna
    # # xy
    hist_dna_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_dna_info_xy['collector'] = np.zeros(hist_dna_info_xy['n_bins'])
    hist_dna_info_xy['collector'] *= 0
    # # xz
    hist_dna_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_dna_info_xz['collector'] = np.zeros(hist_dna_info_xz['n_bins'])
    hist_dna_info_xz['collector'] *= 0
    # # yz
    hist_dna_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_dna_info_yz['collector'] = np.zeros(hist_dna_info_yz['n_bins'])
    hist_dna_info_yz['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        pos_r = np.linalg.norm(crds.positions, axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_crd_info)
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # foci
        pos_r = np.linalg.norm(foci.positions, axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_foci_info)
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # dna
        pos_r = np.linalg.norm(dna.positions, axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_dna_info)
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # bug
        pos_r = np.linalg.norm(bug.positions, axis=1)
        pos_hist = frame_dist_hist(pos_r, r_hist_mon_info)
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # histogram in 2 dimensional space
        # crds
        # # xy
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=hist_crd_info_xy['bin_edges'],
            range=hist_crd_info_xy['range'],
        )
        hist_crd_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=hist_crd_info_xz['bin_edges'],
            range=hist_crd_info_xz['range'],
        )
        hist_crd_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=hist_crd_info_yz['bin_edges'],
            range=hist_crd_info_yz['range'],
        )
        hist_crd_info_yz['collector'] += frame_hist
        # bug
        # # xy
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=hist_mon_info_xy['bin_edges'],
            range=hist_mon_info_xy['range'],
        )
        hist_mon_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=hist_mon_info_xz['bin_edges'],
            range=hist_mon_info_xz['range'],
        )
        hist_mon_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=hist_mon_info_yz['bin_edges'],
            range=hist_mon_info_yz['range'],
        )
        hist_mon_info_yz['collector'] += frame_hist
        # foci
        # # xy
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=hist_foci_info_xy['bin_edges'],
            range=hist_foci_info_xy['range'],
        )
        hist_foci_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=hist_foci_info_xz['bin_edges'],
            range=hist_foci_info_xz['range'],
        )
        hist_foci_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=hist_foci_info_yz['bin_edges'],
            range=hist_foci_info_yz['range'],
        )
        hist_foci_info_yz['collector'] += frame_hist
        # dna
        # # xy
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=hist_dna_info_xy['bin_edges'],
            range=hist_dna_info_xy['range'],
        )
        hist_dna_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=hist_dna_info_xz['bin_edges'],
            range=hist_dna_info_xz['range'],
        )
        hist_dna_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=hist_dna_info_yz['bin_edges'],
            range=hist_dna_info_yz['range'],
        )
        hist_dna_info_yz['collector'] += frame_hist
    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info,
            'Foci': r_hist_foci_info,
            'Dna': r_hist_dna_info
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': hist_crd_info_xy,
            'Mon': hist_mon_info_xy,
            'Foci': hist_foci_info_xy,
            'Dna': hist_dna_info_xy
        },
        'xzHist': {
            'Crd': hist_crd_info_xz,
            'Mon': hist_mon_info_xz,
            'Foci': hist_foci_info_xz,
            'Dna': hist_dna_info_xz
        },
        'yzHist': {
            'Crd': hist_crd_info_yz,
            'Mon': hist_mon_info_yz,
            'Foci': hist_foci_info_yz,
            'Dna': hist_dna_info_yz
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def hns_nucleoid_cub(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'nucleoid' atom
    group in the `geometry` of interest.

    Note
    ----
    In the HNS-DNA-Crowder project, we have a single semi-flexible ring
    polymer (called "bug") and several H-NS proteins (each formed of a core
    and two poles or patches) crowded by soft LJ spheres (crowders) in
    free ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to the
    center of mass of the single polymer (bug).

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = HnsCub(
        trajectory,
        lineage,
        'cubic',
        'nucleoid',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(
        sim_info.mmon * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'ndump', so trajectory dt is:
    sim_real_dt = sim_info.ndump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        unwrap_images=True,  # in_memory=True,
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups:
    bug = cell.select_atoms('resid 1')  # the bug
    hns_patch = cell.select_atoms('type 2')  # the hns patches
    # transformations:
    workflow = [
        mda_trans.unwrap(cell.atoms),
        mda_trans.center_in_box(bug),
        mda_trans.wrap(cell.atoms)
    ]
    cell.trajectory.add_transformations(*workflow)
    # defining collectors
    # bug:
    gyr_t = []
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = []
    shape_parameter_t = []
    # - bond info
    n_bonds = len(bug.bonds.indices)
    bond_lengths = np.zeros((n_bonds, 1), dtype=np.float64)
    cosine_corrs = np.zeros(n_bonds, dtype=np.float64)
    # mon and hns-patch are attracted if their distance <= the below distance:
    m_hpatch_attr_cutoff = 0.5 * (sim_info.dmon + sim_info.dhns_patch)
    m_hpatch_shape = (sim_info.nmon, 2 * sim_info.nhns)
    m_m_shape = (sim_info.nmon, sim_info.nmon)
    # distance matrices
    dist_m_hpatch = np.zeros(m_hpatch_shape, dtype=np.float64)
    dist_m_m = np.zeros(m_m_shape, dtype=np.float64)
    # contact matrices
    dir_contacts_m_hpatch = np.zeros(m_hpatch_shape, dtype=np.int64)
    dir_contacts_m_m = np.zeros(m_m_shape, dtype=np.int64)
    for _ in sliced_trj:
        # bug:
        # various measures of chain size
        gyr_t.append(bug.radius_of_gyration())
        # shape parameters:
        asphericity_t.append(bug.asphericity())
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes()]),
            axis=0
        )
        shape_parameter_t.append(bug.shape_parameter())
        # bond info
        bond_dummy, cosine_dummy = correlations.bond_info(
            bug,
            sim_info.topology
            )
        bond_lengths += bond_dummy
        cosine_corrs += cosine_dummy
        # bug - hns patch:
        # distance matrices
        dummy = mda_dist.distance_array(bug, hns_patch, box=cell.dimensions)
        dist_m_hpatch += dummy
        dummy_m_m = np.matmul(dummy, dummy.T)
        dist_m_m += dummy_m_m
        # contact matrices
        dummy = np.asarray(dummy <= m_hpatch_attr_cutoff, dtype=int)
        dir_contacts_m_hpatch += dummy
        dummy_m_m = np.matmul(dummy, dummy.T)
        dir_contacts_m_m += dummy_m_m
    # Saving collectors to memory
    # bug
    np.save(save_to + sim_name + '-gyrTMon.npy', np.array(gyr_t))
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    np.save(save_to + sim_name + '-asphericityTMon.npy',
            np.array(asphericity_t)
            )
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    bond_lengths = bond_lengths / n_frames
    bonds_per_lag = np.arange(n_bonds, 0, -1)
    cosine_corrs = cosine_corrs / (n_frames * bonds_per_lag)
    bond_lengths = bond_lengths.reshape(n_bonds,)
    np.save(save_to + sim_name + '-bondLengthVecMon.npy', bond_lengths)
    np.save(save_to + sim_name + '-bondCosineCorrVecMon.npy', cosine_corrs)
    # bug hns-patch:
    dist_m_hpatch = dist_m_hpatch / n_frames
    dist_m_m = dist_m_m / n_frames
    dir_contacts_m_hpatch = dir_contacts_m_hpatch / n_frames
    dir_contacts_m_m = dir_contacts_m_m / n_frames
    np.save(
        save_to + sim_name + "-distMatMonPatch.npy", dist_m_hpatch
    )
    np.save(
        save_to + sim_name + "-distMatMonMon.npy", dist_m_m
    )
    np.save(
        save_to + sim_name + "-directContactsMatMonPatch.npy",
        dir_contacts_m_hpatch
    )
    np.save(
        save_to + sim_name + "-directContactsMatTMonMon.npy",
        dir_contacts_m_m
    )
    # Simulation stamps:
    print('done.')


def hns_all_cub(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest,and saves a variety of
    outputs (mostly in the csv format) to the `save_to` directory.

    Note
    ----
    In the HNS-DNA-Crowder project, we have a single semi-flexible ring
    polymer (called "bug") and several H-NS proteins (each formed of a core
    and two poles or patches) crowded by soft LJ spheres (crowders) in
    free ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to the
    center of mass of the single polymer (bug).

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = HnsCub(
        trajectory,
        lineage,
        'cubic',
        'all',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {  # edges for distance r from the box center
            'bin_size':  0.1 * min(
                sim_info.dmon, sim_info.dhns, sim_info.dcrowd),
            'lmin': 0,
            'lmax': 0.5 * sim_info.lcube
            },
        'lEdge': {  # edges in cartesian coordinates
            'bin_size':  0.5 * min(
                sim_info.dmon, sim_info.dhns, sim_info.dcrowd),
            'lmin': -0.5 * sim_info.lcube,
            'lmax': 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(
        sim_info.mmon * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'ndump', so trajectory dt is:
    sim_real_dt = sim_info.ndump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        unwrap_images=True,  # in_memory=True,
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    bug = cell.select_atoms('resid 1')  # the bug
    hns = cell.select_atoms('type 3')  # the hns cores
    crds = cell.select_atoms('resid 0')  # crowders
    # transformations:
    workflow = [
        mda_trans.unwrap(cell.atoms),
        mda_trans.center_in_box(bug),
        mda_trans.wrap(cell.atoms)
    ]
    cell.trajectory.add_transformations(*workflow)
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    r_hist_hns_info = fixedsize_bins(
        sim_name,
        'rEdgeHns',
        bin_edges['rEdge']['bin_size'],
        bin_edges['rEdge']['lmin'],
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd_info['collector'].any() != 0,
            r_hist_hns_info['collector'].any() != 0,
            r_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            r_hist_hns_info['collector_std'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # bin edges and histograms in different directions:
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        'xEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        'yEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        'zEdge',
        bin_edges['lEdge']['bin_size'],
        bin_edges['lEdge']['lmin'],
        bin_edges['lEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            x_hist_info['collector'].any() != 0,
            x_hist_info['collector_std'].any() != 0,
            y_hist_info['collector'].any() != 0,
            y_hist_info['collector_std'].any() != 0,
            z_hist_info['collector'].any() != 0,
            z_hist_info['collector_std'].any() != 0,
            ]):
        raise ValueError(
            "One of the histogram collectors is not empty!")
    # 2D hists
    # crd
    # # xy
    hist_crd_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_crd_info_xy['collector'] = np.zeros(hist_crd_info_xy['n_bins'])
    hist_crd_info_xy['collector'] *= 0
    # # xz
    hist_crd_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_xz['collector'] = np.zeros(hist_crd_info_xz['n_bins'])
    hist_crd_info_xz['collector'] *= 0
    # # yz
    hist_crd_info_yz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_crd_info_yz['collector'] = np.zeros(hist_crd_info_yz['n_bins'])
    hist_crd_info_yz['collector'] *= 0
    # mon
    # # xy
    hist_mon_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_mon_info_xy['collector'] = np.zeros(hist_mon_info_xy['n_bins'])
    hist_mon_info_xy['collector'] *= 0
    # # xz
    hist_mon_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_xz['collector'] = np.zeros(hist_mon_info_xz['n_bins'])
    hist_mon_info_xz['collector'] *= 0
    # # yz
    hist_mon_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_mon_info_yz['collector'] = np.zeros(hist_mon_info_yz['n_bins'])
    hist_mon_info_yz['collector'] *= 0
    # hns
    # # xy
    hist_hns_info_xy = {
        'n_bins': (
            x_hist_info['n_bins'],
            y_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            y_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            y_hist_info['range']
        ]
    }
    hist_hns_info_xy['collector'] = np.zeros(hist_hns_info_xy['n_bins'])
    hist_hns_info_xy['collector'] *= 0
    # # xz
    hist_hns_info_xz = {
        'n_bins': (
            x_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            x_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            x_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_hns_info_xz['collector'] = np.zeros(hist_hns_info_xz['n_bins'])
    hist_hns_info_xz['collector'] *= 0
    # # yz
    hist_hns_info_yz = {
        'n_bins': (
            y_hist_info['n_bins'],
            z_hist_info['n_bins']
        ),
        'bin_edges': [
            y_hist_info['bin_edges'],
            z_hist_info['bin_edges']
        ],
        'range': [
            y_hist_info['range'],
            z_hist_info['range']
        ]
    }
    hist_hns_info_yz['collector'] = np.zeros(hist_hns_info_yz['n_bins'])
    hist_hns_info_yz['collector'] *= 0
    for _ in sliced_trj:
        # com of bug as the origin for histograms
        bug_com = bug.center_of_geometry()
        # crds
        pos = crds.positions - bug_com
        # # r
        dist = np.linalg.norm(pos, axis=1)
        dist_hist = frame_dist_hist(dist, r_hist_crd_info)
        r_hist_crd_info['collector'] += dist_hist
        r_hist_crd_info['collector_std'] += np.square(dist_hist)
        # # xy
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 1],
            bins=hist_crd_info_xy['bin_edges'],
            range=hist_crd_info_xy['range'],
        )
        hist_crd_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 2],
            bins=hist_crd_info_xz['bin_edges'],
            range=hist_crd_info_xz['range'],
        )
        hist_crd_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 1],
            pos[:, 2],
            bins=hist_crd_info_yz['bin_edges'],
            range=hist_crd_info_yz['range'],
        )
        hist_crd_info_yz['collector'] += frame_hist
        # bug
        # # r
        pos = bug.positions - bug_com
        dist = np.linalg.norm(pos, axis=1)
        dist_hist = frame_dist_hist(dist, r_hist_mon_info)
        r_hist_mon_info['collector'] += dist_hist
        r_hist_mon_info['collector_std'] += np.square(dist_hist)
        # # xy
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 1],
            bins=hist_mon_info_xy['bin_edges'],
            range=hist_mon_info_xy['range'],
        )
        hist_mon_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 2],
            bins=hist_mon_info_xz['bin_edges'],
            range=hist_mon_info_xz['range'],
        )
        hist_mon_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 1],
            pos[:, 2],
            bins=hist_mon_info_yz['bin_edges'],
            range=hist_mon_info_yz['range'],
        )
        hist_mon_info_yz['collector'] += frame_hist
        # hns
        # # r
        pos = hns.positions - bug_com
        dist = np.linalg.norm(pos, axis=1)
        dist_hist = frame_dist_hist(dist, r_hist_hns_info)
        r_hist_hns_info['collector'] += dist_hist
        r_hist_hns_info['collector_std'] += np.square(dist_hist)
        # # xy
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 1],
            bins=hist_hns_info_xy['bin_edges'],
            range=hist_hns_info_xy['range'],
        )
        hist_hns_info_xy['collector'] += frame_hist
        # # xz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 0],
            pos[:, 2],
            bins=hist_hns_info_xz['bin_edges'],
            range=hist_hns_info_xz['range'],
        )
        hist_hns_info_xz['collector'] += frame_hist
        # # yz
        frame_hist, _, _ = np.histogram2d(
            pos[:, 1],
            pos[:, 2],
            bins=hist_hns_info_yz['bin_edges'],
            range=hist_hns_info_yz['range'],
        )
        hist_hns_info_yz['collector'] += frame_hist
    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info,
            'Hns': r_hist_hns_info,
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': hist_crd_info_xy,
            'Mon': hist_mon_info_xy,
            'Hns': hist_hns_info_xy,
        },
        'xzHist': {
            'Crd': hist_crd_info_xz,
            'Mon': hist_mon_info_xz,
            'Hns': hist_hns_info_xz,
        },
        'yzHist': {
            'Crd': hist_crd_info_yz,
            'Mon': hist_mon_info_yz,
            'Hns': hist_hns_info_yz,
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def hns_nucleoid_cyl(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'nucleoid' atom
    group in the `geometry` of interest.

    Note
    ----
    In the HNS-DNA-Crowder project, we have a single semi-flexible ring
    polymer (called "bug") and several H-NS proteins (each formed of a core
    and two poles or patches) crowded by soft LJ spheres (crowders) in
    free ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to the
    center of mass of the single polymer (bug).

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        warnings.warn(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trajectories.",
            UserWarning
        )
    print("Setting the name of analyze file...")
    sim_info = HnsCyl(
        trajectory,
        lineage,
        'cylindrical',
        'nucleoid',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(
        sim_info.mmon * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'ndump', so trajectory dt is:
    sim_real_dt = sim_info.ndump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        unwrap_images=True,  # in_memory=True,
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups:
    bug = cell.select_atoms('resid 1')  # the bug
    hns_patch = cell.select_atoms('type 2')  # the hns patches
    hns_core = cell.select_atoms('type 3')  # the hns patches
    # defining collectors
    # bug:
    gyr_t = []
    fsd_t = []
    trans_size_t = []
    principal_axes_t = []
    asphericity_t = []
    shape_parameter_t = []
    # - bond info
    n_bonds = len(bug.bonds.indices)
    bond_lengths = np.zeros((n_bonds, 1), dtype=np.float64)
    cosine_corrs = np.zeros(n_bonds, dtype=np.float64)
    # distance matrices
    dist_m_hpatch_t = []
    dist_m_hcore_t = []
    dist_m_m_t = []
    for _ in sliced_trj:
        # bug:
        # various measures of chain size
        gyr_t.append(bug.radius_of_gyration())
        trans_size_t.append(transverse_size(bug))
        fsd_t.append(fsd(bug.positions))
        # shape parameters:
        asphericity_t.append(bug.asphericity(wrap=False, unwrap=False))
        shape_parameter_t.append(bug.shape_parameter(wrap=False))
        principal_axes_t.append(bug.principal_axes(wrap=False))
        # bond info
        bond_dummy, cosine_dummy = correlations.bond_info(
            bug,
            sim_info.topology
            )
        bond_lengths += bond_dummy
        cosine_corrs += cosine_dummy
        # bug - hns patch:
        # distance matrices
        dummy = mda_dist.distance_array(bug, hns_patch, box=cell.dimensions)
        dist_m_hpatch_t.append(dummy)
        dummy = mda_dist.distance_array(bug, hns_core, box=cell.dimensions)
        dist_m_hcore_t.append(dummy)
        dummy = mda_dist.distance_array(bug, bug, box=cell.dimensions)
        dist_m_m_t.append(dummy)
    # Saving collectors to memory
    # bug
    np.save(save_to + sim_name + '-transSizeTMon.npy', np.array(trans_size_t))
    np.save(save_to + sim_name + '-fsdTMon.npy', np.array(fsd_t))
    np.save(save_to + sim_name + '-gyrTMon.npy', np.array(gyr_t))
    np.save(
        save_to + sim_name + '-asphericityTMon.npy', np.array(asphericity_t)
        )
    np.save(save_to + sim_name + '-shapeTMon.npy', np.array(shape_parameter_t))
    np.save(
        save_to + sim_name + '-principalTMon.npy', np.array(principal_axes_t)
        )
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    # bond info
    bond_lengths = bond_lengths / n_frames
    bonds_per_lag = np.arange(n_bonds, 0, -1)
    cosine_corrs = cosine_corrs / (n_frames * bonds_per_lag)
    bond_lengths = bond_lengths.reshape(n_bonds,)
    np.save(save_to + sim_name + '-bondLengthVecMon.npy', bond_lengths)
    np.save(save_to + sim_name + '-bondCosineCorrVecMon.npy', cosine_corrs)
    # bug hns-patch:
    np.save(
        save_to + sim_name + "-distMatMonPatchT.npy", np.array(dist_m_hpatch_t)
        )
    np.save(save_to + sim_name + "-distMatMonCoreT.npy", np.array(dist_m_m_t))
    np.save(save_to + sim_name + "-distMatMonMonT.npy", np.array(dist_m_m_t))
    # Simulation stamps:
    print('done.')
