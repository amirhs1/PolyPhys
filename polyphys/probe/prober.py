from typing import Optional
import MDAnalysis as mda
import numpy as np
from polyphys.manage.parser import SumRule
from polyphys.manage.organizer import invalid_keyword


def log_outputs(log_file, geometry):
    cell_attrs = SumRule(log_file, geometry)
    output = f"N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}"
    details_out = output + "-log_details.csv"
    with open(details_out, 'w') as detailsfile:
        # neig_modify delay NUM every NUM check YES/NO:
        detailsfile.write('groupname,filename,ens,run_seg,rskin,delay,every,\
            check,epsilon,dcrowd,')
        detailsfile.write('ncrowd,lcyl,dcyl,nmon,total_time_s,cores,timestep,\
            atoms,ts_per_sec,')
        # Section columns: min time, avg time, max time, %varavg, %total"
        # Section rows: Pair, Bond, Neigh, Comm, Output, Modify, Other
        detailsfile.write('pair_avg_s,pair_pct,bond_avg_s,bond_pct,neigh_avg_s,\
            neigh_pct,comm_avg_s,')
        detailsfile.write('comm_pct,output_avg_s,output_pct,modify_avg_s,\
            modify_pct,other_avg_s,other_pct,dangerous\n')
    runtime_out = output + "-log_runtime.csv"
    with open(runtime_out, 'w') as runfile:
        runfile.write('groupname,filename,ncores,natoms,wall_time\n')
    return details_out, runtime_out


def lammps_log_details(log_files, details_out, runtime_out):
    for file in log_files:
        filename = file[0].split('.log')
        filename = filename[0]
        filename = filename.split('/')[-1]
        ens = filename.split('ens')[-1]
        groupname = filename.split('ens')[0]
        with open(file[0], 'r') as log,\
            open(details_out, 'a') as detailsfile,\
                open(runtime_out, 'a') as runfile:
            line = log.readline()
            # The other of while loop are important
            # neigh_modify delay every check page one
            j = 1
            while line:
                while line.startswith('variable'):
                    words = line.split()
                    line = log.readline()
                    if words[1] == 'epsilon1':
                        epsilon = words[3]
                    if words[1] == 'sig2':
                        dcrowd = words[3]
                    if words[1] == 'n_crowd':
                        ncrowd = words[3]
                    if words[1] == 'lz':
                        lcyl = str(2*float(words[3]))
                    if words[1] == 'r':
                        dcyl = str(2*float(words[3]))
                    if words[1] == 'n_bug':
                        nmon = words[3]
                if line.startswith('neighbor'):
                    words = line.split()
                    rskin = words[1].strip()  # rskin
                # neigh_modify delay every check page one
                if line.startswith('neigh_modify'):
                    words = line.split()
                    # picking the NUMs and Yes/No from neigh_modify
                    delay = words[2].strip()
                    every = words[4].strip()
                    check = words[6].strip()
                if line.startswith('Loop time'):
                    detailsfile.write(groupname)
                    detailsfile.write(",")
                    detailsfile.write(filename)
                    detailsfile.write(",")
                    detailsfile.write(ens)
                    detailsfile.write(",")
                    detailsfile.write(str(j))  # total time
                    detailsfile.write(",")
                    j += 1
                    # neighbor and neigh_modify ocurres\
                    # one time but other occures 15 times.
                    detailsfile.write(rskin)  # rskin
                    detailsfile.write(",")
                    detailsfile.write(delay)  # delay
                    detailsfile.write(",")
                    detailsfile.write(every)  # every
                    detailsfile.write(",")
                    detailsfile.write(check)  # check
                    detailsfile.write(",")
                    detailsfile.write(epsilon)  # epsilon
                    detailsfile.write(",")
                    detailsfile.write(dcrowd)  # dcrowd
                    detailsfile.write(",")
                    detailsfile.write(ncrowd)  # ncrowd
                    detailsfile.write(",")
                    detailsfile.write(lcyl)  # lcyl
                    detailsfile.write(",")
                    detailsfile.write(dcyl)  # dcyl
                    detailsfile.write(",")
                    detailsfile.write(nmon)  # nmon
                    detailsfile.write(",")
                    words = line.split()
                    detailsfile.write(words[3].strip())  # total time
                    detailsfile.write(",")
                    ncores = words[5].strip()
                    detailsfile.write(ncores)  # # of cores
                    detailsfile.write(",")
                    detailsfile.write(words[8].strip())  # total timesteps
                    detailsfile.write(",")
                    natoms = words[11].strip()
                    detailsfile.write(natoms)  # total atoms
                    detailsfile.write(",")
                if line.startswith('Performance:'):
                    words = line.split()
                    detailsfile.write(words[3].strip())  # timesteps per second
                    detailsfile.write(",")
                if line.startswith('Section'):
                    _ = log.readline()
                    for i in range(6):  \
                            # Section rows: Pair, Bond, Neigh, Comm, Output,\
                        # Modify, Other
                        # Section columns: min time, avg time, max time,\
                        # %varavg, %total"
                        line = log.readline()
                        sect_min = line.split('|')[2].strip()
                        detailsfile.write(sect_min)
                        detailsfile.write(",")

                        sect_pct = line.split()[-1]  # Pair pct of total time
                        detailsfile.write(sect_pct)
                        detailsfile.write(",")
                    line = log.readline()
                    sect_min = line.split('|')[2].strip()
                    detailsfile.write(sect_min)
                    detailsfile.write(",")
                    sect_pct = line.split()[-1]  # Pair pct of total time
                    detailsfile.write(sect_pct)
                    detailsfile.write(",")
                if line.startswith('Dangerous'):
                    words = line.split()
                    detailsfile.write(str(int(words[-1])))  \
                        # # number of dangerous builds
                    detailsfile.write("\n")
                # runtime files
                if line.startswith('Total wall time'):
                    runfile.write(groupname)
                    runfile.write(",")
                    runfile.write(filename)
                    runfile.write(",")
                    runfile.write(ncores)
                    runfile.write(",")
                    runfile.write(natoms)
                    runfile.write(",")
                    words = line.split()
                    runfile.write(words[-1])  # total wall time
                    runfile.write("\n")
                line = log.readline()


def stamps_report_with_measures(report_name, sim_info, n_frames, **measures):
    """
    writes a summary of stamps (properties and attributes) of a simulation \
    to file.

    `stamps_report` generates a dataset called `report_name` of the \
    values of some attributes of `sim_info`, the number of frames \
    `n_frames`, all the key and value pairs in all the given dictionaries \
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: SumRule object
        A SumRule object that contains information about the name, parents, \
        and physical attributes of a simulation.
    n_frames: int
        Number of frames/snapshots/configurations in a simulation.
    *measures: list of one or more dict
        A list of one or more measures. Each measure is a dictionary in which \
        each key and values pair are the name and value of a physical property.
    """
    with open(report_name, mode='w') as report:
        # write header
        for lineage_name in sim_info.genealogy:
            report.write(f"{lineage_name},")
        for attr_name in sim_info.attributes:
            report.write(f"{attr_name},")
        for measure in measures:  # each measure is a dict
            for stat_name in measure.keys():
                report.write(f"{stat_name},")
        report.write("n_frames\n")
        # write values
        for lineage_name in sim_info.genealogy:
            attr_value = getattr(sim_info, lineage_name)
            report.write(f"{attr_value},")
        for attr_name in sim_info.attributes:
            attr_value = getattr(sim_info, attr_name)
            report.write(f"{attr_value},")
        for measure in measures:
            for stat_value in measure.values():
                report.write(f"{stat_value},")
        report.write(f"{n_frames}")


def stamps_report(report_name, sim_info, n_frames):
    """
    writes a summary of stamps (properties and attributes) of a simulation \
    to file.

    `stamps_report` generates a dataset called `report_name` of the \
    values of some attributes of `sim_info`, the number of frames \
    `n_frames`, all the key and value pairs in all the given dictionaries \
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: SumRule object
        A SumRule object that contains information about the name, parents, \
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


def simple_stats(property_name, array):
    """
    measures the mean, standard deviation, variance, and standard error \
    of the mean (sem) for an array.

    Parameters
    ----------
    property_name : str
        Name of physical property.
    array: numpy array of float
        Values of `property_name`.

    Returns
    -------
    stats: dict
        an dict of mean, std, var, and sem of `property_name`
    """
    # Unbiased std, var, and sem.
    stat = {
        property_name + '_mean': np.mean(array),
        property_name + '_var': np.var(array, ddof=1),
        property_name + '_sem': np.std(array, ddof=1) / np.sqrt(len(array))
    }
    return stat


def end_to_end(positions):
    """
    measures the end-to-end distance of a linear polymer, in the frame of\
    reference located at the polymer's center of geometry.

    `positions` is sorted by atom id, so the end-to-end distance is simply \
    the difference between last and first items in `positions`/

    Parameters
    ----------
    positions : numpy ndarray of dtype float
        Positions (a natoms*ndim  array) of the N atoms within an atom group \
        in a frame or snapshot or time step. `positions` is sorted by atom \
        number form 1 to N.

    Returns
    -------
    edn_to_end: numppy array of dtype float
    """
    # calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    positions = positions[-1] - positions[0]
    return np.linalg.norm(positions)


def max_distance(positions):
    """
    measures the maximum ditance in each of the three Cartesian direction, \
    in the frame of reference located at the polymer's center of geometry.

    The maximum distance is computed by subtracting the max and min of all \
    the particle in an atom group in a given frame or snapshot or time step.

    Parameters
    ----------
    positions: numpy ndarray of dtype float
        Positions (a natoms*ndim  array) of the N atoms within an atom group \
        in a frame or snapshot or time step. `positions` is sorted by atom \
        number form 1 to N.

    Returns
    -------
    [xmax, ymax, zmax]: numpy array of dtype float
    """
    # calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    xmax = np.abs(np.amax(positions[:, 0]) - np.amin(positions[:, 0]))
    ymax = np.abs(np.amax(positions[:, 1]) - np.amin(positions[:, 1]))
    zmax = np.abs(np.amax(positions[:, 2]) - np.amin(positions[:, 2]))
    return np.array([xmax, ymax, zmax])


def fsd(
    positions: np.ndarray,
    axis=2
) -> np.ndarray:
    """
    calculates the average size/diameter of a polymer confined in a \
    cylindrical geometry based on the farthermost distance concept.

    fsd stands for Feret's statistical diameter: other names are the mean \
    span dimension, the farthermost distance, or the mean caliper diameter.

    Parameters
    ----------
    positions: numpy ndarray of dtype float
        Positions (a natoms*ndim  array) of the N atoms within an atom group \
        in a frame or snapshot or time step. `positions` is sorted by atom \
        number form 1 to N.
    axis: int or tuple of int, default 2
        The index of the axis of the cylinder; by default it is the z axis.\
        It can be any integer in the range of the spatial dimension of the \
        system.

    Returns
    -------
    fsd: numpy array of dtype float

    References
    ----------
    "A Theoretical Study of the Separation Principle in Size Exclusion \
    Chromatography", Wang Y Teraoka
    I Hansen FY Peters GH Ole H. Macromolecules 2010, 43, 3, 1651-1659
    https://doi.org/10.1021/ma902377g
    """
    # calculation in the center of geometry of the atom group:
    positions = positions - np.mean(positions, axis=0)
    positions = np.ptp(positions[:, axis])
    return positions


def bin_create(sim_name, edge_name, bin_size, lmin, lmax, save_to):
    """
    generates arrays of bins and histograms

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

    Returns
    -------
    bin_edges : numpy array of dtype float
        The edges to pass into a histogram. Save `bin_edges` to file if \
        `save_to` is not None.
    hist: array of dtype int
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
    bin_type: str = 'ordinaray',
    save_to: Optional[str] = None,
    ):
    """
    generates arrays of bins and histograms, ensuring that the `bin_size` \
    guaranteed. To achieve this, it extend the `lmin` and `lmax` limits.

    To-do List
    ----------
    1. Following the idea used: 
    https://docs.mdanalysis.org/1.1.1/_modules/MDAnalysis/lib/util.html#fixedwidth_bins 
    Makes input array-like so bins can be calculated for 1D data (then all parameters 
    are simple floats) or nD data (then parameters are supplied as arrays, with each entry
    correpsonding to one dimension).
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
    bin_type: {'ordinary', 'nonnegative', 'periodic'}, defualt 'ordinary'
        The type of bin in a given direction in a given coordinate system:

        'ordinary'
            A bounded or unbounded coordinate such as any of the cartesian coordinates or \
            or the polar coordinate in the spherical coordinate system. For such coordinates, \
            the `lmin` and `lmax` limits are equaly extended to ensure the `bin_size`.

        'nonnegative'
            A nonnegative coordinate such as the r direction in the polar or spherical \
            coordinate system. For such coordinates, ONLY `lmax` limit is extended to \
            ensure the `bin_size`.

        'periodic'
            A periodic coordinate such as the azimuthal direction in the spherical coordinate. \
            For such coordinates, the initial `lmin` and `lmax` limits are used. 

    save_to : str, default None
        Whether save outputs to memory as csv files or not.

    Returns
    -------
    bin_edges : numpy array of dtype float
        The edges to pass into a histogram. Save `bin_edges` to file if \
        `save_to` is not None.
    hist: array of dtype int
        An empty histogram

    Reference:
    https://docs.mdanalysis.org/1.1.1/documentation_pages/lib/util.html#MDAnalysis.analysis.density.fixedwidth_bins
    """
    _bin_types = ['ordinary', 'nonnagative', 'periodic']
    if not np.all(lmin < lmax):
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = np.asarray(bin_size, dtype=np.float_)
    _lmin = np.asarray(lmin, dtype=np.float_)
    _lmax = np.asarray(lmax, dtype=np.float_)
    _length = _lmax - _lmin

    if bin_type == 'ordinary':
        nbins = np.ceil(_length / _delta).astype(np.int_)  # number of bins
        dx = 0.5 * (nbins * _delta - _length)  # add half of the excess to each end
        _lmin = _lmin - dx 
        _lmax =  _lmax + dx
    elif bin_type == 'nonnegative':
        nbins = np.ceil(_length / _delta).astype(np.int_)
        dx = 0.5 * (nbins * _delta - _length)
        _lmax =  _lmax + 2 * dx
    elif bin_type == 'periodic':
        periodic_bin_edges = np.arange(_lmin, _lmax + _delta, _delta)
    else:
        invalid_keyword(bin_type, _bin_types)
    # create empty grid with the right dimensions (and get the edges)
    hist_collectors, bin_edges = np.histogram(
        np.zeros(1),
        bins=nbins,            
        range=(_lmin, _lmax)
    )
    if bin_type == 'periodic': # improved this.
        bin_edges = periodic_bin_edges
    hist_collectors *= 0
    np.save(save_to + sim_name + '-' + edge_name + '.npy', bin_edges)
    return bin_edges, hist_collectors


def probe_bug_with_histogram(
    topology,
    trajectory,
    geometry,
    lineage,
    save_to="./",
    continuous=False
) -> None:
    """
    runs various analyses on a `lineage` simulation of a 'bug' atom group in \
    the `geometry` of interest.

    Parameters
    ----------
    toplogy: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}
        Shape of the simulation box.
    lineage: {'segment', 'whole'}
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory \
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        print(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trjectories.")
    print("Setting the name of analyze file...")
    sim_info = SumRule(
        trajectory,
        geometry=geometry,
        group='bug',
        lineage=lineage
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {
            'bin_size':  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            'lmax': 0.5 * sim_info.dcyl
            },
        'zEdge': {
            'bin_size':  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            'lmax': 0.5 * sim_info.lcyl
            },
        'thetaEdge': {
            'bin_size':  np.pi / 36,
            'lmax': np.pi
            },
        'rfloryEdge': {
            'bin_size': 0.5 * sim_info.dmon,
            'lmax': sim_info.nmon * sim_info.dmon
            }
        }
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
    # bin edges and histograms in different directions:
    r_edges, r_hist = fixedsize_bins(
        sim_name, 'rEdge', bin_edges['rEdge']['bin_size'],
        0.0, bin_edges['rEdge']['lmax'], save_to)
    z_edges, z_hist = fixedsize_bins(
        sim_name, 'zEdge', bin_edges['zEdge']['bin_size'],
        -1 * bin_edges['zEdge']['lmax'], bin_edges['zEdge']['lmax'], save_to)
    theta_edges, theta_hist = fixedsize_bins(
        sim_name, 'thetaEdge', bin_edges['thetaEdge']['bin_size'],
        -1 * bin_edges['thetaEdge']['lmax'], bin_edges['thetaEdge']['lmax'],
        save_to)  # in radians
    # distribution of the size of the end-to-end
    rflory_edges, rflory_hist = fixedsize_bins(
        sim_name, 'rfloryEdge', bin_edges['rfloryEdge']['bin_size'],
        0, bin_edges['rfloryEdge']['lmax'], save_to)
    fsd_t = np.empty([0])
    rflory_t = np.empty([0])
    gyr_t = np.empty([0])
    if any([
            r_hist.any() != 0, z_hist.any() != 0,
            theta_hist.any() != 0, rflory_hist.any() != 0
            ]):
        raise ValueError("One of the histogram collectors are not empty!")
    for _ in sliced_trj:
        # various measures of chain size
        fsd_t = np.append(fsd_t, np.array([fsd(bug.positions)]), axis=0)
        gyr_t = np.append(gyr_t, np.array([bug.radius_of_gyration()]), axis=0)
        rms = end_to_end(bug.positions)
        rflory_t = np.append(rflory_t, np.array([rms]), axis=0)
        dummy_hist, _ = np.histogram(rms, rflory_edges)
        # RDF of the end-to-end distance
        rflory_hist = np.add(rflory_hist, dummy_hist)
        # number density in the cell's frame of reference
        # histogram in r direction
        rmon = np.linalg.norm(bug.positions[:, :2], axis=1)
        dummy_hist, _ = np.histogram(rmon, r_edges)
        r_hist = np.add(r_hist, dummy_hist)
        # histogram in z direction
        zmon = bug.positions[:, 2]
        dummy_hist, _ = np.histogram(zmon, z_edges)
        z_hist = np.add(z_hist, dummy_hist)
        # histogram in theta
        theta = np.arctan2(
                    bug.positions[:, 1], bug.positions[:, 0])  # in radians
        dummy_hist, _ = np.histogram(theta, theta_edges)
        theta_hist = np.add(theta_hist, dummy_hist)
    np.save(save_to + sim_name + '-rHistMon.npy', r_hist)
    np.save(save_to + sim_name + '-zHistMon.npy', z_hist)
    np.save(save_to + sim_name + '-thetaHistMon.npy', theta_hist)
    np.save(save_to + sim_name + '-rfloryHistMon.npy', rflory_hist)
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def probe_bug(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './',
    continuous: bool = False
) -> None:
    """
    runs various analyses on a `lineage` simulation of a 'bug' atom group in \
    the `geometry` of interest.

    Parameters
    ----------
    toplogy: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, defualt 'biaxial'
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory \
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        print(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trjectories.")
    print("Setting the name of analyze file...")
    sim_info = SumRule(
        trajectory,
        geometry=geometry,
        group='bug',
        lineage=lineage
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
    # dict of bin edges:
    bin_edges = {
        'rfloryEdge': {
            'bin_size': 0.5 * sim_info.dmon,
            'lmax': sim_info.nmon * sim_info.dmon
            }
        }
    # distribution of the size of the end-to-end
    rflory_edges, rflory_hist = fixedsize_bins(
        sim_name,
        'rfloryEdgeMon',
        bin_edges['rfloryEdge']['bin_size'],
        0,
        bin_edges['rfloryEdge']['lmax'],
        bin_type='nonnegative', 
        save_to=save_to
    )
    fsd_t = np.empty([0])
    rflory_t = np.empty([0])
    gyr_t = np.empty([0])
    if rflory_hist.any() != 0:
        raise ValueError("One of the histogram collectors are not empty!")
    for _ in sliced_trj:
        # various measures of chain size
        fsd_t = np.append(fsd_t, np.array([fsd(bug.positions)]), axis=0)
        gyr_t = np.append(gyr_t, np.array([bug.radius_of_gyration()]), axis=0)
        rms = end_to_end(bug.positions)
        rflory_t = np.append(rflory_t, np.array([rms]), axis=0)
        dummy_hist, _ = np.histogram(rms, rflory_edges)
        # RDF of the end-to-end distance
        rflory_hist = np.add(rflory_hist, dummy_hist)
    np.save(save_to + sim_name + '-rfloryHistMon.npy', rflory_hist)
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def rmsd_bug(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './'
) -> None:
    """
    comptues the rmsd of a 'segment simulation of a 'bug' atom group in the \
    `geometry` of interest, and then saves the output to the `save_to` \
    directory.

    `rmsd_bug` does not support the `continuous` option defined in \
    `probe_bug` and `probe_all`.

    Parameters
    ----------
    toplogy: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
        Type of the input file.
    save_to: str, default './'
        The absolute/relative path of a directory to which output is saved.
    """
    sim_info = SumRule(
        trajectory,
        geometry=geometry,
        group='bug',
        lineage=lineage
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print(f"Doing RMSD analysis on {sim_name} ...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(sim_info.mmon * sim_info.eps_others) \
        # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
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
    cell.transfer_to_memory(step=50, verbose=False)
    _ = mda.analysis.align.AlignTraj(
        cell,
        cell,
        select='resid 1',
        filename=sim_name + '.dcd'
    ).run()
    matrix = mda.analysis.diffusionmap.DistanceMatrix(
        cell,
        select='resid 1'
    ).run()
    np.save(save_to + sim_name + '-rmsdMatrixMon.npy', matrix.dist_matrix)


def probe_all(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """
    runs various analyses on a `lineage` simulation of an 'all' atom \
    group in the `geometry` of interest, and saves a variety of \
    outputs (mostly in the csv format) to the `save_to` directory.

    Parameters
    ----------
    toplogy: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
        Type of the input file.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    continuous: bool, default False
        Whether a `trajectory` file is a part of a sequence of trajectory \
        segments or not.
    """
    if (lineage == 'segment') & (continuous is False):
        print(
            "lineage is "
            f"'{lineage}' "
            "and 'continuous' is "
            f"'{continuous}. "
            "Please ensure the "
            f"'{trajectory}' is NOT part of a sequence of trjectories.")
    print("Setting the name of analyze file...\n")
    sim_info = SumRule(
        trajectory,
        geometry=geometry,
        group='all',
        lineage=lineage
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...")
    # dict of bin edges:
    bin_edges = {
        'rEdge': {
            'bin_size':  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            'lmax': 0.5 * sim_info.dcyl
            },
        'zEdge': {
            'bin_size':  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            'lmax': 0.5 * sim_info.lcyl
            },
        'thetaEdge': {
            'bin_size':  np.pi / 36,
            'lmax': np.pi
            }
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(sim_info.mmon * sim_info.eps_others) \
        # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
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
    # radial direction of the cylindrical coordinate system
    r_edges, r_hist_crd = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges['rEdge']['bin_size'],
        0.0,
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    _, r_hist_mon = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges['rEdge']['bin_size'],
        0.0,
        bin_edges['rEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    # z direction of the cylindrical coordinate system
    z_edges, z_hist_crd = fixedsize_bins(
        sim_name,
        'zEdgeCrd',
        bin_edges['zEdge']['bin_size'],
        -1 * bin_edges['zEdge']['lmax'],
        bin_edges['zEdge']['lmax'],
        bin_type='ordinary',
        save_to=save_to
    )
    _, z_hist_mon = fixedsize_bins(
        sim_name,
        'zEdgeMon',
        bin_edges['zEdge']['bin_size'],
        -1 * bin_edges['zEdge']['lmax'],
        bin_edges['zEdge']['lmax'],
         bin_type='ordinary',
        save_to=save_to
    )
    # theta of the cylindrical coordinate system
    theta_edges, theta_hist_crd = fixedsize_bins(
        sim_name, 'thetaEdgeCrd',
        bin_edges['thetaEdge']['bin_size'],
        -1 * bin_edges['thetaEdge']['lmax'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
    )  # in radians
    _, theta_hist_mon = fixedsize_bins(
        sim_name,
        'thetaEdgeMon',
        bin_edges['thetaEdge']['bin_size'],
        -1 * bin_edges['thetaEdge']['lmax'],
        bin_edges['thetaEdge']['lmax'],
        bin_type='periodic',
        save_to=save_to
    )  # in radians
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd.any() != 0, r_hist_mon.any() != 0,
            z_hist_crd.any() != 0, z_hist_mon.any() != 0,
            theta_hist_crd.any() != 0, theta_hist_mon.any() != 0
            ]):
        raise ValueError(
            "One of the histogram collectors are not empty!")
    for _ in sliced_trj:
        # histogram in r direction
        rpos = np.linalg.norm(crds.positions[:, :2], axis=1)
        dummy_hist, _ = np.histogram(rpos, r_edges)
        r_hist_crd = np.add(r_hist_crd, dummy_hist)
        rpos = np.linalg.norm(bug.positions[:, :2], axis=1)
        dummy_hist, _ = np.histogram(rpos, r_edges)
        r_hist_mon = np.add(r_hist_mon, dummy_hist)
        # histogram in z direction
        zpos = crds.positions[:, 2]
        dummy_hist, _ = np.histogram(zpos, z_edges)
        z_hist_crd = np.add(z_hist_crd, dummy_hist)
        zpos = bug.positions[:, 2]
        dummy_hist, _ = np.histogram(zpos, z_edges)
        z_hist_mon = np.add(z_hist_mon, dummy_hist)
        # histogram in theta
        theta = np.arctan2(
            crds.positions[:, 1],
            crds.positions[:, 0]
        )  # in degrees
        dummy_hist, _ = np.histogram(theta, theta_edges)
        theta_hist_crd = np.add(theta_hist_crd, dummy_hist)
        theta = np.arctan2(
            bug.positions[:, 1],
            bug.positions[:, 0]
        )  # in degrees
        dummy_hist, _ = np.histogram(theta, theta_edges)
        theta_hist_mon = np.add(theta_hist_mon, dummy_hist)
    lastname = 'Crd'
    np.save(save_to + sim_name + '-rHist' + lastname + '.npy', r_hist_crd)
    np.save(save_to + sim_name + '-zHist' + lastname + '.npy', z_hist_crd)
    np.save(
        save_to + sim_name + '-thetaHist' + lastname + '.npy', theta_hist_crd
    )
    lastname = 'Mon'
    np.save(save_to + sim_name + '-rHist' + lastname + '.npy', r_hist_mon)
    np.save(save_to + sim_name + '-zHist' + lastname + '.npy', z_hist_mon)
    np.save(
        save_to + sim_name + '-thetaHist' + lastname + '.npy', theta_hist_mon
    )
    print('done.')
