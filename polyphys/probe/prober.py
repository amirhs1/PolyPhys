from typing import Optional, Tuple, Dict, Type
import MDAnalysis as mda
import numpy as np
import pandas as pd
import os
from itertools import combinations
from polyphys.manage.parser import SumRule, TransFoci
from polyphys.manage.organizer import invalid_keyword


def log_stat(
    logpath: str,
    save_to: str = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Parses a LAMMPS log file loacted at the `logpath` to extract the
    performance information about a simulation in a given project.

    Parameters
    ----------
    logpath: str
        Path to a log file
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    run_stat: pd.DataFrame
        A dataframe in which each row is equivalent to one run loop in a
        LAMMPS input file (a LAMMPS input file can has several "run" commands
        and thus several performance breakdowns. The columns in the dataframe
        are about the parallel efficiency, speed, performance breakdown, and
        dagerousn buils of neighbor lists.
    wall_time_stat: pd.DataFrame
        A dataframe with one row in which the columns are the log name, number
        of atoms, number of cores (cpus), and the total wall time.
    """
    logname, _ = os.path.splitext(logpath)
    logname = logname.split("/")[-1]
    run_stat = {
        "log_name": [],
        "loop_id": [],
        "loop_time_second": [],
        "n_cores": [],
        "loop_steps": [],
        "n_atoms": [],
        "step_per_second": [],
        "mpi_task": [],
        "openmp_threads": [],
        "pair_avg_s": [],
        "pair_pct": [],
        "bond_avg_s": [],
        "bond_pct": [],
        "neigh_avg_s": [],
        "neigh_pct": [],
        "comm_avg_s": [],
        "comm_pct": [],
        "output_avg_s": [],
        "output_pct": [],
        "modify_avg_s": [],
        "modify_pct": [],
        "other_avg_s": [],
        "other_pct": [],
        "dangerous_builds": []
    }
    wall_time_stat = {
        "log_name": [],
        "n_cores": [],
        "n_atoms": [],
        "wall_time": []
    }
    with open(logpath, 'r') as logfile:
        lines = logfile.readlines()
        # neigh_modify delay every check page one
        j = 1  # if j > 1, then there is a loop in the lammmps input file.
        for idx, line in enumerate(lines):
            if line.startswith('Loop time'):
                run_stat["log_name"].append(logname)  # total time
                run_stat["loop_id"].append(j)  # total time
                j += 1
                words = line.split()
                run_stat["loop_time_second"].append(float(words[3].strip()))
                n_cores = int(words[5].strip())
                run_stat["n_cores"].append(n_cores)
                run_stat["loop_steps"].append(int(words[8].strip()))
                n_atoms = int(words[11].strip())
                run_stat["n_atoms"].append(n_atoms)
            if line.startswith('Performance:'):
                words = line.split()
                run_stat["step_per_second"].append(float(words[3].strip()))
                words = lines[idx+1].split()  # next line
                run_stat["mpi_task"].append(int(words[4].strip()))
                run_stat["openmp_threads"].append(int(words[8].strip()))
            if line.startswith("MPI task timing breakdown:"):
                words = lines[idx+3].split("|")
                run_stat["pair_avg_s"].append(float(words[2].strip()))
                run_stat["pair_pct"].append(float(words[5].strip()))
                words = lines[idx+4].split("|")
                run_stat["bond_avg_s"].append(float(words[2].strip()))
                run_stat["bond_pct"].append(float(words[5].strip()))
                words = lines[idx+5].split("|")
                run_stat["neigh_avg_s"].append(float(words[2].strip()))
                run_stat["neigh_pct"].append(float(words[5].strip()))
                words = lines[idx+6].split("|")
                run_stat["comm_avg_s"].append(float(words[2].strip()))
                run_stat["comm_pct"].append(float(words[5].strip()))
                words = lines[idx+7].split("|")
                run_stat["output_avg_s"].append(float(words[2].strip()))
                run_stat["output_pct"].append(float(words[5].strip()))
                words = lines[idx+8].split("|")
                run_stat["modify_avg_s"].append(float(words[2].strip()))
                run_stat["modify_pct"].append(float(words[5].strip()))
                words = lines[idx+9].split("|")
                run_stat["other_avg_s"].append(float(words[2].strip()))
                run_stat["other_pct"].append(float(words[5].strip()))
            if line.startswith('Dangerous'):
                words = line.split()
                run_stat["dangerous_builds"].append(int(words[-1]))
            if line.startswith('Total wall time'):
                wall_time_stat['log_name'].append(logname)
                wall_time_stat['n_cores'].append(n_cores)
                wall_time_stat['n_atoms'].append(n_atoms)
                words = line.split()
                # total wall time:
                wall_time_stat['wall_time'].append(words[-1])
        run_stat = pd.DataFrame.from_dict(run_stat)
        run_stat['filename'] = logname
        wall_time_stat = pd.DataFrame.from_dict(wall_time_stat)
        wall_time_stat['filename'] = logname
        if save_to is not None:
            run_stat.to_csv(logname + '-runStat.csv', index=False)
            wall_time_stat.to_csv(logname + "-runWallTime.csv", index=False)
        return run_stat, wall_time_stat


def thermo_multi(
    logpath: str,
    save_to: str = None
) -> pd.DataFrame:
    """Parses a LAMMPS `log` file written with "multi" thermo style to extract
    the thermodynamic data about a simulation.

    The thermodyanmic information can be related to multiple "run" commands in
    a LAMMPS input file.

    Parameters
    ----------
    logpath: str
        Path to a log file
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    thermo: pd.DataFrame
        A dataframe in which the columns are the various thermodynamic
        quantities of the "multi" thermo style.
    """
    logname, _ = os.path.splitext(logpath)
    logname = logname.split("/")[-1]
    thermo = {
        "Step": [],
        "TotEng": [],
        "KinEng": [],
        "Temp": [],
        "PotEng": [],
        "E_bond": [],
        "E_angle": [],
        "E_dihed": [],
        "E_impro": [],
        "E_vdwl": [],
        "E_coul": [],
        "E_long": [],
        "Press": []
    }
    with open(logpath, 'r') as logfile:
        lines = logfile.readlines()
        for idx, line in enumerate(lines):
            if line.startswith("---------------- Step"):
                words = line.split()
                thermo["Step"].append(int(words[2]))  # Step
                next_line = lines[idx+1]
                # the "els .. break" statement ensures all the 4 lines after
                # this line belong to the same time step.
                if next_line.startswith("TotEng"):
                    words = next_line.split()
                    thermo["TotEng"].append(float(words[2]))  # TotEng
                    thermo["KinEng"].append(float(words[5]))  # KinEng
                    thermo["Temp"].append(float(words[8]))  # Temp
                else:
                    break
                next_line = lines[idx+2]
                if next_line.startswith("PotEng"):
                    words = next_line.split()
                    thermo["PotEng"].append(float(words[2]))  # PotEng
                    thermo["E_bond"].append(float(words[5]))  # E_bond
                    thermo["E_angle"].append(float(words[8]))  # E_angle
                else:
                    break
                next_line = lines[idx+3]
                if next_line.startswith("E_dihed"):
                    words = next_line.split()
                    thermo["E_dihed"].append(float(words[2]))  # E_coul
                    thermo["E_impro"].append(float(words[5]))  # E_impro
                    thermo["E_vdwl"].append(float(words[8]))  # E_vdwl
                else:
                    break
                next_line = lines[idx+4]
                if next_line.startswith("E_coul"):
                    words = next_line.split()
                    thermo["E_coul"].append(float(words[2]))  # E_coul
                    thermo["E_long"].append(float(words[5]))  # E_long
                    thermo["Press"].append(float(words[8]))  # Press
                else:
                    break
    thermo = pd.DataFrame.from_dict(thermo)
    thermo.drop_duplicates(inplace=True)
    thermo.reset_index(inplace=True, drop=True)
    thermo['filename'] = logname
    if save_to is not None:
        thermo.to_csv(logname + '-thermo.csv', index=False)
    return thermo


def thermo_one_and_custom(
    logpath: str,
    save_to: str = None
) -> pd.DataFrame:
    """Parses a LAMMPS `log` file written with "one" or "custom" thermo styles
    to extractthe thermodynamic data about a simulation.

    The thermodyanmic information can be related to multiple "run" commands in
    a LAMMPS input file.

    Parameters
    ----------
    logpath: str
        Path to a log file
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    thermo: pd.DataFrame
        A dataframe in which the columns are the various thermodynamic
        quantities of the "one" thermo style or the custom thermodynamic
        quanitites of the "custom" thermo style.
    """
    logname, _ = os.path.splitext(logpath)
    logname = logname.split("/")[-1]
    thermo = {}
    with open(logpath, 'r') as logfile:
        line = logfile.readline()
        loop_counter = 1  # counts the # of runs/loops in a log
        while(line):
            if line.startswith("Step"):
                if loop_counter == 1:
                    headers = [word.strip() for word in line.split()]
                    thermo = {header: [] for header in headers}
                loop_counter += 1
                line = logfile.readline()
                while(line):
                    if line.startswith("Loop time of "):
                        break
                    try:
                        words = [float(i) for i in line.split()]
                        for key, value in zip(thermo.keys(), words):
                            if key == 'Step':
                                thermo[key].append(int(value))  # Step
                            else:
                                thermo[key].append(value)
                    except ValueError:
                        # Some text has snuck into the thermo output,
                        # this text is ignored.
                        pass
                    line = logfile.readline()
            else:
                line = logfile.readline()
    thermo = pd.DataFrame.from_dict(thermo)
    thermo.drop_duplicates(inplace=True)
    thermo['filename'] = logname
    if save_to is not None:
        thermo.to_csv(logname + '-thermo.csv', index=False)
    return thermo


def thermo_one(
    logpath: str,
    save_to: str = None
) -> pd.DataFrame:
    """Parses a LAMMPS `log` file written with "one" or "custom" thermo styles
    to extractthe thermodynamic data about a simulation.

    The thermodyanmic information can be related to multiple "run" commands in
    a LAMMPS input file.

    Parameters
    ----------
    logpath: str
        Path to a log file
    save_to : str, default None
        Absolute/relative path of a directory to which outputs are saved.

    Return
    ------
    thermo: pd.DataFrame
        A dataframe in which the columns are the various thermodynamic
        quantities of the "one" thermo style.
    """
    logname, _ = os.path.splitext(logpath)
    logname = logname.split("/")[-1]
    thermo = {
        "Step": [],
        "Temp": [],
        "E_pair": [],
        "E_mol": [],
        "TotEng": [],
        "Press": []
    }
    with open(logpath, 'r') as logfile:
        line = logfile.readline()
        while(line):
            if line.startswith("Step Temp E_pair E_mol TotEng Press"):
                line = logfile.readline()
                while(line):
                    if line.startswith("Loop time of "):
                        break
                    try:
                        words = [float(i) for i in line.strip().split()]
                        thermo["Step"].append(int(words[0]))  # Step
                        thermo["Temp"].append(float(words[1]))  # Temp
                        thermo["E_pair"].append(float(words[2]))  # E_pair
                        thermo["E_mol"].append(float(words[3]))  # E_mol
                        thermo["TotEng"].append(float(words[4]))  # TotEng
                        thermo["Press"].append(float(words[5]))  # Press
                    except ValueError:
                        # Some text has snuck into the thermo output,
                        # this text is ignored.
                        pass
                    line = logfile.readline()
            else:
                line = logfile.readline()
    thermo = pd.DataFrame.from_dict(thermo)
    thermo.drop_duplicates(inplace=True)
    thermo['filename'] = logname
    if save_to is not None:
        thermo.to_csv(logfile + '-thermo.csv', index=False)
    return thermo


def stamps_report_with_measures(
    report_name: str,
    sim_info: Type[SumRule],
    n_frames: int,
    measures: Dict[str, float],
) -> None:
    """Writes a summary of stamps (properties and attributes) of a simulation
    to file.

    `stamps_report` generates a dataset called `report_name` of the
    values of some attributes of `sim_info`, the number of frames
    `n_frames`, all the key and value pairs in all the given dictionaries
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: SumRule object
        A SumRule object that contains information about the name, parents,
        and physical attributes of a simulation.
    n_frames: int
        Number of frames/snapshots/configurations in a simulation.
    *measures: dict
        A dictionary of measures where a key and value pair is the name and
        value of a physical property.
    """
    with open(report_name, mode='w') as report:
        # write header
        for lineage_name in sim_info.genealogy:
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


def stamps_report(report_name, sim_info, n_frames):
    """Writes a summary of stamps (properties and attributes) of a simulation
    to file.

    `stamps_report` generates a dataset called `report_name` of the
    values of some attributes of `sim_info`, the number of frames
    `n_frames`, all the key and value pairs in all the given dictionaries
    `measures`.

    Parameters
    ----------
    report_name: str
        Name of the report.
    sim_info: SumRule object
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


def simple_stats(property_name, array):
    """Measures the mean, standard deviation, variance, and standard error
    of the mean (sem) for an array.

    Parameters
    ----------
    property_name : str
        Name of physical property.
    array: numpy array of float
        Values of `property_name`.

    Return
    ------
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
    """Measures the end-to-end distance of a linear polymer, in the frame of
    reference located at the polymer's center of geometry.

    `positions` is sorted by atom id, so the end-to-end distance is simply
    the difference between last and first items in `positions`

    Parameters
    ----------
    positions : numpy array of float
        Positions (a n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.

    Return
    ------
    end_to_end: numpy array of float
    """
    # calculation in the center of geometry of the atom group.
    positions = positions - np.mean(positions, axis=0)
    positions = positions[-1] - positions[0]
    return np.linalg.norm(positions)


def max_distance(positions):
    """Measures the maximum distance in each of the three Cartesian direction,
    in the frame of reference located at the polymer's center of geometry.

    The maximum distance is computed by subtracting the max and min of all
    the particle in an atom group in a given frame or snapshot or time step.

    Parameters
    ----------
    positions: numpy array of dtype float
        Positions (a n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.

    Return
    ------
    [xmax, ymax, zmax]: numpy array of  float
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
    """Calculates the average size/diameter of a polymer confined in a
    cylindrical geometry based on the farthermost distance concept.

    fsd stands for Feret's statistical diameter: other names are the mean
    span dimension, the farthermost distance, or the mean caliper diameter.

    Parameters
    ----------
    positions: numpy array of  float
        Positions (a n_atoms*n_dim  array) of the N atoms within an atom
        group in a frame or snapshot or time step. `positions` is sorted
        by atom number form 1 to N.
    axis: int or tuple of int, default 2
        The index of the axis of the cylinder; by default it is the z axis.
        It can be any integer in the range of the spatial dimension of the
        system.

    Return
    ------
    fsd: numpy array of  float

    References
    ----------
    "A Theoretical Study of the Separation Principle in Size Exclusion
    Chromatography", Wang Y Teraoka
    I Hansen FY Peters GH Ole H. Macromolecules 2010, 43, 3, 1651-1659
    https://doi.org/10.1021/ma902377g
    """
    # calculation in the center of geometry of the atom group:
    positions = positions - np.mean(positions, axis=0)
    positions = np.ptp(positions[:, axis])
    return positions


def bin_create(sim_name, edge_name, bin_size, lmin, lmax, save_to):
    """Generates arrays of bins and histograms

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
) -> dict:
    """Generates arrays of bins and histograms, ensuring that the `bin_size`
    guaranteed. To achieve this, it extend the `lmin` and `lmax` limits.

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
        Whether save outputs to memory as csv files or not.

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
    _delta = np.asarray(bin_size, dtype=np.float_)
    _lmin = np.asarray(lmin, dtype=np.float_)
    _lmax = np.asarray(lmax, dtype=np.float_)
    _length = _lmax - _lmin
    if bin_type == 'ordinary':
        n_bins = np.ceil(_length / _delta).astype(np.int_)  # number of bins
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
        n_bins = np.ceil(_length / _delta).astype(np.int_)
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
        n_bins = np.ceil(_length / _delta).astype(np.int_)  # number of bins
        print(
            f" number of bins '{n_bins}'"
            " is more than or equal to the acutal number of bins in "
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
    hist_collectors *= 0
    hist_collectors_std = hist_collectors.copy()
    np.save(save_to + sim_name + '-' + edge_name + '.npy', bin_edges)
    results = {
        'n_bins': n_bins,
        'bin_edges': bin_edges,
        'collector': hist_collectors,
        'collector_std': hist_collectors_std,
        'range': (_lmin, _lmax)
    }
    return results


def sum_rule_bug(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
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
    fsd_t = np.empty(0)
    rflory_t = np.empty(0)
    gyr_t = np.empty(0)
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = np.empty(0)
    shape_parameter_t = np.empty(0)
    for _ in sliced_trj:
        # various measures of chain size
        fsd_t = np.append(fsd_t, np.array([fsd(bug.positions)]), axis=0)
        gyr_t = np.append(gyr_t, np.array([bug.radius_of_gyration()]), axis=0)
        rms = end_to_end(bug.positions)
        rflory_t = np.append(rflory_t, np.array([rms]), axis=0)
        # shape parameters:
        asphericity_t = np.append(
            asphericity_t,
            np.array([bug.asphericity(pbc=False, unwrap=False)]),
            axis=0
        )
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(pbc=False)]),
            axis=0
        )
        shape_parameter_t = np.append(
            shape_parameter_t,
            np.array([bug.shape_parameter(pbc=False)]),
            axis=0
        )

    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def sum_rule_bug_flory_hist(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
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
            'lmin': 0,
            'lmax': sim_info.nmon * sim_info.dmon
            }
        }
    # distribution of the size of the end-to-end
    rflory_hist_info = fixedsize_bins(
        sim_name,
        'rfloryEdgeMon',
        bin_edges['rfloryEdge']['bin_size'],
        bin_edges['rfloryEdge']['lmin'],
        bin_edges['rfloryEdge']['lmax'],
        bin_type='nonnegative',
        save_to=save_to
    )
    fsd_t = np.empty(0)
    rflory_t = np.empty(0)
    gyr_t = np.empty(0)
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = np.empty(0)
    shape_parameter_t = np.empty(0)
    if any([
            rflory_hist_info['collector'].any() != 0,
            rflory_hist_info['collector_std'].any() != 0,
            ]):
        raise ValueError("One of the histogram collectors are not empty!")
    for _ in sliced_trj:
        # various measures of chain size
        fsd_t = np.append(fsd_t, np.array([fsd(bug.positions)]), axis=0)
        gyr_t = np.append(gyr_t, np.array([bug.radius_of_gyration()]), axis=0)
        rms = end_to_end(bug.positions)
        rflory_t = np.append(rflory_t, np.array([rms]), axis=0)
        frame_hist, _ = np.histogram(
            rms,
            bins=rflory_hist_info['bin_edges'],
            range=rflory_hist_info['range']
        )
        # RDF of the end-to-end distance
        rflory_hist_info['collector'] += frame_hist
        rflory_hist_info['collector_std'] += np.square(frame_hist)
        # shape parameters:
        asphericity_t = np.append(
            asphericity_t,
            np.array([bug.asphericity(pbc=False, unwrap=False)]),
            axis=0
        )
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(pbc=False)]),
            axis=0
        )
        shape_parameter_t = np.append(
            shape_parameter_t,
            np.array([bug.shape_parameter(pbc=False)]),
            axis=0
        )

    np.save(
        save_to + sim_name + '-rfloryHistMon.npy',
        rflory_hist_info['collector']
    )
    np.save(
        save_to + sim_name + '-rfloryHistStdMon.npy',
        rflory_hist_info['collector_std']
    )
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


# noinspection PyUnresolvedReferences
def sum_rule_bug_rmsd(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './'
) -> None:
    """Computes the rmsd of a 'segment simulation of a 'bug' atom group in the
    `geometry` of interest, and then saves the output to the `save_to`
    directory.

    `rmsd_bug` does not support the `continuous` option defined in
    `probe_bug` and `probe_all`.

    Parameters
    ----------
    topology: str
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


def sum_rule_all(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest, and saves a variety of
    outputs (mostly in the csv format) to the `save_to` directory.

    Parameters
    ----------
    topology: str
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
        )  # in radians betwene [-np.pi, np.pi]
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

    lastname = 'Crd'
    np.save(
        save_to + sim_name + '-rHist' + lastname + '.npy',
        r_hist_crd_info['collector']
    )
    np.save(
        save_to + sim_name + '-rHistStd' + lastname + '.npy',
        r_hist_crd_info['collector_std']
    )
    np.save(
        save_to + sim_name + '-zHist' + lastname + '.npy',
        z_hist_crd_info['collector']
    )
    np.save(
        save_to + sim_name + '-zHistStd' + lastname + '.npy',
        z_hist_crd_info['collector_std']
    )
    np.save(
        save_to + sim_name + '-thetaHist' + lastname + '.npy',
        theta_hist_crd_info['collector']
    )
    np.save(
        save_to + sim_name + '-thetaHistStd' + lastname + '.npy',
        theta_hist_crd_info['collector_std']
    )
    lastname = 'Mon'
    np.save(
        save_to + sim_name + '-rHist' + lastname + '.npy',
        r_hist_mon_info['collector']
    )
    np.save(
        save_to + sim_name + '-rHistStd' + lastname + '.npy',
        r_hist_mon_info['collector_std']
    )
    np.save(
        save_to + sim_name + '-zHist' + lastname + '.npy',
        z_hist_mon_info['collector']
    )
    np.save(
        save_to + sim_name + '-zHistStd' + lastname + '.npy',
        z_hist_mon_info['collector_std']
    )
    np.save(
        save_to + sim_name + '-thetaHist' + lastname + '.npy',
        theta_hist_mon_info['collector']
    )
    np.save(
        save_to + sim_name + '-thetaHistStd' + lastname + '.npy',
        theta_hist_mon_info['collector_std']
    )
    print('done.')


def trans_fuci_bug(
    topology: str,
    trajectory: str,
    geometry: str = 'biaxial',
    lineage: str = 'segment',
    save_to: str = './',
    continuous: bool = False
) -> None:
    """Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Parameters
    ----------
    topology: str
        Name of the topology file.
    trajectory: str
        Name of the trajectory file.
    geometry : {'biaxial', 'slit', 'box'}, default 'biaxial'
        Shape of the simulation box.
    lineage: {'segment', 'whole'}, default 'segment'
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
    sim_info = TransFoci(
        trajectory,
        geometry=geometry,
        group='bug',
        lineage=lineage
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    lj_nstep = sim_info.bdump  # Sampling steps via dump command in Lammps
    lj_dt = sim_info.dt
    sim_real_dt = lj_nstep * lj_dt * time_unit
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
    bug = cell.select_atoms('resid 1')  # the bug/polymer
    # defining collectors
    foci_pairs = list(combinations(range(len(foci)), 2))  # number of foci
    foci_t = np.empty([0, 5])
    fsd_t = np.empty(0)
    rflory_t = np.empty(0)
    gyr_t = np.empty(0)
    principal_axes_t = np.empty([0, 3, 3])
    asphericity_t = np.empty(0)
    shape_parameter_t = np.empty(0)
    for _ in sliced_trj:
        # various measures of chain size
        fsd_t = np.append(fsd_t, np.array([fsd(bug.positions)]), axis=0)
        gyr_t = np.append(gyr_t, np.array([bug.radius_of_gyration()]), axis=0)
        rms = end_to_end(bug.positions)
        rflory_t = np.append(rflory_t, np.array([rms]), axis=0)
        # shape parameters:
        asphericity_t = np.append(
            asphericity_t,
            np.array([bug.asphericity(pbc=False, unwrap=False)]),
            axis=0
        )
        principal_axes_t = np.append(
            principal_axes_t,
            np.array([bug.principal_axes(pbc=False)]),
            axis=0
        )
        shape_parameter_t = np.append(
            shape_parameter_t,
            np.array([bug.shape_parameter(pbc=False)]),
            axis=0
        )
        for (i, j) in foci_pairs:
            pair_dist = np.array([
                i,
                j,
                foci.atoms[i].id,
                foci.atoms[j].id,
                np.linalg.norm(
                    foci.atoms[i].position - foci.atoms[j].position)]
            )
            foci_t = np.append(foci_t, np.array([pair_dist]), axis=0)

    np.save(save_to + sim_name + '-distTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')
