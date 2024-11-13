import warnings
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any, Tuple
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mda_dist
from polyphys.manage.typer import (
    ParserInstance,
    SumRuleCylInstance,
    TwoMonDepInstance)
from polyphys.manage.organizer import invalid_keyword
from polyphys.analyze import clusters, correlations
from polyphys.analyze.measurer import transverse_size, fsd, end_to_end


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
    bin_type: str = "ordinary",
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
    bin_type: {"ordinary", "nonnegative", "periodic"}, default "ordinary"
        The type of bin in a given direction in a given coordinate system:

        "ordinary"
            A bounded or unbounded coordinate such as any of the cartesian
            coordinates or the polar coordinate in the spherical coordinate
            system. For such coordinates, the `lmin` and `lmax` limits are
            equally extended to ensure the `bin_size`.

        "nonnegative"
            A nonnegative coordinate such as the r direction in the polar
            or spherical coordinate system. For such coordinates, ONLY `lmax`
            limit is extended to ensure the `bin_size`. `lmin` is either 0.0
            or a positive number smaller than `lmax`.

        "periodic"
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
    bin_types = ["ordinary", 'nonnagative', "periodic"]
    if lmin >= lmax:
        raise ValueError('Boundaries are not sane: should be xmin < xmax.')
    _delta = bin_size
    _lmin = lmin
    _lmax = lmax
    _length = _lmax - _lmin
    n_bins: int = 0
    if bin_type == "ordinary":
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
    elif bin_type == "nonnegative":
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
    elif bin_type == "periodic":  # Assuming that the _length=period:
        n_bins = int(np.ceil(_length / _delta))
        warnings.warn(
            f"Number of bins (n_bins='{n_bins}')"
            " is more than or equal to the actual number of bins in "
            f""periodic" bin type because the 'period=lmax-min={_length}'"
            f"and delta='{_delta}'"
            ",not 'n_bins', are used to created 'bin_edges'.",
            UserWarning
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


def write_hists(
    hist_infos: Dict[str, Any],
    sim_name: str,
    save_to: str,
    std: bool = False
) -> None:
    """
    Writes histogram per entity per direction to file.

    Parameters
    ----------
    hist_infos : Dict[str, Any]
        A dict of dicts that contains the information about direction, entities,
         and histograms.
    sim_name: str
        The name of simulation file to which the `hist_infos` belongs.
    save_to: str
        The absolute/relative path of a directory to which outputs are saved.
    std : bool, default False
        _description_, by default False
    """
    for dir_ in hist_infos.keys():
        for entity, hist in hist_infos[dir_].items():
            np.save(
                save_to + sim_name + '-' + dir_ + entity + '.npy',
                hist['collector']
            )
            if std is True:
                np.save(
                    save_to + sim_name + '-' + dir_ + 'Std' + entity + '.npy',
                    hist['collector_std']
                )
        # end of loop
    # end of loop


class ProberBase(ABC):
    """_summary_

    Parameters
    ----------
    ABC : _type_
        _description_
    """
    def __init__(
        self,
        topology: str,
        trajectory: str,
        parser: ParserInstance,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        if (self._parser.lineage == 'segment') & (continuous is False):
            warnings.warn(
                "lineage is "
                f"'{self.parser.lineage}' "
                "and 'continuous' is "
                f"'{continuous}. "
                "Please ensure the "
                f"'{trajectory}' is NOT part of a sequence of trajectories.",
                RuntimeWarning)
        self._topology = topology
        self._trajectory = trajectory
        self._parser = parser
        self._continuous = continuous
        self._save_to = save_to
        self._sim_name: str = f"{self._parser.name}-{self._parser.group}"
        self._n_frames, self._sliced_trj = self._set_trj_slice()
        self._atom_groups: Dict[str, mda.AtomGroup] = {}
        self._define_atom_groups()
        self._collectors: Dict[str, Any] = {}
        self._bin_edges: Dict[str, Dict[str, float]] = {}
        self._setup_collectors()

    @property
    def topology(self) -> str:
        """
        Returns the topology file associated with the molecular dynamics
        system.
        """
        return self._topology

    @property
    def trajectory(self) -> str:
        """
        Returns the trajectory file associated with the molecular dynamics
        system.
        """
        return self._trajectory

    @property
    def continuous(self) -> bool:
        """
        Shows whether the trajectory is whole or segment.
        """
        return self._continuous

    @property
    def parser(self) -> ParserInstance:
        """
        Returns a parser object, containing all the information a molecular
        dynamics system, extracted from `trjectory` filename.
        """
        return self._parser

    @property
    def save_to(self) -> str:
        """
        Returns the directory to which the artifacts, i.e. the files associated
        with the `results` of probing a molecular dynamics system, are written.
        """
        return self._save_to

    @property
    def sim_name(self) -> str:
        """
        Returns the simulation name.
        """
        return self._sim_name

    @property
    def sliced_trj(self) -> mda.coordinates.base.FrameIteratorSliced:
        """
        Returns the sliced trajectory over which probing is performed.
        """

    @property
    def n_frames(self) -> int:
        """
        Returns the total number of time frames in a molecular dynamics system.
        """

    @property
    def atom_groups(self) -> Dict[str, mda.AtomGroup]:
        """
        Returns the dict of atom groups and their associated MDAnalysis
        atomgroup instances.
        """
        return self._atom_groups

    @property
    def bin_edges(self) -> Dict[str, Dict[str, float]]:
        """
        Returns the dictionary of directions, bin sizes, and limits for bin
        edges.
        """
        return self._bin_edges

    @property
    @abstractmethod
    def _universe(self) -> mda.Universe:
        """
        an MDAnalysis Universe object, containing all the information
        describing a molecular dynamics system.
        """

    @property
    def universe(self) -> mda.Universe:
        """
        Returns an MDAnalysis Universe object, containing all the information
        describing a molecular dynamics system.
        """
        return self._universe

    @property
    @abstractmethod
    def _damping_time(self) -> float:
        """
        The time difference between two consecutive time frames/snapshots in
        the trajectory file of a molecular dynamics system.
        """

    @property
    def damping_time(self) -> float:
        """
        Returns the time unit of a molecular dynamics system.
        """
        return self._damping_time

    @property
    def collectors(self) -> Dict[str, np.ndarray]:
        """
        Returns a list of physical properties, i.e. physical measurement,
        extracted from a molecular dynamics trajectory by the prober.
        """
        return self._collectors

    def _set_trj_slice(
        self
    ) -> Tuple[int, mda.coordinates.base.FrameIteratorSliced]:
        if self._continuous is True:
            return self._universe.trajectory.n_frames - 1, \
                self._universe.trajectory[0: -1]
        return self._universe.trajectory.n_frames, self._universe.trajectory

    def simulation_report(self) -> None:
        """
        Writes a moecular dynamics simulation details to a file.

        `stamps_report` generates a dataset called `report_name` of the
        values of some attributes of `sim_info`, the number of frames
        `n_frames`, all the key and value pairs in all the given dictionaries
        `measures`.
        """
        report_name = f"{self._sim_name}-stamps.csv"
        with open(report_name, mode="w", encoding="utf-8") as report:
            # write header
            for attr in self._parser.attributes:
                report.write(f"{attr},")
            report.write("n_frames\n")
            # write values
            for lineage_name in self._parser.attributes:
                attr_value = getattr(self._parser, lineage_name)
                report.write(f"{attr_value},")
            report.write(f"{self.n_frames}")
        print("Simulation report written.")

    @abstractmethod
    def _define_atom_groups(self) -> None:
        """
        Defines atom groups in a molecular dynamics project.
        """

    @abstractmethod
    def _setup_collectors(self) -> None:
        """
        Defines physical properties and initialize their defaults values.
        """

    @abstractmethod
    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single frame for specific atom group and update collectors.

        Parameters
        ----------
        idx: int
            Index a trajectory frame (snapashot).
        """

    def run(self) -> None:
        """
        Probes all the frames selected via 

        Parameters
        ----------
        atom_groups : dict
            _description_
        """
        print(f"Probing '{self._parser.name}' ...")
        for idx, _ in enumerate(self._sliced_trj):
            self._probe_frame(idx)
        print("done.")

    def save_artifacts(self):
        """
        Save all relevant analysis data for the current simulation
        """
        for prop, artifact in self._collectors.items():
            filename = f"{self._save_to}{self._sim_name}-{prop}.npy"
            np.save(filename, artifact)
        print("Artifacts saved.")


class TwoMonDepCubBugProber(ProberBase):
    """_summary_

    Parameters
    ----------
    ProberBase : _type_
        _description_
    """
    def __init__(
        self,
        topology: str,
        trajectory: str,
        parser: TwoMonDepInstance,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, parser, save_to,
                         continuous=continuous)
        self._damping_time = \
            getattr(self._parser, "bdump") * getattr(self._parser, "dt")
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id type x y z",
            dt=self._damping_time
        )

    def _define_atom_groups(self) -> None:
        # the two monomers:
        self._atom_groups['bug'] = self._universe.select_atoms('type 1')

    def _setup_collectors(self) -> None:
        self._collectors = {
            "gyrTMon": np.zeros(self._n_frames),
            "dxTMon": np.zeros(self._n_frames),
            "dyTMon": np.zeros(self._n_frames),
            "dzTMon": np.zeros(self._n_frames)
        }

    def _probe_frame(self, idx) -> None:
        self._collectors["gyrTMon"][idx] = \
            self._atom_groups['bug'].radius_of_gyration()
        self._collectors["dxTMon"][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=0)
        self._collectors["dyTMon"][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=1)
        self._collectors["dzTMon"][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=2)


class SumRuleCylBugProber(ProberBase):
    """
    Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, LAMMPS recenter is used to restrict the center of
    mass of "bug" (monomers) to the center of simulation box; and consequently,
    coordinates of all the particles in a trajectory or topology file is
    recentered to fulfill this constraint.

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    def __init__(
        self,
        topology: str,
        trajectory: str,
        parser: SumRuleCylInstance,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, parser, save_to,
                         continuous=continuous)
        self._damping_time = \
            getattr(self._parser, "bdump") * getattr(self._parser, "dt")
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self._damping_time
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['bug'] = self._universe.select_atoms('resid 1')

    def _setup_collectors(self) -> None:
        self._collectors = {
            "transSizeTMon": np.zeros(self._n_frames),
            "fsdTMon": np.zeros(self._n_frames),
            "gyrTMon": np.zeros(self._n_frames),
            "rfloryTMon": np.zeros(self._n_frames),
            "asphericityTMon": np.zeros(self._n_frames),
            "shapeTMon": np.zeros(self._n_frames),
            "principalTMon": np.zeros([self._n_frames, 3, 3])
        }

    def _probe_frame(self, idx) -> None:
        self._collectors["transSizeTMon"][idx] = \
            transverse_size(self._atom_groups['bug'].positions, axis=2)
        self._collectors["fsdTMon"][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=2)
        self._collectors["gyrTMon"][idx] = \
            self._atom_groups['bug'].radius_of_gyration()
        self._collectors["rfloryTMon"][idx] = \
            end_to_end(self._atom_groups['bug'].positions)
        self._collectors["asphericityTMon"][idx] = \
            self._atom_groups['bug'].asphericity(wrap=False, unwrap=False)
        self._collectors["shapeTMon"][idx] = \
            self._atom_groups['bug'].shape_parameter(wrap=False)
        self._collectors["principalTMon"][idx] = \
            self._atom_groups['bug'].principal_axes(wrap=False)


class SumRuleCylAllProber(ProberBase):
    """
    Runs various analyses on a `lineage` simulation of a 'bug' atom group in
    the `geometry` of interest.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, LAMMPS recenter is used to restrict the center of
    mass of "bug" (monomers) to the center of simulation box; and consequently,
    coordinates of all the particles in a trajectory or topology file is
    recentered to fulfill this constraint.

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    def __init__(
        self,
        topology: str,
        trajectory: str,
        parser: ParserInstance,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, parser, save_to,
                         continuous=continuous)
        self._damping_time = \
            getattr(self._parser, "adump") * getattr(self._parser, "dt")
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self._damping_time
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['crds'] = self._universe.select_atoms('resid 0')
        self._atom_groups['bug'] = self._universe.select_atoms('resid 1')

    def _setup_bins(self) -> None:
        self._bin_edges = {
            "rEdge": {
                "bin_size":  0.1 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                "lmin": 0,
                "lmax": 0.5 * getattr(self._parser, 'dcyl'),
                },
            "zEdge": {
                "bin_size":  0.5 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                "lmin": -0.5 * getattr(self._parser, 'lcyl'),
                "lmax": 0.5 * getattr(self._parser, 'lcyl'),
                },
            "thetaEdge": {
                "bin_size":  np.pi / 36,
                "lmin": -1 * np.pi,
                "lmax": np.pi
                },
            "lEdge": {
                # edges for 2d hist in x and y direction
                "bin_size":  0.1 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                "lmin": -0.5 * getattr(self._parser, 'dcyl'),
                "lmax": 0.5 * getattr(self._parser, 'dcyl')
                }
        }
        entities = ["Mon", "Crd"]
        hist1d_bin_type = \
            {"r": "nonnegative", "z": "ordinary", "theta": "periodic"}
        hist2d_planes_dirs = ["xy": "z", "xz": "y", "yz": "x"]
        for ent in entities:
            for dir, edge_t in hist1d_bin_type.items():
                collector_info = fixedsize_bins(
                    self._sim_name,
                    f"{dir}Edge",
                    self._bin_edges[f"{dir}Edge"]["bin_size"],
                    self._bin_edges[f"{dir}Edge"]["lmin"],
                    self._bin_edges[f"{dir}Edge"]["lmax"],
                    bin_type=edge_t,
                    save_to=self._save_to
                )
                self._collectors[f"{dir}Hist{ent}"] = \
                    np.zeros(collector_info['n_bins'])
            for plane, dir in hist2d_planes_dirs.items():
                collector_info = fixedsize_bins(
                    self._sim_name,
                    f"{dir}Edge",
                    self._bin_edges[f"lEdge"]["bin_size"],
                    self._bin_edges[f"lEdge"]["lmin"],
                    self._bin_edges[f"lEdge"]["lmax"],
                    bin_type='ordinary',
                    save_to=self._save_to
                )
                self._collectors[f"{plane}Hist{ent}"] = 0


    def _probe_frame(self, idx) -> None:
        self._collectors["transSizeTMon"][idx] = \
            transverse_size(self._atom_groups['bug'].positions, axis=2)
        self._collectors["fsdTMon"][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=2)
        self._collectors["gyrTMon"][idx] = \
            self._atom_groups['bug'].radius_of_gyration()
        self._collectors["rfloryTMon"][idx] = \
            end_to_end(self._atom_groups['bug'].positions)
        self._collectors["asphericityTMon"][idx] = \
            self._atom_groups['bug'].asphericity(wrap=False, unwrap=False)
        self._collectors["shapeTMon"][idx] = \
            self._atom_groups['bug'].shape_parameter(wrap=False)
        self._collectors["principalTMon"][idx] = \
            self._atom_groups['bug'].principal_axes(wrap=False)

def sum_rule_cyl_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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

    # radial direction of the cylindrical coordinate system
    # longitudinal direction of the cylindrical coordinate system
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    # y direction of the cartesian coordinate system
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
            theta_hist_crd_info['collector_std'].any() != 0,
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    for _ in sliced_trj:
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions[:, :2], axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
        )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            crds.positions[:, 2],
            bins=z_hist_crd_info['bin_edges'],
            range=z_hist_crd_info['range']
        )
        z_hist_crd_info['collector'] += pos_hist
        z_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # theta in radians
        pos_hist, _ = np.histogram(
            np.arctan2(crds.positions[:, 1], crds.positions[:, 0]),
            bins=theta_hist_crd_info['bin_edges'],
            range=theta_hist_crd_info['range']
        )
        theta_hist_crd_info['collector'] += pos_hist
        theta_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions[:, :2], axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
        )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            bug.positions[:, 2],
            bins=z_hist_mon_info['bin_edges'],
            range=z_hist_mon_info['range']
        )
        z_hist_mon_info['collector'] += pos_hist
        z_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # theta in radian
        pos_hist, _ = np.histogram(
            np.arctan2(bug.positions[:, 1], bug.positions[:, 0]),
            bins=theta_hist_mon_info['bin_edges'],
            range=theta_hist_mon_info['range']
        )
        theta_hist_mon_info['collector'] += pos_hist
        theta_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def trans_foci_cyl_bug(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # foci
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # bug: small & large mon
    # defining collectors
    # -bug:
    trans_size_t = np.zeros(n_frames)
    fsd_t = np.zeros(n_frames)
    gyr_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # -foci:
    foci_t = np.zeros((n_frames, sim_info.nmon_large, sim_info.nmon_large))
    dir_contacts_t = np.zeros(
        (n_frames, sim_info.nmon_large, sim_info.nmon_large), dtype=int
    )
    bonds_t = np.zeros((n_frames, sim_info.nmon_large), dtype=int)
    clusters_t = np.zeros((n_frames, sim_info.nmon_large+1), dtype=int)
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain size
        trans_size_t[idx] = transverse_size(bug)
        fsd_t[idx] = fsd(bug.positions)
        gyr_t[idx] = bug.radius_of_gyration()
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=False)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        # foci:
        dist_mat = clusters.self_dist_array(foci.positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t[idx] = foci_pair_dist
        dir_contacts = clusters.find_direct_contacts(dist_mat, cluster_cutoff)
        dir_contacts_t[idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        bonds_t[idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        clusters_t[idx] = clusters_stat
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-transSizeTMon.npy', trans_size_t)
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def trans_foci_cyl_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
        "rEdge": {
            "bin_size":  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.dcyl
            },
        "zEdge": {
            "bin_size":  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcyl,
            "lmax": 0.5 * sim_info.lcyl
            },
        "thetaEdge": {
            "bin_size":  np.pi / 36,
            "lmin": -1 * np.pi,
            "lmax": np.pi
            },
        "lEdge": {
            # edges for 2d hist in x and y direction
            "bin_size":  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.dcyl,
            "lmax": 0.5 * sim_info.dcyl
            }
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # polymer/monomers
    dna: mda.AtomGroup = cell.select_atoms('type 1')  # small monomers
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # large monomers
    # bin edges and histograms in different directions:
    # radial direction of the cylindrical coordinate system
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # z direction of the cylindrical coordinate system
    z_hist_crd_info = fixedsize_bins(
        sim_name,
        'zEdgeCrd',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    z_hist_mon_info = fixedsize_bins(
        sim_name,
        'zEdgeMon',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    z_hist_foci_info = fixedsize_bins(
        sim_name,
        'zEdgeFoci',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    z_hist_dna_info = fixedsize_bins(
        sim_name,
        'zEdgeDna',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # theta of the cylindrical coordinate system
    theta_hist_crd_info = fixedsize_bins(
        sim_name,
        'thetaEdgeCrd',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    theta_hist_mon_info = fixedsize_bins(
        sim_name,
        'thetaEdgeMon',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    theta_hist_foci_info = fixedsize_bins(
        sim_name,
        'thetaEdgeFoci',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    theta_hist_dna_info = fixedsize_bins(
        sim_name,
        'thetaEdgeDna',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
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
            theta_hist_mon_info['collector_std'].any() != 0,
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # foci
    # # xy
    xy_hist_foci_info = {
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
    xy_hist_foci_info['collector'] = np.zeros(xy_hist_foci_info['n_bins'])
    xy_hist_foci_info['collector'] *= 0
    # # xz
    xz_hist_foci_info = {
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
    xz_hist_foci_info['collector'] = np.zeros(xz_hist_foci_info['n_bins'])
    xz_hist_foci_info['collector'] *= 0
    # # yz
    yz_hist_foci_info = {
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
    yz_hist_foci_info['collector'] = np.zeros(yz_hist_foci_info['n_bins'])
    yz_hist_foci_info['collector'] *= 0
    # dna
    # # xy
    xy_hist_dna_info = {
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
    xy_hist_dna_info['collector'] = np.zeros(xy_hist_dna_info['n_bins'])
    xy_hist_dna_info['collector'] *= 0
    # # xz
    xz_hist_dna_info = {
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
    xz_hist_dna_info['collector'] = np.zeros(xz_hist_dna_info['n_bins'])
    xz_hist_dna_info['collector'] *= 0
    # # yz
    yz_hist_dna_info = {
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
    yz_hist_dna_info['collector'] = np.zeros(yz_hist_dna_info['n_bins'])
    yz_hist_dna_info['collector'] *= 0
    for _ in sliced_trj:
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions[:, :2], axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
        )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            crds.positions[:, 2],
            bins=z_hist_crd_info['bin_edges'],
            range=z_hist_crd_info['range']
        )
        z_hist_crd_info['collector'] += pos_hist
        z_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # theta in radians
        pos_hist, _ = np.histogram(
            np.arctan2(crds.positions[:, 1], crds.positions[:, 0]),
            bins=theta_hist_crd_info['bin_edges'],
            range=theta_hist_crd_info['range']
        )
        theta_hist_crd_info['collector'] += pos_hist
        theta_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions[:, :2], axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
        )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            bug.positions[:, 2],
            bins=z_hist_mon_info['bin_edges'],
            range=z_hist_mon_info['range']
        )
        z_hist_mon_info['collector'] += pos_hist
        z_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # theta in radian
        pos_hist, _ = np.histogram(
            np.arctan2(bug.positions[:, 1], bug.positions[:, 0]),
            bins=theta_hist_mon_info['bin_edges'],
            range=theta_hist_mon_info['range']
        )
        theta_hist_mon_info['collector'] += pos_hist
        theta_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # foci
        pos_hist, _ = np.histogram(
            np.linalg.norm(foci.positions[:, :2], axis=1),
            bins=r_hist_foci_info['bin_edges'],
            range=r_hist_foci_info['range']
        )
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            foci.positions[:, 2],
            bins=z_hist_foci_info['bin_edges'],
            range=z_hist_foci_info['range']
        )
        z_hist_foci_info['collector'] += pos_hist
        z_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # theta in radians
        pos_hist, _ = np.histogram(
            np.arctan2(foci.positions[:, 1], foci.positions[:, 0]),
            bins=theta_hist_foci_info['bin_edges'],
            range=theta_hist_foci_info['range']
        )
        theta_hist_foci_info['collector'] += pos_hist
        theta_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=xy_hist_foci_info['bin_edges'],
            range=xy_hist_foci_info['range'],
        )
        xy_hist_foci_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=xz_hist_foci_info['bin_edges'],
            range=xz_hist_foci_info['range'],
        )
        xz_hist_foci_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=yz_hist_foci_info['bin_edges'],
            range=yz_hist_foci_info['range'],
        )
        yz_hist_foci_info['collector'] += pos_hist
        # dna
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(dna.positions[:, :2], axis=1),
            bins=r_hist_dna_info['bin_edges'],
            range=r_hist_dna_info['range']
        )
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            dna.positions[:, 2],
            bins=z_hist_dna_info['bin_edges'],
            range=z_hist_dna_info['range']
        )
        z_hist_dna_info['collector'] += pos_hist
        z_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # theta in radians
        pos_hist, _ = np.histogram(
            np.arctan2(dna.positions[:, 1], dna.positions[:, 0]),
            bins=theta_hist_dna_info['bin_edges'],
            range=theta_hist_dna_info['range']
        )
        theta_hist_dna_info['collector'] += pos_hist
        theta_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=xy_hist_dna_info['bin_edges'],
            range=xy_hist_dna_info['range'],
        )
        xy_hist_dna_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=xz_hist_dna_info['bin_edges'],
            range=xz_hist_dna_info['range'],
        )
        xz_hist_dna_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=yz_hist_dna_info['bin_edges'],
            range=yz_hist_dna_info['range'],
        )
        yz_hist_dna_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Foci': xy_hist_foci_info,
            'Dna': xy_hist_dna_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Foci': xz_hist_foci_info,
            'Dna': xz_hist_dna_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Foci': yz_hist_foci_info,
            'Dna': yz_hist_dna_info
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def trans_foci_cub_bug(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # foci: large monomers
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # bug: small & large mon
    # defining collectors
    # -bug:
    gyr_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # -foci:
    foci_t = np.zeros((n_frames, sim_info.nmon_large, sim_info.nmon_large))
    dir_contacts_t = np.zeros(
        (n_frames, sim_info.nmon_large, sim_info.nmon_large), dtype=int
    )
    bonds_t = np.zeros((n_frames, sim_info.nmon_large), dtype=int)
    clusters_t = np.zeros((n_frames, sim_info.nmon_large+1), dtype=int)
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain size
        gyr_t[idx] = bug.radius_of_gyration()
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=True)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        # foci:
        dist_mat = clusters.self_dist_array(foci.positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t[idx] = foci_pair_dist
        dir_contacts = clusters.find_direct_contacts(dist_mat, cluster_cutoff)
        dir_contacts_t[idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        bonds_t[idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        clusters_t[idx] = clusters_stat
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def trans_foci_cub_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
        "rEdge": {  # edges for distance r from the box center
            "bin_size":  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.lcube
            },
        "lEdge": {  # edges in cartesian coordinates
            "bin_size":  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcube,
            "lmax": 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    dna: mda.AtomGroup = cell.select_atoms('type 1')  # small monomers
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # large monomers
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # foci
    # # xy
    xy_hist_foci_info = {
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
    xy_hist_foci_info['collector'] = np.zeros(xy_hist_foci_info['n_bins'])
    xy_hist_foci_info['collector'] *= 0
    # # xz
    xz_hist_foci_info = {
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
    xz_hist_foci_info['collector'] = np.zeros(xz_hist_foci_info['n_bins'])
    xz_hist_foci_info['collector'] *= 0
    # # yz
    yz_hist_foci_info = {
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
    yz_hist_foci_info['collector'] = np.zeros(yz_hist_foci_info['n_bins'])
    yz_hist_foci_info['collector'] *= 0
    # dna
    # # xy
    xy_hist_dna_info = {
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
    xy_hist_dna_info['collector'] = np.zeros(xy_hist_dna_info['n_bins'])
    xy_hist_dna_info['collector'] *= 0
    # # xz
    xz_hist_dna_info = {
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
    xz_hist_dna_info['collector'] = np.zeros(xz_hist_dna_info['n_bins'])
    xz_hist_dna_info['collector'] *= 0
    # # yz
    yz_hist_dna_info = {
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
    yz_hist_dna_info['collector'] = np.zeros(yz_hist_dna_info['n_bins'])
    yz_hist_dna_info['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions, axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
            )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions, axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
            )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # foci
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(foci.positions, axis=1),
            bins=r_hist_foci_info['bin_edges'],
            range=r_hist_foci_info['range']
            )
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=xy_hist_foci_info['bin_edges'],
            range=xy_hist_foci_info['range'],
        )
        xy_hist_foci_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=xz_hist_foci_info['bin_edges'],
            range=xz_hist_foci_info['range'],
        )
        xz_hist_foci_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=yz_hist_foci_info['bin_edges'],
            range=yz_hist_foci_info['range'],
        )
        yz_hist_foci_info['collector'] += pos_hist
        # dna
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(dna.positions, axis=1),
            bins=r_hist_dna_info['bin_edges'],
            range=r_hist_dna_info['range']
            )
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=xy_hist_dna_info['bin_edges'],
            range=xy_hist_dna_info['range'],
        )
        xy_hist_dna_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=xz_hist_dna_info['bin_edges'],
            range=xz_hist_dna_info['range'],
        )
        xz_hist_dna_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=yz_hist_dna_info['bin_edges'],
            range=yz_hist_dna_info['range'],
        )
        yz_hist_dna_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Foci': xy_hist_foci_info,
            'Dna': xy_hist_dna_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Foci': xz_hist_foci_info,
            'Dna': xz_hist_dna_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Foci': yz_hist_foci_info,
            'Dna': yz_hist_dna_info
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def sum_rule_hetero_ring_cub_bug(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = SumRuleCubHeteroRing(
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
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # foci: large monomers
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # bug: small & large mon
    # defining collectors
    # -bug:
    gyr_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # -foci:
    foci_t = np.zeros((n_frames, sim_info.nmon_large, sim_info.nmon_large))
    dir_contacts_t = np.zeros(
        (n_frames, sim_info.nmon_large, sim_info.nmon_large), dtype=int
    )
    bonds_t = np.zeros((n_frames, sim_info.nmon_large), dtype=int)
    clusters_t = np.zeros((n_frames, sim_info.nmon_large+1), dtype=int)
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain size
        gyr_t[idx] = bug.radius_of_gyration()
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=True)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        # foci:
        dist_mat = clusters.self_dist_array(foci.positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t[idx] = foci_pair_dist
        dir_contacts = clusters.find_direct_contacts(dist_mat, cluster_cutoff)
        dir_contacts_t[idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        bonds_t[idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        clusters_t[idx] = clusters_stat
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def sum_rule_hetero_ring_cub_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = SumRuleCubHeteroRing(
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
        "rEdge": {  # edges for distance r from the box center
            "bin_size":  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.lcube
            },
        "lEdge": {  # edges in cartesian coordinates
            "bin_size":  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcube,
            "lmax": 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    dna: mda.AtomGroup = cell.select_atoms('type 1')  # small monomers
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # large monomers
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # foci
    # # xy
    xy_hist_foci_info = {
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
    xy_hist_foci_info['collector'] = np.zeros(xy_hist_foci_info['n_bins'])
    xy_hist_foci_info['collector'] *= 0
    # # xz
    xz_hist_foci_info = {
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
    xz_hist_foci_info['collector'] = np.zeros(xz_hist_foci_info['n_bins'])
    xz_hist_foci_info['collector'] *= 0
    # # yz
    yz_hist_foci_info = {
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
    yz_hist_foci_info['collector'] = np.zeros(yz_hist_foci_info['n_bins'])
    yz_hist_foci_info['collector'] *= 0
    # dna
    # # xy
    xy_hist_dna_info = {
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
    xy_hist_dna_info['collector'] = np.zeros(xy_hist_dna_info['n_bins'])
    xy_hist_dna_info['collector'] *= 0
    # # xz
    xz_hist_dna_info = {
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
    xz_hist_dna_info['collector'] = np.zeros(xz_hist_dna_info['n_bins'])
    xz_hist_dna_info['collector'] *= 0
    # # yz
    yz_hist_dna_info = {
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
    yz_hist_dna_info['collector'] = np.zeros(yz_hist_dna_info['n_bins'])
    yz_hist_dna_info['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions, axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
            )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions, axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
            )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # foci
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(foci.positions, axis=1),
            bins=r_hist_foci_info['bin_edges'],
            range=r_hist_foci_info['range']
            )
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=xy_hist_foci_info['bin_edges'],
            range=xy_hist_foci_info['range'],
        )
        xy_hist_foci_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=xz_hist_foci_info['bin_edges'],
            range=xz_hist_foci_info['range'],
        )
        xz_hist_foci_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=yz_hist_foci_info['bin_edges'],
            range=yz_hist_foci_info['range'],
        )
        yz_hist_foci_info['collector'] += pos_hist
        # dna
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(dna.positions, axis=1),
            bins=r_hist_dna_info['bin_edges'],
            range=r_hist_dna_info['range']
            )
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=xy_hist_dna_info['bin_edges'],
            range=xy_hist_dna_info['range'],
        )
        xy_hist_dna_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=xz_hist_dna_info['bin_edges'],
            range=xz_hist_dna_info['range'],
        )
        xz_hist_dna_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=yz_hist_dna_info['bin_edges'],
            range=yz_hist_dna_info['range'],
        )
        yz_hist_dna_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Foci': xy_hist_foci_info,
            'Dna': xy_hist_dna_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Foci': xz_hist_foci_info,
            'Dna': xz_hist_dna_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Foci': yz_hist_foci_info,
            'Dna': yz_hist_dna_info
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def sum_rule_hetero_linear_cub_bug(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = SumRuleCubHeteroLinear(
        trajectory,
        lineage,
        'cubic',
        'bug',
        'linear'
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
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # foci: large monomers
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # bug: small & large mon
    # defining collectors
    # -bug:
    gyr_t = np.zeros(n_frames)
    rflory_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # -foci:
    foci_t = np.zeros((n_frames, sim_info.nmon_large, sim_info.nmon_large))
    dir_contacts_t = np.zeros(
        (n_frames, sim_info.nmon_large, sim_info.nmon_large), dtype=int
    )
    bonds_t = np.zeros((n_frames, sim_info.nmon_large), dtype=int)
    clusters_t = np.zeros((n_frames, sim_info.nmon_large+1), dtype=int)
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain size
        gyr_t[idx] = bug.radius_of_gyration()
        rflory_t[idx] = end_to_end(bug.positions)
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=True)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        # foci:
        dist_mat = clusters.self_dist_array(foci.positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        # keep atom ids on the diag
        np.fill_diagonal(foci_pair_dist, foci.atoms.ids)
        foci_t[idx] = foci_pair_dist
        dir_contacts = clusters.find_direct_contacts(dist_mat, cluster_cutoff)
        dir_contacts_t[idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        bonds_t[idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        clusters_t[idx] = clusters_stat
    # Saving collectors to memory
    # -bug
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-rfloryTMon.npy', rflory_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    # -foci
    np.save(save_to + sim_name + '-distMatTFoci.npy', foci_t)
    np.save(save_to + sim_name + '-directContactsMatTFoci.npy', dir_contacts_t)
    np.save(save_to + sim_name + '-bondsHistTFoci.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistTFoci.npy', clusters_t)
    # Simulation stamps:
    outfile = save_to + sim_name + "-stamps.csv"
    stamps_report(outfile, sim_info, n_frames)
    print('done.')


def sum_rule_hetero_linear_cub_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = SumRuleCubHeteroLinear(
        trajectory,
        lineage,
        'cubic',
        'all',
        'linear'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        "rEdge": {  # edges for distance r from the box center
            "bin_size":  0.1 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.lcube
            },
        "lEdge": {  # edges in cartesian coordinates
            "bin_size":  0.5 * min(sim_info.dmon_small, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcube,
            "lmax": 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon_small * np.sqrt(
        sim_info.mmon_small * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    dna: mda.AtomGroup = cell.select_atoms('type 1')  # small monomers
    foci: mda.AtomGroup = cell.select_atoms('type 2')  # large monomers
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_foci_info = fixedsize_bins(
        sim_name,
        'rEdgeFoci',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_dna_info = fixedsize_bins(
        sim_name,
        'rEdgeDna',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # foci
    # # xy
    xy_hist_foci_info = {
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
    xy_hist_foci_info['collector'] = np.zeros(xy_hist_foci_info['n_bins'])
    xy_hist_foci_info['collector'] *= 0
    # # xz
    xz_hist_foci_info = {
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
    xz_hist_foci_info['collector'] = np.zeros(xz_hist_foci_info['n_bins'])
    xz_hist_foci_info['collector'] *= 0
    # # yz
    yz_hist_foci_info = {
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
    yz_hist_foci_info['collector'] = np.zeros(yz_hist_foci_info['n_bins'])
    yz_hist_foci_info['collector'] *= 0
    # dna
    # # xy
    xy_hist_dna_info = {
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
    xy_hist_dna_info['collector'] = np.zeros(xy_hist_dna_info['n_bins'])
    xy_hist_dna_info['collector'] *= 0
    # # xz
    xz_hist_dna_info = {
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
    xz_hist_dna_info['collector'] = np.zeros(xz_hist_dna_info['n_bins'])
    xz_hist_dna_info['collector'] *= 0
    # # yz
    yz_hist_dna_info = {
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
    yz_hist_dna_info['collector'] = np.zeros(yz_hist_dna_info['n_bins'])
    yz_hist_dna_info['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions, axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
            )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions, axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
            )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # foci
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(foci.positions, axis=1),
            bins=r_hist_foci_info['bin_edges'],
            range=r_hist_foci_info['range']
            )
        r_hist_foci_info['collector'] += pos_hist
        r_hist_foci_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 1],
            bins=xy_hist_foci_info['bin_edges'],
            range=xy_hist_foci_info['range'],
        )
        xy_hist_foci_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 0],
            foci.positions[:, 2],
            bins=xz_hist_foci_info['bin_edges'],
            range=xz_hist_foci_info['range'],
        )
        xz_hist_foci_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            foci.positions[:, 1],
            foci.positions[:, 2],
            bins=yz_hist_foci_info['bin_edges'],
            range=yz_hist_foci_info['range'],
        )
        yz_hist_foci_info['collector'] += pos_hist
        # dna
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(dna.positions, axis=1),
            bins=r_hist_dna_info['bin_edges'],
            range=r_hist_dna_info['range']
            )
        r_hist_dna_info['collector'] += pos_hist
        r_hist_dna_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 1],
            bins=xy_hist_dna_info['bin_edges'],
            range=xy_hist_dna_info['range'],
        )
        xy_hist_dna_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 0],
            dna.positions[:, 2],
            bins=xz_hist_dna_info['bin_edges'],
            range=xz_hist_dna_info['range'],
        )
        xz_hist_dna_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            dna.positions[:, 1],
            dna.positions[:, 2],
            bins=yz_hist_dna_info['bin_edges'],
            range=yz_hist_dna_info['range'],
        )
        yz_hist_dna_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Foci': xy_hist_foci_info,
            'Dna': xy_hist_dna_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Foci': xz_hist_foci_info,
            'Dna': xz_hist_dna_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Foci': yz_hist_foci_info,
            'Dna': yz_hist_dna_info
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
    polymer (called "bug") and several H-NS proteins (each formed of a core and
    two poles or patches) crowded by soft LJ spheres (crowders) in free
    ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to
    the center of mass of the single polymer (bug).

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups:
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # the bug
    hns_patch = cell.select_atoms('type 2')  # the hns patches
    # defining collectors
    # bug:
    gyr_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # bond info
    n_bonds = len(bug.bonds.indices)
    bond_lengths = np.zeros((n_bonds, 1), dtype=np.float64)
    cosine_corrs = np.zeros(n_bonds, dtype=np.float64)
    # H-NS binding:
    cis_threshold = 4
    dist_m_hpatch = []
    lj_cut = 2**(1/6)
    r_cutoff = np.round(
        0.5 * lj_cut * (sim_info.dmon + sim_info.dhns_patch), 3
        )
    binding_stats_t = {
        'n_m_hpatch_bound': [],
        'n_hpatch_free': [],
        'n_hpatch_engaged': [],
        'n_hcore_free': [],
        'n_hcore_bridge': [],
        'n_hcore_dangle': [],
        'n_hcore_cis': [],
        'n_hcore_trans': []
    }
    if sim_info.topology == 'linear':
        loop_length_hist_t = np.zeros(sim_info.nmon, dtype=int)
    elif sim_info.topology == 'ring':
        loop_length_hist_t = np.zeros((sim_info.nmon//2)+1, dtype=int)
    else:
        raise ValueError(
            f"The genomic distance is not defined for '{topology}' topology"
        )
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain sizes
        gyr_t[idx] = bug.radius_of_gyration()
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=True)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        # bond info
        bond_dummy, cosine_dummy = correlations.bond_info(
            bug,
            sim_info.topology
            )
        bond_lengths += bond_dummy
        cosine_corrs += cosine_dummy
        # bug - hns patch:
        # distance matrices
        dummy = mda_dist.distance_array(bug, hns_patch)
        dist_m_hpatch.append(dummy)
        d_contact_m_hpatch = clusters.find_direct_contacts(
            dummy, r_cutoff, inclusive=False
            )
        binding_stats_t, loop_length_hist_t = clusters.hns_binding(
                d_contact_m_hpatch,
                sim_info.topology,
                cis_threshold=cis_threshold,
                binding_stats=binding_stats_t,
                loop_length_hist=loop_length_hist_t
                )
    # Saving collectors to memory
    # bug
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
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
    # H-NS binding stats:
    binding_stats_names = {
        'n_m_hpatch_bound': 'nBoundTHnsPatch',
        'n_hpatch_free': 'nFreeTHnsPatch',
        'n_hpatch_engaged': 'nEngagedTHnsPatch',
        'n_hcore_free': 'nFreeTHnsCore',
        'n_hcore_bridge': 'nBridgeTHnsCore',
        'n_hcore_dangle': 'nDangleTHnsCore',
        'n_hcore_cis': 'nCisTHnsCore',
        'n_hcore_trans': 'nTransTHnsCore',
    }
    for key, value in binding_stats_t.items():
        np.save(
            save_to + sim_name + '-' + binding_stats_names[key] + '.npy',
            np.array(value)
        )
    np.save(
        save_to + sim_name + '-loopLengthHistMon.npy',
        np.array(loop_length_hist_t)
        )
    # distance matirx
    np.save(
        save_to + sim_name + '-distMatTMonHnsPatch.npy',
        np.array(dist_m_hpatch)
        )
    print('done.')


def hns_nucleoid_cub_dis_hc_hc_cluster(
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
    polymer (called "bug") and several H-NS proteins (each formed of a core and
    two poles or patches) crowded by soft LJ spheres (crowders) in free
    ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to
    the center of mass of the single polymer (bug).

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    cluster_cutoff = sim_info.dhns + sim_info.dcrowd
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups:
    hns_core: mda.AtomGroup = cell.select_atoms('type 3')  # the hns cores
    # defining collectors
    bonds_t = np.zeros((n_frames, sim_info.nhns), dtype=int)
    clusters_t = np.zeros((n_frames, sim_info.nhns + 1), dtype=int)
    for idx, _ in enumerate(sliced_trj):
        # distance matrices
        dist_mat = clusters.self_dist_array(hns_core.positions)
        # keep atom ids on the diag
        dir_contacts = clusters.find_direct_contacts(dist_mat, cluster_cutoff)
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        bonds_t[idx] = bonds_stat
        contacts = clusters.generate_contact_matrix_new(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        clusters_t[idx] = clusters_stat
    # distance matirx
    np.save(save_to + sim_name + '-bondsHistDirDepTHnsCore.npy', bonds_t)
    np.save(save_to + sim_name + '-clustersHistDirDepTHnsCore.npy', clusters_t)
    print('done.')


def hns_cub_all(
    topology: str,
    trajectory: str,
    lineage: str,
    save_to: str = "./",
    continuous: Optional[bool] = False
) -> None:
    """Runs various analyses on a `lineage` simulation of an 'all' atom
    group in the `geometry` of interest,and saves a variety of outputs (mostly
    in the csv format) to the `save_to` directory.

    Note
    ----
    In the HNS-DNA-Crowder project, we have a single semi-flexible ring
    polymer (called "bug") and several H-NS proteins (each formed of a core and
    two poles or patches) crowded by soft LJ spheres (crowders) in free
    ("cubic" geometry) or confined ("cylindrical" geometry) space.

    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, the coordinated are recentered with respect to the
    center of mass of the single polymer (bug).

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
        "rEdge": {  # edges for distance r from the box center
            "bin_size":  0.1 * min(
                sim_info.dmon, sim_info.dhns, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.lcube
            },
        "lEdge": {
            # edges for 2d hist in x and y direction
            "bin_size":  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcube,
            "lmax": 0.5 * sim_info.lcube
            }
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(
        sim_info.mmon * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    hns: mda.AtomGroup = cell.select_atoms('type 3')  # hns cores
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_hns_info = fixedsize_bins(
        sim_name,
        'rEdgeHns',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # hns
    # # xy
    xy_hist_hns_info = {
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
    xy_hist_hns_info['collector'] = np.zeros(xy_hist_hns_info['n_bins'])
    xy_hist_hns_info['collector'] *= 0
    # # xz
    xz_hist_hns_info = {
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
    xz_hist_hns_info['collector'] = np.zeros(xz_hist_hns_info['n_bins'])
    xz_hist_hns_info['collector'] *= 0
    # # yz
    yz_hist_hns_info = {
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
    yz_hist_hns_info['collector'] = np.zeros(yz_hist_hns_info['n_bins'])
    yz_hist_hns_info['collector'] *= 0
    for _ in sliced_trj:
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions, axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
        )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions, axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
        )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # hns
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(hns.positions, axis=1),
            bins=r_hist_hns_info['bin_edges'],
            range=r_hist_hns_info['range']
        )
        r_hist_hns_info['collector'] += pos_hist
        r_hist_hns_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 0],
            hns.positions[:, 1],
            bins=xy_hist_hns_info['bin_edges'],
            range=xy_hist_hns_info['range'],
        )
        xy_hist_hns_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 0],
            hns.positions[:, 2],
            bins=xz_hist_hns_info['bin_edges'],
            range=xz_hist_hns_info['range'],
        )
        xz_hist_hns_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 1],
            hns.positions[:, 2],
            bins=yz_hist_hns_info['bin_edges'],
            range=yz_hist_hns_info['range'],
        )
        yz_hist_hns_info['collector'] += pos_hist

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
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Hns': xy_hist_hns_info,
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Hns': xz_hist_hns_info,
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Hns': yz_hist_hns_info,
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info: HnsCyl = HnsCyl(
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
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
        n_frames = cell.trajectory.n_frames - 1
    else:
        sliced_trj = cell.trajectory
        n_frames = cell.trajectory.n_frames
    # selecting atom groups:
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    hns_patch = cell.select_atoms('type 2')  # hns patches
    # defining collectors
    # bug:
    gyr_t = np.zeros(n_frames)
    fsd_t = np.zeros(n_frames)
    trans_size_t = np.zeros(n_frames)
    principal_axes_t = np.empty([n_frames, 3, 3])
    asphericity_t = np.zeros(n_frames)
    shape_parameter_t = np.zeros(n_frames)
    # bond info
    n_bonds = len(bug.bonds.indices)
    bond_lengths = np.zeros((n_bonds, 1), dtype=np.float64)
    cosine_corrs = np.zeros(n_bonds, dtype=np.float64)
    # H-NS binding:
    cis_threshold = 4
    dist_m_hpatch = []
    lj_cut = 2**(1/6)
    r_cutoff = np.round(
        0.5 * lj_cut * (sim_info.dmon + sim_info.dhns_patch), 3
        )
    binding_stats_t = {
        'n_m_hpatch_bound': [],
        'n_hpatch_free': [],
        'n_hpatch_engaged': [],
        'n_hcore_free': [],
        'n_hcore_bridge': [],
        'n_hcore_dangle': [],
        'n_hcore_cis': [],
        'n_hcore_trans': []
    }
    if sim_info.topology == 'linear':
        loop_length_hist_t = np.zeros(sim_info.nmon, dtype=int)
    elif sim_info.topology == 'ring':
        loop_length_hist_t = np.zeros((sim_info.nmon//2)+1, dtype=int)
    else:
        raise ValueError(
            f"The genomic distance is not defined for '{topology}' topology"
        )
    for idx, _ in enumerate(sliced_trj):
        # bug:
        # various measures of chain size
        gyr_t[idx] = bug.radius_of_gyration()
        trans_size_t[idx] = transverse_size(bug)
        fsd_t[idx] = fsd(bug.positions)
        # shape parameters:
        asphericity_t[idx] = bug.asphericity(wrap=False, unwrap=False)
        shape_parameter_t[idx] = bug.shape_parameter(wrap=False)
        principal_axes_t[idx] = bug.principal_axes(wrap=False)
        # bond info
        bond_dummy, cosine_dummy = correlations.bond_info(
            bug,
            sim_info.topology
            )
        bond_lengths += bond_dummy
        cosine_corrs += cosine_dummy
        # bug - hns patch:
        # distance matrices
        dummy = mda_dist.distance_array(bug, hns_patch)
        dist_m_hpatch.append(dummy)
        d_contact_m_hpatch = clusters.find_direct_contacts(
            dummy, r_cutoff, inclusive=False
            )
        binding_stats_t, loop_length_hist_t = clusters.hns_binding(
                d_contact_m_hpatch,
                sim_info.topology,
                cis_threshold=cis_threshold,
                binding_stats=binding_stats_t,
                loop_length_hist=loop_length_hist_t
                )

    # Saving collectors to memory
    # bug
    np.save(save_to + sim_name + '-transSizeTMon.npy', trans_size_t)
    np.save(save_to + sim_name + '-fsdTMon.npy', fsd_t)
    np.save(save_to + sim_name + '-gyrTMon.npy', gyr_t)
    np.save(save_to + sim_name + '-asphericityTMon.npy', asphericity_t)
    np.save(save_to + sim_name + '-shapeTMon.npy', shape_parameter_t)
    np.save(save_to + sim_name + '-principalTMon.npy', principal_axes_t)
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
    # H-NS binding stats:
    binding_stats_names = {
        'n_m_hpatch_bound': 'nBoundTHnsPatch',
        'n_hpatch_free': 'nFreeTHnsPatch',
        'n_hpatch_engaged': 'nEngagedTHnsPatch',
        'n_hcore_free': 'nFreeTHnsCore',
        'n_hcore_bridge': 'nBridgeTHnsCore',
        'n_hcore_dangle': 'nDangleTHnsCore',
        'n_hcore_cis': 'nCisTHnsCore',
        'n_hcore_trans': 'nTransTHnsCore',
    }
    for key, value in binding_stats_t.items():
        np.save(
            save_to + sim_name + '-' + binding_stats_names[key] + '.npy',
            np.array(value)
        )
    np.save(
        save_to + sim_name + '-loopLengthHistMon.npy',
        np.array(loop_length_hist_t)
        )
    # distance matirx
    np.save(
        save_to + sim_name + '-distMatTMonHnsPatch.npy',
        np.array(dist_m_hpatch)
        )
    print('done.')


def hns_cyl_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = HnsCyl(
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
        "rEdge": {
            "bin_size":  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.dcyl
            },
        "zEdge": {
            "bin_size":  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcyl,
            "lmax": 0.5 * sim_info.lcyl
            },
        "thetaEdge": {
            "bin_size":  np.pi / 36,
            "lmin": -1 * np.pi,
            "lmax": np.pi
            },
        "lEdge": {
            # edges for 2d hist in x and y direction
            "bin_size":  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.dcyl,
            "lmax": 0.5 * sim_info.dcyl
            }
    }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dmon * np.sqrt(
        sim_info.mmon * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt,
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('resid 0')  # crowders
    bug: mda.AtomGroup = cell.select_atoms('resid 1')  # chain/monomers
    hns: mda.AtomGroup = cell.select_atoms('type 3')  # hns cores
    # bin edges and histograms in different directions:
    # radial direction of the cylindrical coordinate system
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_hns_info = fixedsize_bins(
        sim_name,
        'rEdgeHns',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # longitudinal direction of the cylindrical coordinate system
    z_hist_crd_info = fixedsize_bins(
        sim_name,
        'zEdgeCrd',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    z_hist_mon_info = fixedsize_bins(
        sim_name,
        'zEdgeMon',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    z_hist_hns_info = fixedsize_bins(
        sim_name,
        'zEdgeHns',
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # theta of the cylindrical coordinate system
    theta_hist_crd_info = fixedsize_bins(
        sim_name,
        'thetaEdgeCrd',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    theta_hist_mon_info = fixedsize_bins(
        sim_name,
        'thetaEdgeMon',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    theta_hist_hns_info = fixedsize_bins(
        sim_name,
        'thetaEdgeHns',
        bin_edges["thetaEdge"]["bin_size"],
        bin_edges["thetaEdge"]["lmin"],
        bin_edges["thetaEdge"]["lmax"],
        bin_type="periodic",
        save_to=save_to
        )  # in radians
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["zEdge"]["bin_size"],
        bin_edges["zEdge"]["lmin"],
        bin_edges["zEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd_info['collector'].any() != 0,
            r_hist_hns_info['collector'].any() != 0,
            r_hist_mon_info['collector'].any() != 0,
            z_hist_crd_info['collector'].any() != 0,
            z_hist_hns_info['collector'].any() != 0,
            z_hist_mon_info['collector'].any() != 0,
            theta_hist_crd_info['collector'].any() != 0,
            theta_hist_hns_info['collector'].any() != 0,
            theta_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            r_hist_hns_info['collector_std'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
            z_hist_crd_info['collector_std'].any() != 0,
            z_hist_hns_info['collector_std'].any() != 0,
            z_hist_mon_info['collector_std'].any() != 0,
            theta_hist_crd_info['collector_std'].any() != 0,
            theta_hist_hns_info['collector_std'].any() != 0,
            theta_hist_mon_info['collector_std'].any() != 0,
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    # hns
    # # xy
    xy_hist_hns_info = {
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
    xy_hist_hns_info['collector'] = np.zeros(xy_hist_hns_info['n_bins'])
    xy_hist_hns_info['collector'] *= 0
    # # xz
    xz_hist_hns_info = {
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
    xz_hist_hns_info['collector'] = np.zeros(xz_hist_hns_info['n_bins'])
    xz_hist_hns_info['collector'] *= 0
    # # yz
    yz_hist_hns_info = {
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
    yz_hist_hns_info['collector'] = np.zeros(yz_hist_hns_info['n_bins'])
    yz_hist_hns_info['collector'] *= 0
    for _ in sliced_trj:
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions[:, :2], axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
        )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            crds.positions[:, 2],
            bins=z_hist_crd_info['bin_edges'],
            range=z_hist_crd_info['range']
        )
        z_hist_crd_info['collector'] += pos_hist
        z_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # theta in radians
        pos_hist, _ = np.histogram(
            np.arctan2(crds.positions[:, 1], crds.positions[:, 0]),
            bins=theta_hist_crd_info['bin_edges'],
            range=theta_hist_crd_info['range']
        )
        theta_hist_crd_info['collector'] += pos_hist
        theta_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions[:, :2], axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
        )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            bug.positions[:, 2],
            bins=z_hist_mon_info['bin_edges'],
            range=z_hist_mon_info['range']
        )
        z_hist_mon_info['collector'] += pos_hist
        z_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # theta in radian
        pos_hist, _ = np.histogram(
            np.arctan2(bug.positions[:, 1], bug.positions[:, 0]),
            bins=theta_hist_mon_info['bin_edges'],
            range=theta_hist_mon_info['range']
        )
        theta_hist_mon_info['collector'] += pos_hist
        theta_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist
        # hns
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(hns.positions[:, :2], axis=1),
            bins=r_hist_hns_info['bin_edges'],
            range=r_hist_hns_info['range']
        )
        r_hist_hns_info['collector'] += pos_hist
        r_hist_hns_info['collector_std'] += np.square(pos_hist)
        # # z
        pos_hist, _ = np.histogram(
            hns.positions[:, 2],
            bins=z_hist_hns_info['bin_edges'],
            range=z_hist_hns_info['range']
        )
        z_hist_hns_info['collector'] += pos_hist
        z_hist_hns_info['collector_std'] += np.square(pos_hist)
        # # theta in radian
        pos_hist, _ = np.histogram(
            np.arctan2(hns.positions[:, 1], hns.positions[:, 0]),
            bins=theta_hist_hns_info['bin_edges'],
            range=theta_hist_hns_info['range']
        )
        theta_hist_hns_info['collector'] += pos_hist
        theta_hist_hns_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 0],
            hns.positions[:, 1],
            bins=xy_hist_hns_info['bin_edges'],
            range=xy_hist_hns_info['range'],
        )
        xy_hist_hns_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 0],
            hns.positions[:, 2],
            bins=xz_hist_hns_info['bin_edges'],
            range=xz_hist_hns_info['range'],
        )
        xz_hist_hns_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            hns.positions[:, 1],
            hns.positions[:, 2],
            bins=yz_hist_hns_info['bin_edges'],
            range=yz_hist_hns_info['range'],
        )
        yz_hist_hns_info['collector'] += pos_hist

    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info,
            'Hns': r_hist_hns_info
        },
        'zHist': {
            'Crd': z_hist_crd_info,
            'Mon': z_hist_mon_info,
            'Hns': z_hist_hns_info
        },
        'thetaHist': {
            'Crd': theta_hist_crd_info,
            'Mon': theta_hist_mon_info,
            'Hns': theta_hist_hns_info
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info,
            'Hns': xy_hist_hns_info,
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info,
            'Hns': xz_hist_hns_info,
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info,
            'Hns': yz_hist_hns_info,
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')


def two_mon_dep_cub_all(
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

    In MDAnalysis, selections by `universe.select_atoms` always return an
    AtomGroup with atoms sorted according to their index in the topology.
    This feature is used below to measure the end-to-end distance (Flory
    radius), genomic distance (index differnce along the backbone), and any
    other measurement that needs the sorted indices of atoms,bonds, angles, and
    any other attribute of an atom.

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
    sim_info = TwoMonDep(
        trajectory,
        lineage,
        'cubic',
        'all',
        'atom'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    print("\n" + sim_name + " is analyzing...\n")
    # dict of bin edges:
    bin_edges = {
        "rEdge": {  # edges for distance r from the box center
            "bin_size":  0.1 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": 0,
            "lmax": 0.5 * sim_info.lcube
            },
        "lEdge": {  # edges in cartesian coordinates
            "bin_size":  0.5 * min(sim_info.dmon, sim_info.dcrowd),
            "lmin": -0.5 * sim_info.lcube,
            "lmax": 0.5 * sim_info.lcube
            },
        }
    # LJ time difference between two consecutive frames:
    time_unit = sim_info.dcrowd * np.sqrt(
        sim_info.mcrowd * sim_info.eps_others)  # LJ time unit
    # Sampling via LAMMPS dump every 'adump', so trajectory dt is:
    sim_real_dt = sim_info.adump * sim_info.dt * time_unit
    cell = mda.Universe(
        topology, trajectory, topology_format='DATA',
        format='LAMMPSDUMP', lammps_coordinate_convention='unscaled',
        atom_style="id resid type x y z", dt=sim_real_dt
        )
    # slicing trajectory based the continuous condition
    if continuous:
        sliced_trj = cell.trajectory[0: -1]
    else:
        sliced_trj = cell.trajectory
    # selecting atom groups
    crds: mda.AtomGroup = cell.select_atoms('type 2')  # crowders
    print(len(crds.atoms))
    bug: mda.AtomGroup = cell.select_atoms('type 1')  # chain/monomers
    print(len(bug.atoms))
    # bin edges and histograms in different directions:
    # distance from the box center
    r_hist_crd_info = fixedsize_bins(
        sim_name,
        'rEdgeCrd',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    r_hist_mon_info = fixedsize_bins(
        sim_name,
        'rEdgeMon',
        bin_edges["rEdge"]["bin_size"],
        bin_edges["rEdge"]["lmin"],
        bin_edges["rEdge"]["lmax"],
        bin_type="nonnegative",
        save_to=save_to
    )
    # x direction of the cartesian coordinate system
    x_hist_info = fixedsize_bins(
        sim_name,
        "xEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # y direction of the cartesian coordinate system
    y_hist_info = fixedsize_bins(
        sim_name,
        "yEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # z direction of the cartesian coordinate system
    z_hist_info = fixedsize_bins(
        sim_name,
        "zEdge",
        bin_edges["lEdge"]["bin_size"],
        bin_edges["lEdge"]["lmin"],
        bin_edges["lEdge"]["lmax"],
        bin_type="ordinary",
        save_to=save_to
    )
    # check if any of the histograms are empty or not.
    if any([
            r_hist_crd_info['collector'].any() != 0,
            r_hist_mon_info['collector'].any() != 0,
            r_hist_crd_info['collector_std'].any() != 0,
            r_hist_mon_info['collector_std'].any() != 0,
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
    xy_hist_crd_info = {
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
    xy_hist_crd_info['collector'] = np.zeros(xy_hist_crd_info['n_bins'])
    xy_hist_crd_info['collector'] *= 0
    # # xz
    xz_hist_crd_info = {
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
    xz_hist_crd_info['collector'] = np.zeros(xz_hist_crd_info['n_bins'])
    xz_hist_crd_info['collector'] *= 0
    # # yz
    yz_hist_crd_info = {
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
    yz_hist_crd_info['collector'] = np.zeros(yz_hist_crd_info['n_bins'])
    yz_hist_crd_info['collector'] *= 0
    # mon
    # # xy
    xy_hist_mon_info = {
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
    xy_hist_mon_info['collector'] = np.zeros(xy_hist_mon_info['n_bins'])
    xy_hist_mon_info['collector'] *= 0
    # # xz
    xz_hist_mon_info = {
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
    xz_hist_mon_info['collector'] = np.zeros(xz_hist_mon_info['n_bins'])
    xz_hist_mon_info['collector'] *= 0
    # # yz
    yz_hist_mon_info = {
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
    yz_hist_mon_info['collector'] = np.zeros(yz_hist_mon_info['n_bins'])
    yz_hist_mon_info['collector'] *= 0
    for _ in sliced_trj:
        # histogram in r direction
        # crds
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(crds.positions, axis=1),
            bins=r_hist_crd_info['bin_edges'],
            range=r_hist_crd_info['range']
            )
        r_hist_crd_info['collector'] += pos_hist
        r_hist_crd_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 1],
            bins=xy_hist_crd_info['bin_edges'],
            range=xy_hist_crd_info['range'],
        )
        xy_hist_crd_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 0],
            crds.positions[:, 2],
            bins=xz_hist_crd_info['bin_edges'],
            range=xz_hist_crd_info['range'],
        )
        xz_hist_crd_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            crds.positions[:, 1],
            crds.positions[:, 2],
            bins=yz_hist_crd_info['bin_edges'],
            range=yz_hist_crd_info['range'],
        )
        yz_hist_crd_info['collector'] += pos_hist
        # bug
        # # r
        pos_hist, _ = np.histogram(
            np.linalg.norm(bug.positions, axis=1),
            bins=r_hist_mon_info['bin_edges'],
            range=r_hist_mon_info['range']
            )
        r_hist_mon_info['collector'] += pos_hist
        r_hist_mon_info['collector_std'] += np.square(pos_hist)
        # # xy
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 1],
            bins=xy_hist_mon_info['bin_edges'],
            range=xy_hist_mon_info['range'],
        )
        xy_hist_mon_info['collector'] += pos_hist
        # # xz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 0],
            bug.positions[:, 2],
            bins=xz_hist_mon_info['bin_edges'],
            range=xz_hist_mon_info['range'],
        )
        xz_hist_mon_info['collector'] += pos_hist
        # # yz
        pos_hist, _, _ = np.histogram2d(
            bug.positions[:, 1],
            bug.positions[:, 2],
            bins=yz_hist_mon_info['bin_edges'],
            range=yz_hist_mon_info['range'],
        )
        yz_hist_mon_info['collector'] += pos_hist

    # end of loop
    hist_1d_groups = {
        'rHist': {
            'Crd': r_hist_crd_info,
            'Mon': r_hist_mon_info
        }
    }
    write_hists(hist_1d_groups, sim_name, save_to, std=True)
    hist_2d_groups = {
        'xyHist': {
            'Crd': xy_hist_crd_info,
            'Mon': xy_hist_mon_info
        },
        'xzHist': {
            'Crd': xz_hist_crd_info,
            'Mon': xz_hist_mon_info
        },
        'yzHist': {
            'Crd': yz_hist_crd_info,
            'Mon': yz_hist_mon_info
        }
    }
    write_hists(hist_2d_groups, sim_name, save_to, std=False)
    print('done.')
