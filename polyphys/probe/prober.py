import warnings
import pathlib
from abc import ABC, abstractmethod
from typing import Dict, Any, Tuple, Literal, Optional, List
import numpy as np
import MDAnalysis as mda
from polyphys.manage.parser import SumRuleCyl, TwoMonDep
from polyphys.manage.typer import ParserInstance
from polyphys.analyze.measurer import (
    transverse_size,
    fsd,
    end_to_end,
    fixedsize_bins,
    axial_histogram,
    radial_cyl_histogram,
    azimuth_cyl_histogram,
    planar_cartesian_histogram
)


class ProberBase(ABC):
    """
    Base class for probing molecular dynamics simulations, providing a
    framework for extracting physical properties and generating reports from
    trajectory and topology files.

    Subclasses must implement methods for defining atom groups, setting up
    collectors, and probing individual frames to extract specific properties.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., files containing physical properties)
        will be saved. Must exist prior to instantiation.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Attributes
    ----------
    topology : str
        Path to the topology file provided during initialization.
    trajectory : str
        Path to the trajectory file provided during initialization.
    lineage : str
        The type of lineage associated with the trajectory.
    save_to : str
        Directory where probing results will be saved.
    continuous : bool
        Indicates whether the trajectory is part of a continuous sequence.
    sim_name : str
        A unique name derived from the parser and group attributes for the
        current simulation.
    sliced_trj : mda.coordinates.base.FrameIteratorSliced
        The sliced trajectory object based on the `continuous` attribute.
    n_frames : int
        The total number of time frames in the trajectory.
    atom_groups : dict of str -> MDAnalysis.AtomGroup
        Dictionary mapping names to MDAnalysis AtomGroup instances for the
        defined atom groups.
    collectors : dict of str -> numpy.ndarray
        Dictionary storing collected physical properties extracted during
        probing.

    Abstract Attributes
    -------------------
    _universe : MDAnalysis.Universe
        MDAnalysis Universe object describing the molecular dynamics system.
    _parser : ParserInstance
        Parser object containing metadata about the simulation.
    _damping_time : float
        Time difference between two consecutive frames in the trajectory.

    Abstract Methods
    ----------------
    _define_atom_groups() -> None
        Defines atom groups specific to the molecular dynamics project.
        Subclasses must set up the `self._atom_groups` dictionary with keys
        representing group names (str) and values as MDAnalysis AtomGroup
        objects.
    _prepare() -> None
        Sets up collectors for storing physical properties.
        Subclasses must implement this method to initialize `self._collectors`.
    _probe_frame(idx: int) -> None
        Probes a single frame for specific atom groups and updates collectors.
        Subclasses must implement this method.

    Methods
    -------
    _set_trj_slice() -> Tuple[int, mda.coordinates.base.FrameIteratorSliced]
        Slices the trajectory based on the `continuous` attribute and returns
        the number of frames and the sliced trajectory object.
    simulation_report() -> None
        Generates and saves a CSV report containing key attributes and frame
        count to `save_to` directory.
    run() -> None
        Iterates over the trajectory frames, applying the `_probe_frame`
        method for data collection.
    save_artifacts() -> None
        Saves all collected analysis data to the `save_to` directory as `.npy`
        files.

    Raises
    ------
    ValueError
        If the directory specified in `save_to` does not exist.
    RuntimeWarning
        If `lineage` is "segment" and `continuous` is False, indicating that
        the trajectory might be part of a sequence.
    AttributeError
        If accessing uninitialized properties `_parser`, `_universe`, or
        `_damping_time`.
    """
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: Literal["segment", "whole"],
        save_to: str,
        continuous: bool = False,
    ) -> None:
        if lineage == "segment" and not continuous:
            warnings.warn(
                f"Lineage is '{self.lineage}' and 'continuous' is False. "
                f"Ensure the trajectory '{trajectory}' is NOT part of a "
                "sequence of trajectories.",
                RuntimeWarning,
            )
        self._topology = topology
        self._trajectory = trajectory
        self._lineage = lineage
        if not pathlib.Path(save_to).is_dir():
            raise ValueError(f"The directory '{save_to}' does not exist.")
        self._save_to = save_to
        self._continuous = continuous
        self._damping_time: Optional[float] = None
        self._parser: Optional[ParserInstance] = None
        self._universe: Optional[mda.Universe] = None
        self._sim_name: str = f"{self.parser.name}-{self.parser.group}"
        self._n_frames, self._sliced_trj = self._set_trj_slice()
        self._atom_groups: Dict[str, mda.AtomGroup] = {}
        self._define_atom_groups()
        self._collectors: Dict[str, Any] = {}
        self._prepare()

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
    def lineage(self) -> str:
        """
        Returns the current lineage.
        """
        return self._lineage

    @property
    def continuous(self) -> bool:
        """
        Indicates whether the trajectory is part of a continuous sequence.
        """
        return self._continuous

    @property
    def save_to(self) -> str:
        """
        Returns the directory where probe results will be saved.
        """
        return self._save_to

    @property
    def sim_name(self) -> str:
        """
        Returns the simulation name derived from the parser and group.
        """
        return self._sim_name

    @property
    def sliced_trj(self) -> mda.coordinates.base.FrameIteratorSliced:
        """
        Returns the sliced trajectory over which the analysis is performed.
        """
        return self._sliced_trj

    @property
    def n_frames(self) -> int:
        """
        Returns the total number of time frames in the trajectory.
        """
        return self._n_frames

    @property
    def atom_groups(self) -> Dict[str, mda.AtomGroup]:
        """
        Returns a dictionary of atom groups and their associated MDAnalysis
        AtomGroup instances.
        """
        return self._atom_groups

    @property
    def universe(self) -> mda.Universe:
        """
        Returns the MDAnalysis Universe object describing the molecular
        dynamics system.
        """
        if self._universe is None:
            raise AttributeError("'_universe' has not been initialized.")
        return self._universe

    @property
    def parser(self) -> ParserInstance:
        """
        Returns the parser object containing metadata about the simulation.
        """
        if self._parser is None:
            raise AttributeError("'_parser' has not been initialized.")
        return self._parser

    @property
    def damping_time(self) -> Optional[float]:
        """
        Returns the time difference between two consecutive trajectory frames.
        """
        if self._damping_time is None:
            raise AttributeError("'_damping_time' has not been initialized.")
        return self._damping_time

    @property
    def collectors(self) -> Dict[str, np.ndarray]:
        """
        Returns a dictionary of collected physical properties extracted from
        the trajectory.
        """
        return self._collectors

    def _set_trj_slice(
        self
    ) -> Tuple[int, mda.coordinates.base.FrameIteratorSliced]:
        """
        Slices the trajectory based on the `continuous` attribute.

        Returns
        -------
        Tuple[int, mda.coordinates.base.FrameIteratorSliced]
            The number of frames and the sliced trajectory.
        """
        if self.continuous:  # Using the public property
            return self.universe.trajectory.n_frames - 1, \
               self.universe.trajectory[0: -1]  # Using the public property
        return self.universe.trajectory.n_frames, self.universe.trajectory

    def simulation_report(self) -> None:
        """
        Generates and saves a simulation report containing key attributes and
        frame count as a CSV file.
        """
        report_name = f"{self.save_to}/{self.sim_name}-stamps.csv"
        with open(report_name, mode="w", encoding="utf-8") as report:
            # Write header
            report.write(",".join(self.parser.attributes) + ",n_frames\n")
            # Write values
            values = [
                getattr(self._parser, attr) for attr in self.parser.attributes
                ]
            report.write(",".join(map(str, values)) + f",{self.n_frames}\n")
        print(f"Simulation report saved to '{report_name}'.")

    @abstractmethod
    def _define_atom_groups(self) -> None:
        """
        Defines atom groups specific to the molecular dynamics project.

        Subclasses must set up the `self._atom_groups` dictionary with keys
        representing group names (str) and values as MDAnalysis AtomGroup
        objects.
        """

    @abstractmethod
    def _prepare(self) -> None:
        """
        Sets up collectors to store physical properties.
        Subclasses must implement this method.
        """

    @abstractmethod
    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single frame for specific atom groups and updates collectors.

        Parameters
        ----------
        idx : int
            Index of the trajectory frame to probe.
        """

    def run(self) -> None:
        """
        Probes all selected trajectory frames and collects data.
        """
        print(f"Probing '{self.parser.name}'...")
        for idx, _ in enumerate(self.sliced_trj):
            self._probe_frame(idx)
        print("Probing complete.")

    def save_artifacts(self) -> None:
        """
        Saves all collected analysis data to the specified output directory.
        """
        for prop, artifact in self.collectors.items():
            filename = f"{self.save_to}{self.sim_name}-{prop}.npy"
            np.save(filename, artifact)
        print("All artifacts saved.")


class TwoMonDepCubBugProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *TwoMonDepCub*
    molecular dynamics project to extract specific physical properties.

    Physical Properties Extracted
    -----------------------------
    - `gyrTMon` : numpy.ndarray
        Radius of gyration of the polymer for each frame.
    - `dxTMon` : numpy.ndarray
        Framewise standard deviation of polymer's position along the x-axis.
    - `dyTMon` : numpy.ndarray
        Framewise standard deviation of polymer's position along the y-axis.
    - `dzTMon` : numpy.ndarray
        Framewise standard deviation of polymer's position along the z-axis.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., physical property files) will be
        saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Notes
    -----
    - Atoms of type `1` in the LAMMPS dump format represent the `bug` atom
      group in this project.
    - Coordinates are wrapped and unscaled.

    Examples
    --------
    Creating an instance of `TwoMonDepCubBugProber` for a specific simulation:

    >>> prober = TwoMonDepCubBugProber(
    ...     topology="topology.data",
    ...     trajectory="trajectory.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=True
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon']

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: Literal["segment", "whole"],
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology,
                         trajectory,
                         lineage,
                         save_to,
                         continuous=continuous)

        self._parser = TwoMonDep(trajectory, lineage, 'bug')
        self._damping_time = \
            getattr(self._parser, 'bdump') * getattr(self._parser, 'dt')
        self._universe = mda.Universe(
            self.topology,
            self.trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id type x y z",
            dt=self.damping_time
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['bug'] = self.universe.select_atoms('type 1')

    def _prepare(self) -> None:
        self._collectors = {
            'gyrTMon': np.zeros(self._n_frames),
            'dxTMon': np.zeros(self._n_frames),
            'dyTMon': np.zeros(self._n_frames),
            'dzTMon': np.zeros(self._n_frames)
        }

    def _probe_frame(self, idx) -> None:
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['bug'].radius_of_gyration()
        self._collectors['dxTMon'][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=0)
        self._collectors['dyTMon'][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=1)
        self._collectors['dzTMon'][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=2)


class SumRuleCylBugProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *SumRuleCyl*
    molecular dynamics project to extract specific physical properties.

    Physical Properties Extracted
    -----------------------------
    - `transSizeTMon` : numpy.ndarray
        Transverse size of the polymer for each frame in the trajectory.
    - `fsdTMon` : numpy.ndarray
        Framewise standard deviation of polymer's position.
    - `gyrTMon` : numpy.ndarray
        Radius of gyration of the polymer for each frame.
    - `rfloryTMon` : numpy.ndarray
        End-to-end distance of the polymer for each frame.
    - `asphericityTMon` : numpy.ndarray
        Asphericity of the polymer's configuration for each frame.
    - `shapeTMon` : numpy.ndarray
        Shape parameter of the polymer for each frame.
    - `principalTMon` : numpy.ndarray
        Principal axes of the polymer for each frame.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., physical property files) will be
        saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Notes
    -----
    - Atoms with `resid 1` in the LAMMPS dump format represent the `bug` atom
      group in this project.
    - Coordinates are wrapped and unscaled. The polymer's center of mass is
      recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `SumRuleCylBugProber` for a specific simulation:

    >>> prober = SumRuleCylBugProber(
    ...     topology="topology.data",
    ...     trajectory="trajectory.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Crd']

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: Literal["segment", "whole"],
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        self._parser = SumRuleCyl(trajectory, lineage, 'bug')
        self._damping_time = (
            getattr(self._parser, 'bdump') * getattr(self._parser, 'dt')
        )
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
        self._atom_groups['bug'] = self.universe.select_atoms('resid 1')

    def _prepare(self) -> None:
        self._collectors = {
            'transSizeTMon': np.zeros(self._n_frames),
            'fsdTMon': np.zeros(self._n_frames),
            'gyrTMon': np.zeros(self._n_frames),
            'rfloryTMon': np.zeros(self._n_frames),
            'asphericityTMon': np.zeros(self._n_frames),
            'shapeTMon': np.zeros(self._n_frames),
            'principalTMon': np.zeros([self._n_frames, 3, 3])
        }

    def _probe_frame(self, idx) -> None:
        self._collectors['transSizeTMon'][idx] = \
            transverse_size(self._atom_groups['bug'].positions, axis=2)
        self._collectors['fsdTMon'][idx] = \
            fsd(self._atom_groups['bug'].positions, axis=2)
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['bug'].radius_of_gyration()
        self._collectors['rfloryTMon'][idx] = \
            end_to_end(self._atom_groups['bug'].positions)
        self._collectors['asphericityTMon'][idx] = \
            self._atom_groups['bug'].asphericity(wrap=False, unwrap=False)
        self._collectors['shapeTMon'][idx] = \
            self._atom_groups['bug'].shape_parameter(wrap=False)
        self._collectors['principalTMon'][idx] = \
            self._atom_groups['bug'].principal_axes(wrap=False)


class SumRuleCylAllProber(ProberBase):
    """
    Runs various analyses on a `lineage` simulation of the LAMMPS 'all' atom
    group in the 'SumRuleCyl' molecular dynamics project.

    Note
    ----
    In this project, coordinates are wrapped and unscaled in a trajectory or
    topology file; moreover, LAMMPS recenter is used to restrict the center of
    mass of 'bug' (monomers) to the center of simulation box; and consequently,
    coordinates of all the particles in a trajectory or topology file is
    recentered to fulfill this constraint.
    """
    _entities = ['Mon', 'Crd']
    _hist1d_directions = ['r', 'z', 'theta']
    _hist2d_planes = ['xy', 'yz', 'zx']

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: Literal["segment", "whole"],
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        self._parser = SumRuleCyl(trajectory, lineage, 'all')
        self._damping_time = (
            getattr(self._parser, 'bdump') * getattr(self._parser, 'dt')
        )
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
        self._atom_groups['Crd'] = self.universe.select_atoms("resid 0")
        self._atom_groups['Mon'] = self.universe.select_atoms("resid 1")

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, List[Any]]] = {}
        bin_edges = {
            'rEdge': {
                'bin_size':  0.1 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                'lmin': 0,
                'lmax': 0.5 * getattr(self._parser, 'lcyl'),
                },
            'zEdge': {
                'bin_size':  0.5 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                'lmin': -0.5 * getattr(self._parser, 'lcyl'),
                'lmax': 0.5 * getattr(self._parser, 'lcyl'),
                },
            'thetaEdge': {
                'bin_size':  np.pi / 36,
                'lmin': -1 * np.pi,
                'lmax': np.pi
                },
            'xEdge': {
                # edges for 2d hist in x and y direction
                'bin_size':  0.1 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                'lmin': -0.5 * getattr(self._parser, 'dcyl'),
                'lmax': 0.5 * getattr(self._parser, 'dcyl')
                },
            'yEdge': {
                # edges for 2d hist in x and y direction
                'bin_size':  0.1 * min(getattr(self._parser, 'dmon'),
                                       getattr(self._parser, 'dcrowd')),
                'lmin': -0.5 * getattr(self._parser, 'dcyl'),
                'lmax': 0.5 * getattr(self._parser, 'dcyl')
                }
        }
        hist1d_bin_type = \
            {'r': 'nonnegative', 'z': 'ordinary', 'theta': 'periodic'}
        hist2d_planes_dirs = {'xy': 'z', 'zx': 'y', 'yz': 'x'}
        for ent in self._entities:
            for dir_, edge_t in hist1d_bin_type.items():
                output = f"{self._save_to}{self._sim_name}-{dir_}Edge{ent}"
                hist_inits = fixedsize_bins(
                    bin_edges[f'{dir_}Edge']['bin_size'],
                    bin_edges[f'{dir_}Edge']['lmin'],
                    bin_edges[f'{dir_}Edge']['lmax'],
                    bin_type=edge_t,
                    save_bin_edges=output
                )
                self._hist1d_props[f'{dir_}Hist{ent}'] = {
                    'n_bins': hist_inits['n_bins'],
                    'bin_edges': hist_inits['bin_edges'],
                    'range': hist_inits['range']
                }
                self._collectors[f'{dir_}Hist{ent}'] = hist_inits['collector']
                self._collectors[f'{dir_}HistStd{ent}'] = \
                    hist_inits['collector_std']
            for plane in hist2d_planes_dirs:
                # the z binning scheme is the same as the one used in 1D hist
                # along z direction above. Thus, we unfortunatley,
                # rewrite the 'zEdge' data to storage.
                self._hist2d_props[f'{plane}Hist{ent}'] = {
                    'n_bins': [],
                    'bin_edges': [],
                    'range': []
                }
                for t_dir in plane:
                    output = \
                        f"{self._save_to}{self._sim_name}-{t_dir}Edge{ent}"
                    hist_inits = fixedsize_bins(
                        bin_edges[f'{t_dir}Edge']['bin_size'],
                        bin_edges[f'{t_dir}Edge']['lmin'],
                        bin_edges[f'{t_dir}Edge']['lmax'],
                        bin_type='ordinary',
                        save_bin_edges=output
                    )
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins'].append(
                        hist_inits['n_bins'])
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins'].append(
                        hist_inits['bin_edges'])
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins'].append(
                        hist_inits['range'])
                self._collectors[f'{plane}Hist{ent}'] = \
                    np.zeros(self._hist2d_props[f'{plane}Hist{ent}']['n_bins'])

    def _probe_frame(self, idx) -> None:
        # r
        for ent in self._entities:
            pos_hist, _ = radial_cyl_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._hist1d_props[f'rHist{ent}']['bin_edges'],
                self._hist1d_props[f'rHist{ent}']['range'],
                2
            )
            self._collectors[f'rHist{ent}'] += pos_hist
            self._collectors[f'rHistStd{ent}'] += np.square(pos_hist)

        # z
        for ent in self._entities:
            pos_hist, _ = axial_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._hist1d_props[f'zHist{ent}']['bin_edges'],
                self._hist1d_props[f'zHist{ent}']['range'],
                2
            )
            self._collectors[f'zHist{ent}']['collector'] += pos_hist
            self._collectors[f'zHistStd{ent}'] += np.square(pos_hist)
        # theta
        for ent in self._entities:
            pos_hist, _ = azimuth_cyl_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._hist1d_props[f'thetaHist{ent}']['bin_edges'],
                self._hist1d_props[f'thetaHist{ent}']['range'],
                2
            )
            self._collectors[f'thetaHist{ent}'] += pos_hist
            self._collectors[f'thetaHistStd{ent}'] += np.square(pos_hist)
        # xy:
        for ent in self._entities:
            pos_hist, _ = planar_cartesian_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._hist2d_props[f'xyHist{ent}']['bin_edges'],
                self._hist2d_props[f'xyHist{ent}']['range'],
                2
            )
            self._collectors[f'xyHist{ent}'] += pos_hist
        # yz
        for ent in self._entities:
            pos_hist, _ = planar_cartesian_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._hist2d_props[f'yzHist{ent}']['bin_edges'],
                self._hist2d_props[f'yzHist{ent}']['range'],
                0
            )
            self._collectors[f'yzHist{ent}'] += pos_hist
        # zx
        for ent in self._entities:
            pos_hist, _ = planar_cartesian_histogram(
                self._atom_groups[f'{ent}'].positions,
                self._collectors[f'zxHist{ent}']['bin_edges'],
                self._collectors[f'zxHist{ent}']['range'],
                1
            )
            self._collectors[f'zxHist{ent}'] += pos_hist
