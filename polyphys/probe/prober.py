import warnings
import pathlib
from abc import ABC, abstractmethod
from typing import Dict, Any, Tuple, Literal, Optional, List
import numpy as np
import MDAnalysis as mda
from polyphys.analyze import clusters, correlations
from polyphys.manage.parser import (
    SumRuleCyl, 
    TransFociCub,
    TransFociCyl,
    SumRuleCubHeteroLinear,
    SumRuleCubHeteroRing,
    HnsCub,
    HnsCyl,
    TwoMonDepCub
)
from polyphys.manage.typer import (
    ParserInstance,
    AxisT,
    BinT,
    DirectionT,
    PlaneT,
    EntityT,
    SumRuleCylInstance,
    TwoMonDepCubInstance,
    TransFociCylInstance,
    TransFociCubInstance,
    SumRuleCubHeteroRingInstance,
    SumRuleCubHeteroLinearInstance,
    HnsCubInstance,
    HnsCylInstance
)
from polyphys.analyze.measurer import (
    transverse_size,
    fsd,
    end_to_end,
    fixedsize_bins,
    radial_histogram,
    radial_cyl_histogram,
    axial_histogram,
    azimuth_cyl_histogram,
    planar_cartesian_histogram
)


class CylindricalHistogramMixIn(ABC):
    """
    Mixin for histograms in cylindrical coordinate systems. Handles histograms
    for `r`, `z`, and `theta`.
    """
    _hist1d_bin_types: Dict[str, str] = {
        'r': 'nonnegative',
        'z': 'ordinary',
        'theta': 'periodic'
    }

    @abstractmethod
    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
        """
        Must be implemented in subclasses to supply cylindrical-specific bin configurations.
        """

    def _initialize_collectors_1d(
        self,
        bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes collectors for 1D histograms in cylindrical coordinates.

        Parameters
        ----------
        bin_edges : Dict[str, Dict[str, Any]]
            Bin edge configurations for cylindrical directions.
        """
        for ent in self._entities:
            for direction, bin_type in self._hist1d_bin_types.items():
                edge_key = f"{direction}Edge"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type=bin_type
                )
                self._hist1d_props[f"{direction}Hist{ent}"] = {
                    'n_bins': hist_data['n_bins'],
                    'bin_edges': hist_data['bin_edges'],
                    'range': hist_data['range'],
                }
                self._collectors[f"{direction}Hist{ent}"] = hist_data['collector']

    def _update_histogram_1d(
        self, direction: str,
        histogram_func,
        entity: str,
        axis: int
    ) -> None:
        """
        Updates a 1D histogram in cylindrical coordinates.

        Parameters
        ----------
        direction : str
            The spatial direction (`r`, `z`, or `theta`).
        histogram_func : Callable
            Function to compute the histogram.
        entity : str
            Particle entity (e.g., `Mon` or `Crd`).
        axis : int
            Spatial axis corresponding to the direction.
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f"{direction}Hist{entity}"]['bin_edges'],
            self._hist1d_props[f"{direction}Hist{entity}"]['range'],
            axis
        )
        self._collectors[f"{direction}Hist{entity}"] += pos_hist


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
        self._save_to = save_to if save_to.endswith('/') else f"{save_to}/"
        self._continuous = continuous
        self._damping_time: Optional[float] = None
        self._parser: Optional[ParserInstance] = None
        self._universe: Optional[mda.Universe] = None
        self._sim_name: str = f"{self.parser.name}-{self.parser.group}"
        self._n_frames, self._sliced_trj = self._set_trj_slice()
        self._atom_groups: Dict[EntityT, mda.AtomGroup] = {}
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
        report_name = f"{self.save_to}{self.sim_name}-stamps.csv"
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
    For the particle entity 'Mon' (representing monomers):
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
    - Atoms with `type 1` in the LAMMPS dump format represent `Mon` (monomers).
    - Coordinates are wrapped and unscaled.

    Examples
    --------
    Creating an instance of `TwoMonDepCubBugProber` for a specific simulation:

    >>> prober = TwoMonDepCubBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=True
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon']

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

        self._parser: TwoMonDepCubInstance = \
            TwoMonDepCub(trajectory, lineage, 'bug')
        self._damping_time = \
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
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
        self._atom_groups['Mon'] = \
            self.universe.select_atoms('type 1')  # monomers

    def _prepare(self) -> None:
        self._collectors = {
            'gyrTMon': np.zeros(self._n_frames),
            'dxTMon': np.zeros(self._n_frames),
            'dyTMon': np.zeros(self._n_frames),
            'dzTMon': np.zeros(self._n_frames)
        }

    def _probe_frame(self, idx) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['Mon'].radius_of_gyration()
        self._collectors['dxTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=0)
        self._collectors['dyTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=1)
        self._collectors['dzTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=2)


class TwoMonDepCubAllProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the *TwoMonDepCub*
    molecular dynamics project to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> ('Mon' for  monomersand 'Crd' for crowders):
    - `rHist<entity>` : numpy.ndarray
        Radial histogram in spherical coordinates.
    - `xyHist<entity>` : numpy.ndarray
        2D histogram in the xy-plane.
    - `yzHist<entity>` : numpy.ndarray
        2D histogram in the yz-plane.
    - `zxHist<entity>` : numpy.ndarray
        2D histogram in the zx-plane.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., histogram files) will be saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Methods
    -------
    _initialize_bin_edges() -> Dict[str, Dict[str, Any]]
        Computes bin edge configurations for histograms.
    _initialize_collectors(bin_edges: Dict[str, Dict[str, Any]]) -> None
        Initializes histogram properties and collectors.
    _update_histogram(direction: DirectionT, histogram_func, entity: EntityT,
                    axis: int) -> None
        Updates a one-dimensional histogram for a specific direction.
    _update_histogram_2d(plane: PlaneT, entity: EntityT) -> None
        Updates a two-dimensional histogram for a specific plane.

    Notes
    -----
    - Atoms with `resid 0` in the LAMMPS dump format represent `Crd`
    (crowders), while atoms with `resid 1` represent `Mon` (monomers).
    - Coordinates are wrapped and unscaled. The center of mass of the `bug`
    group is recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `TwoMonDepCubAllProber` for a specific simulation:

    >>> prober = TwoMonDepCubAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon', 'Crd']
    _hist1d_bin_types: Dict[DirectionT, BinT] = {'r': 'nonnegative'}
    _hist2d_planes_dirs: Dict[PlaneT, DirectionT] = \
        {'xy': 'z', 'yz': 'x', 'zx': 'y'}
    _hist2d_planes_axes: Dict[PlaneT, AxisT] = {'xy': 2, 'yz': 0, 'zx': 1}

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

        self._parser: TwoMonDepCubInstance = \
            TwoMonDepCub(trajectory, lineage, 'all')
        self._damping_time = \
            getattr(self.parser, 'adump') * getattr(self.parser, 'dt')
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id type x y z",
            dt=self.damping_time,
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['Crd'] = \
            self.universe.select_atoms("resid 0")  # crowders
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # small and large monomers

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges = self._initialize_bin_edges()
        self._initialize_collectors(bin_edges)

    def _initialize_bin_edges(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations based on simulation parameters.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """
        dmon = getattr(self._parser, 'dmon')
        dcrowd = getattr(self._parser, 'dcrowd')
        lcube = getattr(self._parser, 'lcube')
        return {
            'rEdge': {'bin_size': 0.1 * min(dmon, dcrowd), 
                      'lmin': 0,
                      'lmax': 0.5 * lcube},
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
            'xEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
            'yEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
        }

    def _initialize_collectors(
            self,
            bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes histogram properties and collectors.

        Parameters
        ----------
        bin_edges : Dict[str, Dict[str, Any]]
            A dictionary of bin edge configurations for spatial dimensions.
        """
        for ent in self._entities:
            # 1D histograms
            for direction, bin_type in self._hist1d_bin_types.items():
                edge_key = f"{direction}Edge"
                output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type=bin_type,
                    save_bin_edges=output,
                )
                self._hist1d_props[f'{direction}Hist{ent}'] = {
                    'n_bins': hist_data['n_bins'],
                    'bin_edges': hist_data['bin_edges'],
                    'range': hist_data['range']
                    }
                self._collectors[f'{direction}Hist{ent}'] = \
                    hist_data['collector']
                self._collectors[f'{direction}HistStd{ent}'] = \
                    hist_data['collector_std']

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._hist2d_props[f'{plane}Hist{ent}'] = \
                    {'n_bins': [], 'bin_edges': [], 'range': []}
                for axis in plane:
                    edge_key = f"{axis}Edge"
                    output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                    hist_data = fixedsize_bins(
                        bin_edges[edge_key]['bin_size'],
                        bin_edges[edge_key]['lmin'],
                        bin_edges[edge_key]['lmax'],
                        bin_type='ordinary',
                        save_bin_edges=output,
                    )
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']\
                        .append(hist_data['n_bins'])
                    self._hist2d_props[f'{plane}Hist{ent}']['bin_edges']\
                        .append(hist_data['bin_edges'])
                    self._hist2d_props[f'{plane}Hist{ent}']['range']\
                        .append(hist_data['range'])

                self._collectors[f'{plane}Hist{ent}'] = np.zeros(
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']
                )
            print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _update_histogram(
        self,
        direction: DirectionT,
        histogram_func,
        entity: EntityT,
    ) -> None:
        """
        Updates a one-dimensional histogram for a specific direction.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction ('r', 'z', or 'theta') for the histogram.
        histogram_func : Callable
            The histogram computation function to use (e.g., radial histogram).
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f'{direction}Hist{entity}']['bin_edges'],
            self._hist1d_props[f'{direction}Hist{entity}']['range']
        )
        self._collectors[f'{direction}Hist{entity}'] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)

    def _update_histogram_2d(
        self,
        plane: PlaneT,
        entity: EntityT,
    ) -> None:
        """
        Updates a two-dimensional histogram for a Cartesian plane.

        Parameters
        ----------
        plane : PlaneT
            A Cartesian plane ('xy', 'yz', or 'zx') for the histogram.
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = planar_cartesian_histogram(
            self._atom_groups[entity].positions,
            self._hist2d_props[f'{plane}Hist{entity}']['bin_edges'],
            self._hist2d_props[f'{plane}Hist{entity}']['range'],
            self._hist2d_planes_axes[plane],
        )
        self._collectors[f'{plane}Hist{entity}'] += pos_hist

    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        for ent in self._entities:
            # 1D histograms
            self._update_histogram('r', radial_histogram, ent)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)



class SumRuleCylBugProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *SumRuleCyl*
    molecular dynamics project to extract specific physical properties.

    Physical Properties Extracted
    -----------------------------
    For the particle entity 'Mon' (representing monomers):
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
    - Atoms with `resid 1` in the LAMMPS dump format represent `Mon`(monomers).
    - Coordinates are wrapped and unscaled. The polymer's center of mass is
      recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `SumRuleCylBugProber` for a specific simulation:

    >>> prober = SumRuleCylBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon']

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
        self._parser: SumRuleCylInstance = \
            SumRuleCyl(trajectory, lineage, 'bug')
        self._damping_time = \
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['Mon'] = \
            self.universe.select_atoms('resid 1')  # monomers

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
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collectors['transSizeTMon'][idx] = \
            transverse_size(self._atom_groups['Mon'].positions, axis=2)
        self._collectors['fsdTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=2)
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['Mon'].radius_of_gyration()
        self._collectors['rfloryTMon'][idx] = \
            end_to_end(self._atom_groups['Mon'].positions)
        self._collectors['asphericityTMon'][idx] = \
            self._atom_groups['Mon'].asphericity(wrap=False, unwrap=False)
        self._collectors['shapeTMon'][idx] = \
            self._atom_groups['Mon'].shape_parameter(wrap=False)
        self._collectors['principalTMon'][idx] = \
            self._atom_groups['Mon'].principal_axes(wrap=False)


class SumRuleCylAllProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the *SumRuleCyl*
    molecular dynamics project to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> ('Mon' for monomers and 'Crd' for crowders):
    - `rHist<entity>` : numpy.ndarray
        Radial histogram in cylindrical coordinates.
    - `zHist<entity>` : numpy.ndarray
        Axial histogram in cylindrical coordinates.
    - `thetaHist<entity>` : numpy.ndarray
        Azimuthal histogram in cylindrical coordinates.
    - `xyHist<entity>` : numpy.ndarray
        2D histogram in the xy-plane.
    - `yzHist<entity>` : numpy.ndarray
        2D histogram in the yz-plane.
    - `zxHist<entity>` : numpy.ndarray
        2D histogram in the zx-plane.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., histogram files) will be saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Methods
    -------
    _initialize_bin_edges() -> Dict[str, Dict[str, Any]]
        Computes bin edge configurations for histograms.
    _initialize_collectors(bin_edges: Dict[str, Dict[str, Any]]) -> None
        Initializes histogram properties and collectors.
    _update_histogram(direction: DirectionT, histogram_func, entity: EntityT,
                    axis: int) -> None
        Updates a one-dimensional histogram for a specific direction.
    _update_histogram_2d(plane: PlaneT, entity: EntityT) -> None
        Updates a two-dimensional histogram for a specific plane.

    Notes
    -----
    - Atoms with `resid 0` in the LAMMPS dump format represent `Crd`
    (crowders), while atoms with `resid 1` represent `Mon` (monomers).
    - Coordinates are wrapped and unscaled. The center of mass of the `bug`
    group is recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `SumRuleCylAllProber` for a specific simulation:

    >>> prober = SumRuleCylAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon', 'Crd']
    _hist1d_bin_types: Dict[DirectionT, BinT] = \
        {'r': 'nonnegative', 'z': 'ordinary', 'theta': 'periodic'}
    _hist2d_planes_dirs: Dict[PlaneT, DirectionT] = \
        {'xy': 'z', 'yz': 'x', 'zx': 'y'}
    _hist2d_planes_axes: Dict[PlaneT, AxisT] = {'xy': 2, 'yz': 0, 'zx': 1}

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

        self._parser: SumRuleCylInstance = \
            SumRuleCyl(trajectory, lineage, 'all')
        self._damping_time = \
            getattr(self.parser, 'adump') * getattr(self.parser, 'dt')
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time,
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['Crd'] = \
            self.universe.select_atoms("resid 0")  # crowders
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # monomers

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges = self._initialize_bin_edges()
        self._initialize_collectors(bin_edges)

    def _initialize_bin_edges(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations based on simulation parameters.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """
        dmon = getattr(self._parser, 'dmon')
        dcrowd = getattr(self._parser, 'dcrowd')
        lcyl = getattr(self._parser, 'lcyl')
        dcyl = getattr(self._parser, 'dcyl')
        return {
            'rEdge': {'bin_size': 0.1 * min(dmon, dcrowd), 'lmin': 0,
                      'lmax': 0.5 * lcyl},
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd), 'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'thetaEdge': {'bin_size': np.pi / 36, 'lmin': -np.pi,
                          'lmax': np.pi},
            'xEdge': {'bin_size': 0.1 * min(dmon, dcrowd), 'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
            'yEdge': {'bin_size': 0.1 * min(dmon, dcrowd), 'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
        }

    def _initialize_collectors(
            self,
            bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes histogram properties and collectors.

        Parameters
        ----------
        bin_edges : Dict[str, Dict[str, Any]]
            A dictionary of bin edge configurations for spatial dimensions.
        """
        for ent in self._entities:
            # 1D histograms
            for direction, bin_type in self._hist1d_bin_types.items():
                edge_key = f"{direction}Edge"
                output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type=bin_type,
                    save_bin_edges=output,
                )
                self._hist1d_props[f'{direction}Hist{ent}'] = {
                    'n_bins': hist_data['n_bins'],
                    'bin_edges': hist_data['bin_edges'],
                    'range': hist_data['range']
                    }
                self._collectors[f'{direction}Hist{ent}'] = \
                    hist_data['collector']
                self._collectors[f'{direction}HistStd{ent}'] = \
                    hist_data['collector_std']

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._hist2d_props[f'{plane}Hist{ent}'] = \
                    {'n_bins': [], 'bin_edges': [], 'range': []}
                for axis in plane:
                    edge_key = f"{axis}Edge"
                    output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                    hist_data = fixedsize_bins(
                        bin_edges[edge_key]['bin_size'],
                        bin_edges[edge_key]['lmin'],
                        bin_edges[edge_key]['lmax'],
                        bin_type='ordinary',
                        save_bin_edges=output,
                    )
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']\
                        .append(hist_data['n_bins'])
                    self._hist2d_props[f'{plane}Hist{ent}']['bin_edges']\
                        .append(hist_data['bin_edges'])
                    self._hist2d_props[f'{plane}Hist{ent}']['range']\
                        .append(hist_data['range'])

                self._collectors[f'{plane}Hist{ent}'] = np.zeros(
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']
                )
            print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _update_histogram(
        self,
        direction: DirectionT,
        histogram_func,
        entity: EntityT,
        axis: int
    ) -> None:
        """
        Updates a one-dimensional histogram for a specific direction.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction ('r', 'z', or 'theta') for the histogram.
        histogram_func : Callable
            The histogram computation function to use (e.g., radial histogram).
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        axis : int
            The spatial axis corresponding to the direction.
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f'{direction}Hist{entity}']['bin_edges'],
            self._hist1d_props[f'{direction}Hist{entity}']['range'],
            axis,
        )
        self._collectors[f'{direction}Hist{entity}'] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)

    def _update_histogram_2d(
        self,
        plane: PlaneT,
        entity: EntityT,
    ) -> None:
        """
        Updates a two-dimensional histogram for a Cartesian plane.

        Parameters
        ----------
        plane : PlaneT
            A Cartesian plane ('xy', 'yz', or 'zx') for the histogram.
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = planar_cartesian_histogram(
            self._atom_groups[entity].positions,
            self._hist2d_props[f'{plane}Hist{entity}']['bin_edges'],
            self._hist2d_props[f'{plane}Hist{entity}']['range'],
            self._hist2d_planes_axes[plane],
        )
        self._collectors[f'{plane}Hist{entity}'] += pos_hist

    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        for ent in self._entities:
            # 1D histograms
            self._update_histogram('r', radial_cyl_histogram, ent, axis=2)
            self._update_histogram('z', axial_histogram, ent, axis=2)
            self._update_histogram('theta', azimuth_cyl_histogram, ent, axis=2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class TransFociCylBugProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *TransFociCyl*
    molecular dynamics project to extract specific physical properties.

    Physical Properties Extracted
    -----------------------------
    The particle entity 'Mon' represents both small and large monomers, while 
    the particle entity 'Foci' represents large monomers. For these entities, 
    the following physical properties are extracted:

    For the `Mon` entity:
    - `transSizeTMon` : numpy.ndarray
        Transverse size of the polymer for each frame in the trajectory.
    - `fsdTMon` : numpy.ndarray
        Framewise standard deviation of the polymer's position.
    - `gyrTMon` : numpy.ndarray
        Radius of gyration of the polymer for each frame.
    - `asphericityTMon` : numpy.ndarray
        Asphericity of the polymer's configuration for each frame.
    - `shapeTMon` : numpy.ndarray
        Shape parameter of the polymer for each frame.
    - `principalTMon` : numpy.ndarray
        Principal axes of the polymer for each frame.

    For the `Foci` entity:
    - `distMatTFoci` : numpy.ndarray
        Pairwise distance matrix for large monomers in the `Foci` entity.
    - `directContactsMatTFoci` : numpy.ndarray
        Direct contact matrix for large monomers based on the cut-off distance.
    - `bondsHistTFoci` : numpy.ndarray
        Histogram of the number of direct contacts for each large monomer.
    - `clustersHistTFoci` : numpy.ndarray
        Histogram of cluster sizes for large monomers in each frame.

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
    - Atoms with `resid 1` in the LAMMPS dump format represent `Mon`(monomers).
    - Coordinates are wrapped and unscaled. The polymer's center of mass is
      recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `SumRuleCylBugProber` for a specific simulation:

    >>> prober = TransFociCylBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Foci', 'Mon']

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
        self._parser: TransFociCylInstance = TransFociCyl(trajectory, lineage, 'bug')
        self._damping_time = (
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
        )
        self._cluster_cutoff = (
            getattr(self.parser, 'dmon_large') + getattr(self.parser, 'dcrowd')
        ) 
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time
        )

    @property
    def cluster_cutoff(self) -> float:
        """
        Return the cut-off distance for considering two large monomers in
        direct contact.
        """
        return self._cluster_cutoff

    def _define_atom_groups(self) -> None:
        self._atom_groups['Mon'] = \
            self.universe.select_atoms('resid 1')  # small and large monomers
        self._atom_groups['Foci'] = \
            self.universe.select_atoms('type 2')  # large monomers

    def _prepare(self) -> None:
        nmon_large = getattr(self.parser, 'dmon_large')
        self._collectors = {
            'transSizeTMon': np.zeros(self._n_frames),
            'fsdTMon': np.zeros(self._n_frames),
            'gyrTMon': np.zeros(self._n_frames),
            'asphericityTMon': np.zeros(self._n_frames),
            'shapeTMon': np.zeros(self._n_frames),
            'principalTMon': np.zeros([self._n_frames, 3, 3]),
            'distMatTFoci': np.zeros((self.n_frames, nmon_large, nmon_large)),
            'directContactsMatTFoci': \
                np.zeros((self.n_frames, nmon_large, nmon_large), dtype=int),
            'bondsHistTFoci': np.zeros((self.n_frames, nmon_large), dtype=int),
            'clustersHistTFoci': \
                np.zeros((self.n_frames, nmon_large+1), dtype=int),
        }

    def _probe_frame(self, idx) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collectors['transSizeTMon'][idx] = \
            transverse_size(self._atom_groups['Mon'].positions, axis=2)
        self._collectors['fsdTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=2)
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['Mon'].radius_of_gyration()
        self._collectors['asphericityTMon'][idx] = \
            self._atom_groups['Mon'].asphericity(wrap=False, unwrap=False)
        self._collectors['shapeTMon'][idx] = \
            self._atom_groups['Mon'].shape_parameter(wrap=False)
        self._collectors['principalTMon'][idx] = \
            self._atom_groups['Mon'].principal_axes(wrap=False)
        dist_mat = \
            clusters.self_dist_array(self._atom_groups['Foci'].positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        np.fill_diagonal(foci_pair_dist, self._atom_groups['Foci'].atoms.ids)
        self._collectors['distMatTFoci'][idx] = foci_pair_dist
        dir_contacts = \
            clusters.find_direct_contacts(dist_mat, self.cluster_cutoff)
        self._collectors['directContactsMatTFoci'][idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        self._collectors['bondsHistTFoci'][idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        self._collectors['clustersHistTFoci'][idx] = clusters_stat

class TransFociCylAllProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the *TransFociCyl*
    molecular dynamics project to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> (`Mon` for both small and large monomers, `Dna` for small 
    monomers, `Foci` for large monomers, and `Crd` for crowders), the following 
    spatial distributions are extracted:

    - `rHist<entity>` : numpy.ndarray
        Radial histogram in cylindrical coordinates.
    - `zHist<entity>` : numpy.ndarray
        Axial histogram in cylindrical coordinates.
    - `thetaHist<entity>` : numpy.ndarray
        Azimuthal histogram in cylindrical coordinates.
    - `xyHist<entity>` : numpy.ndarray
        2D histogram in the xy-plane.
    - `yzHist<entity>` : numpy.ndarray
        2D histogram in the yz-plane.
    - `zxHist<entity>` : numpy.ndarray
        2D histogram in the zx-plane.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., histogram files) will be saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Methods
    -------
    _initialize_bin_edges() -> Dict[str, Dict[str, Any]]
        Computes bin edge configurations for histograms.
    _initialize_collectors(bin_edges: Dict[str, Dict[str, Any]]) -> None
        Initializes histogram properties and collectors.
    _update_histogram(direction: DirectionT, histogram_func, entity: EntityT,
                    axis: int) -> None
        Updates a one-dimensional histogram for a specific direction.
    _update_histogram_2d(plane: PlaneT, entity: EntityT) -> None
        Updates a two-dimensional histogram for a specific plane.

    Notes
    -----
    - Atoms with `resid 0` in the LAMMPS dump format represent `Crd`
    (crowders), while atoms with `resid 1` represent `Mon` (monomers).
    - Coordinates are wrapped and unscaled. The center of mass of the `bug`
    group is recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `TransFociCylAllProber` for a specific simulation:

    >>> prober = TransFociCylAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon', 'Dna', 'Foci', 'Crd']
    _hist1d_bin_types: Dict[DirectionT, BinT] = \
        {'r': 'nonnegative', 'z': 'ordinary', 'theta': 'periodic'}
    _hist2d_planes_dirs: Dict[PlaneT, DirectionT] = \
        {'xy': 'z', 'yz': 'x', 'zx': 'y'}
    _hist2d_planes_axes: Dict[PlaneT, AxisT] = {'xy': 2, 'yz': 0, 'zx': 1}

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

        self._parser: TransFociCylInstance = \
            TransFociCyl(trajectory, lineage, 'all')
        self._damping_time = \
            getattr(self.parser, 'adump') * getattr(self.parser, 'dt')
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time,
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['Crd'] = self.universe.select_atoms("resid 0")
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # small and large monomers
        self._atom_groups['Dna'] = \
            self.universe.select_atoms("type 1")  # small monomers
        self._atom_groups['Foci'] = \
            self.universe.select_atoms("type 2")  # large monomers

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges = self._initialize_bin_edges()
        self._initialize_collectors(bin_edges)

    def _initialize_bin_edges(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations based on simulation parameters.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """
        dmon_small = getattr(self._parser, 'dmon_small')
        dcrowd = getattr(self._parser, 'dcrowd')
        lcyl = getattr(self._parser, 'lcyl')
        dcyl = getattr(self._parser, 'dcyl')
        return {
            'rEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd), 
                      'lmin': 0,
                      'lmax': 0.5 * lcyl},
            'zEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'thetaEdge': {'bin_size': np.pi / 36,
                          'lmin': -np.pi,
                          'lmax': np.pi},
            'xEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
            'yEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
        }

    def _initialize_collectors(
            self,
            bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes histogram properties and collectors.

        Parameters
        ----------
        bin_edges : Dict[str, Dict[str, Any]]
            A dictionary of bin edge configurations for spatial dimensions.
        """
        for ent in self._entities:
            # 1D histograms
            for direction, bin_type in self._hist1d_bin_types.items():
                edge_key = f"{direction}Edge"
                output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type=bin_type,
                    save_bin_edges=output,
                )
                self._hist1d_props[f'{direction}Hist{ent}'] = {
                    'n_bins': hist_data['n_bins'],
                    'bin_edges': hist_data['bin_edges'],
                    'range': hist_data['range']
                    }
                self._collectors[f'{direction}Hist{ent}'] = \
                    hist_data['collector']
                self._collectors[f'{direction}HistStd{ent}'] = \
                    hist_data['collector_std']

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._hist2d_props[f'{plane}Hist{ent}'] = \
                    {'n_bins': [], 'bin_edges': [], 'range': []}
                for axis in plane:
                    edge_key = f"{axis}Edge"
                    output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                    hist_data = fixedsize_bins(
                        bin_edges[edge_key]['bin_size'],
                        bin_edges[edge_key]['lmin'],
                        bin_edges[edge_key]['lmax'],
                        bin_type='ordinary',
                        save_bin_edges=output,
                    )
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']\
                        .append(hist_data['n_bins'])
                    self._hist2d_props[f'{plane}Hist{ent}']['bin_edges']\
                        .append(hist_data['bin_edges'])
                    self._hist2d_props[f'{plane}Hist{ent}']['range']\
                        .append(hist_data['range'])

                self._collectors[f'{plane}Hist{ent}'] = np.zeros(
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']
                )
            print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _update_histogram(
        self,
        direction: DirectionT,
        histogram_func,
        entity: EntityT,
        axis: int
    ) -> None:
        """
        Updates a one-dimensional histogram for a specific direction.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction ('r', 'z', or 'theta') for the histogram.
        histogram_func : Callable
            The histogram computation function to use (e.g., radial histogram).
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        axis : int
            The spatial axis corresponding to the direction.
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f'{direction}Hist{entity}']['bin_edges'],
            self._hist1d_props[f'{direction}Hist{entity}']['range'],
            axis,
        )
        self._collectors[f'{direction}Hist{entity}'] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)

    def _update_histogram_2d(
        self,
        plane: PlaneT,
        entity: EntityT,
    ) -> None:
        """
        Updates a two-dimensional histogram for a Cartesian plane.

        Parameters
        ----------
        plane : PlaneT
            A Cartesian plane ('xy', 'yz', or 'zx') for the histogram.
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = planar_cartesian_histogram(
            self._atom_groups[entity].positions,
            self._hist2d_props[f'{plane}Hist{entity}']['bin_edges'],
            self._hist2d_props[f'{plane}Hist{entity}']['range'],
            self._hist2d_planes_axes[plane],
        )
        self._collectors[f'{plane}Hist{entity}'] += pos_hist

    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        for ent in self._entities:
            # 1D histograms
            self._update_histogram('r', radial_cyl_histogram, ent, axis=2)
            self._update_histogram('z', axial_histogram, ent, axis=2)
            self._update_histogram('theta', azimuth_cyl_histogram, ent, axis=2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class TransFociCubBugProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *TransFociCub*
    molecular dynamics project to extract specific physical properties.

    Physical Properties Extracted
    -----------------------------
    The particle entity 'Mon' represents both small and large monomers, while 
    the particle entity 'Foci' represents large monomers. For these entities, 
    the following physical properties are extracted:

    For the `Mon` entity:
    - `gyrTMon` : numpy.ndarray
        Radius of gyration of the polymer for each frame.
    - `asphericityTMon` : numpy.ndarray
        Asphericity of the polymer's configuration for each frame.
    - `shapeTMon` : numpy.ndarray
        Shape parameter of the polymer for each frame.
    - `principalTMon` : numpy.ndarray
        Principal axes of the polymer for each frame.

    For the `Foci` entity:
    - `distMatTFoci` : numpy.ndarray
        Pairwise distance matrix for large monomers in the `Foci` entity.
    - `directContactsMatTFoci` : numpy.ndarray
        Direct contact matrix for large monomers based on the cut-off distance.
    - `bondsHistTFoci` : numpy.ndarray
        Histogram of the number of direct contacts for each large monomer.
    - `clustersHistTFoci` : numpy.ndarray
        Histogram of cluster sizes for large monomers in each frame.

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
    - Atoms with `resid 1` in the LAMMPS dump format represent `Mon`(monomers).
    - Coordinates are wrapped and unscaled. The polymer's center of mass is
      recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `TransFociCubBugProber` for a specific simulation:

    >>> prober = TransFociCubBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Foci', 'Mon']

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
        self._parser: TransFociCubInstance = TransFociCub(trajectory, lineage, 'bug')
        self._damping_time = (
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
        )
        self._cluster_cutoff = (
            getattr(self.parser, 'dmon_large') + getattr(self.parser, 'dcrowd')
        ) 
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time
        )

    @property
    def cluster_cutoff(self) -> float:
        """
        Return the cut-off distance for considering two large monomers in
        direct contact.
        """
        return self._cluster_cutoff

    def _define_atom_groups(self) -> None:
        self._atom_groups['Mon'] = \
            self.universe.select_atoms('resid 1')  # small and large monomers
        self._atom_groups['Foci'] = \
            self.universe.select_atoms('type 2')  # large monomers

    def _prepare(self) -> None:
        nmon_large = getattr(self.parser, 'dmon_large')
        self._collectors = {
            'gyrTMon': np.zeros(self._n_frames),
            'asphericityTMon': np.zeros(self._n_frames),
            'shapeTMon': np.zeros(self._n_frames),
            'principalTMon': np.zeros([self._n_frames, 3, 3]),
            'distMatTFoci': np.zeros((self.n_frames, nmon_large, nmon_large)),
            'directContactsMatTFoci': \
                np.zeros((self.n_frames, nmon_large, nmon_large), dtype=int),
            'bondsHistTFoci': np.zeros((self.n_frames, nmon_large), dtype=int),
            'clustersHistTFoci': \
                np.zeros((self.n_frames, nmon_large+1), dtype=int),
        }

    def _probe_frame(self, idx) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['Mon'].radius_of_gyration()
        self._collectors['asphericityTMon'][idx] = \
            self._atom_groups['Mon'].asphericity(wrap=False, unwrap=True)
        self._collectors['shapeTMon'][idx] = \
            self._atom_groups['Mon'].shape_parameter(wrap=False)
        self._collectors['principalTMon'][idx] = \
            self._atom_groups['Mon'].principal_axes(wrap=False)
        dist_mat = \
            clusters.self_dist_array(self._atom_groups['Foci'].positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        np.fill_diagonal(foci_pair_dist, self._atom_groups['Foci'].atoms.ids)
        self._collectors['distMatTFoci'][idx] = foci_pair_dist
        dir_contacts = \
            clusters.find_direct_contacts(dist_mat, self.cluster_cutoff)
        self._collectors['directContactsMatTFoci'][idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        self._collectors['bondsHistTFoci'][idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        self._collectors['clustersHistTFoci'][idx] = clusters_stat


class TransFociCubAllProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the *TransFociCub*
    molecular dynamics project to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> (`Mon` for both small and large monomers, `Dna` for small 
    monomers, `Foci` for large monomers, and `Crd` for crowders), the following 
    spatial distributions are extracted:

    - `rHist<entity>` : numpy.ndarray
        Radial histogram in spherical coordinates.
    - `xyHist<entity>` : numpy.ndarray
        2D histogram in the xy-plane.
    - `yzHist<entity>` : numpy.ndarray
        2D histogram in the yz-plane.
    - `zxHist<entity>` : numpy.ndarray
        2D histogram in the zx-plane.

    Parameters
    ----------
    topology : str
        Path to the topology file of the molecular system.
    trajectory : str
        Path to the trajectory file of the molecular system.
    lineage : {'segment', 'whole'}
        Specifies the type of lineage associated with the trajectory file.
    save_to : str
        Directory where artifacts (e.g., histogram files) will be saved.
    continuous : bool, default False
        Indicates whether the trajectory is part of a continuous sequence.

    Methods
    -------
    _initialize_bin_edges() -> Dict[str, Dict[str, Any]]
        Computes bin edge configurations for histograms.
    _initialize_collectors(bin_edges: Dict[str, Dict[str, Any]]) -> None
        Initializes histogram properties and collectors.
    _update_histogram(direction: DirectionT, histogram_func, entity: EntityT,
                    axis: int) -> None
        Updates a one-dimensional histogram for a specific direction.
    _update_histogram_2d(plane: PlaneT, entity: EntityT) -> None
        Updates a two-dimensional histogram for a specific plane.

    Notes
    -----
    - Atoms with `resid 0` in the LAMMPS dump format represent `Crd`
    (crowders), while atoms with `resid 1` represent `Mon` (monomers).
    - Coordinates are wrapped and unscaled. The center of mass of the `bug`
    group is recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `TransFociCubAllProber` for a specific simulation:

    >>> prober = TransFociCubAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon', 'Dna', 'Foci', 'Crd']
    _hist1d_bin_types: Dict[DirectionT, BinT] = {'r': 'nonnegative'}
    _hist2d_planes_dirs: Dict[PlaneT, DirectionT] = \
        {'xy': 'z', 'yz': 'x', 'zx': 'y'}
    _hist2d_planes_axes: Dict[PlaneT, AxisT] = {'xy': 2, 'yz': 0, 'zx': 1}

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

        self._parser: TransFociCubInstance = \
            TransFociCub(trajectory, lineage, 'all')
        self._damping_time = \
            getattr(self.parser, 'adump') * getattr(self.parser, 'dt')
        self._universe = mda.Universe(
            topology,
            trajectory,
            topology_format='DATA',
            format='LAMMPSDUMP',
            lammps_coordinate_convention='unscaled',
            atom_style="id resid type x y z",
            dt=self.damping_time,
        )

    def _define_atom_groups(self) -> None:
        self._atom_groups['Crd'] = \
            self.universe.select_atoms("resid 0")  # crowders
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # small and large monomers
        self._atom_groups['Dna'] = \
            self.universe.select_atoms("type 1")  # small monomers
        self._atom_groups['Foci'] = \
            self.universe.select_atoms("type 2")  # large monomers

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges = self._initialize_bin_edges()
        self._initialize_collectors(bin_edges)

    def _initialize_bin_edges(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations based on simulation parameters.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """
        dmon_small = getattr(self._parser, 'dmon_small')
        dcrowd = getattr(self._parser, 'dcrowd')
        lcube = getattr(self._parser, 'lcube')
        return {
            'rEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd), 
                      'lmin': 0,
                      'lmax': 0.5 * lcube},
            'zEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
            'xEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
            'yEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcube,
                      'lmax': 0.5 * lcube},
        }

    def _initialize_collectors(
            self,
            bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes histogram properties and collectors.

        Parameters
        ----------
        bin_edges : Dict[str, Dict[str, Any]]
            A dictionary of bin edge configurations for spatial dimensions.
        """
        for ent in self._entities:
            # 1D histograms
            for direction, bin_type in self._hist1d_bin_types.items():
                edge_key = f"{direction}Edge"
                output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type=bin_type,
                    save_bin_edges=output,
                )
                self._hist1d_props[f'{direction}Hist{ent}'] = {
                    'n_bins': hist_data['n_bins'],
                    'bin_edges': hist_data['bin_edges'],
                    'range': hist_data['range']
                    }
                self._collectors[f'{direction}Hist{ent}'] = \
                    hist_data['collector']
                self._collectors[f'{direction}HistStd{ent}'] = \
                    hist_data['collector_std']

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._hist2d_props[f'{plane}Hist{ent}'] = \
                    {'n_bins': [], 'bin_edges': [], 'range': []}
                for axis in plane:
                    edge_key = f"{axis}Edge"
                    output = f"{self.save_to}{self.sim_name}-{edge_key}{ent}"
                    hist_data = fixedsize_bins(
                        bin_edges[edge_key]['bin_size'],
                        bin_edges[edge_key]['lmin'],
                        bin_edges[edge_key]['lmax'],
                        bin_type='ordinary',
                        save_bin_edges=output,
                    )
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']\
                        .append(hist_data['n_bins'])
                    self._hist2d_props[f'{plane}Hist{ent}']['bin_edges']\
                        .append(hist_data['bin_edges'])
                    self._hist2d_props[f'{plane}Hist{ent}']['range']\
                        .append(hist_data['range'])

                self._collectors[f'{plane}Hist{ent}'] = np.zeros(
                    self._hist2d_props[f'{plane}Hist{ent}']['n_bins']
                )
            print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _update_histogram(
        self,
        direction: DirectionT,
        histogram_func,
        entity: EntityT,
    ) -> None:
        """
        Updates a one-dimensional histogram for a specific direction.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction ('r', 'z', or 'theta') for the histogram.
        histogram_func : Callable
            The histogram computation function to use (e.g., radial histogram).
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f'{direction}Hist{entity}']['bin_edges'],
            self._hist1d_props[f'{direction}Hist{entity}']['range']
        )
        self._collectors[f'{direction}Hist{entity}'] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)

    def _update_histogram_2d(
        self,
        plane: PlaneT,
        entity: EntityT,
    ) -> None:
        """
        Updates a two-dimensional histogram for a Cartesian plane.

        Parameters
        ----------
        plane : PlaneT
            A Cartesian plane ('xy', 'yz', or 'zx') for the histogram.
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = planar_cartesian_histogram(
            self._atom_groups[entity].positions,
            self._hist2d_props[f'{plane}Hist{entity}']['bin_edges'],
            self._hist2d_props[f'{plane}Hist{entity}']['range'],
            self._hist2d_planes_axes[plane],
        )
        self._collectors[f'{plane}Hist{entity}'] += pos_hist

    def _probe_frame(self, idx: int) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        for ent in self._entities:
            # 1D histograms
            self._update_histogram('r', radial_histogram, ent, axis=2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)

class SumRuleHeteroLinearCubBugProber(ProberBase):


class SumRuleHeteroLinearCubAllProber(ProberBase):


class SumRuleHeteroRingCubBugProber(ProberBase):


class SumRuleHeteroRingCubAllProber(ProberBase):


class HnsCylNucleoidProber(ProberBase):


class HnsCylAllProber(ProberBase):


class HnsCubNucleoidProber(ProberBase):


class HnsCubAllProber(ProberBase):