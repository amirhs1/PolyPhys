import warnings
import pathlib
from abc import ABC, abstractmethod
from typing import (Dict, Any, Tuple, Literal, Optional, List, Callable,
                    ClassVar)
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mda_dist
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
    PrimitiveLineageT,
    ParserType,
    SumRuleCylInstance,
    TwoMonDepCubInstance,
    TransFociCylInstance,
    TransFociCubT,
    SumRuleCubHeteroRingT,
    SumRuleCubHeteroLinearT,
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
    _hist1d_bin_types: ClassVar[Dict[DirectionT, BinT]] = \
        {'r': 'nonnegative', 'z': 'ordinary', 'theta': 'periodic'}

    @abstractmethod
    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations for calculating one dimensional histograms based on simulation parameters.

        Must be implemented in subclasses to supply cylindrical-specific bin configurations.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """

    def _initialize_hist_collectors_1d(
        self,
        entity: EntityT,
        bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes collectors for 1D histograms in cylindrical coordinates.

        Parameters
        ----------
        entity : EntityT
            Particle entity (e.g., `Mon` or `Crd`).
        bin_edges : Dict[str, Dict[str, Any]]
            Bin edge configurations for cylindrical directions.
        """
        for direction, bin_type in self._hist1d_bin_types.items():
            edge_key = f"{direction}Edge"
            output = f"{self.save_to}{self.sim_name}-{edge_key}{entity}"
            hist_data = fixedsize_bins(
                bin_edges[edge_key]['bin_size'],
                bin_edges[edge_key]['lmin'],
                bin_edges[edge_key]['lmax'],
                bin_type=bin_type,
                save_bin_edges=output
            )
            self._hist1d_props[f"{direction}Hist{entity}"] = {
                'n_bins': hist_data['n_bins'],
                'bin_edges': hist_data['bin_edges'],
                'range': hist_data['range'],
            }
            self._collectors[f"{direction}Hist{entity}"] = \
                hist_data['collector']
            self._collectors[f'{direction}HistStd{entity}'] = \
                hist_data['collector_std']

    def _update_histogram_1d(
        self,
        direction: DirectionT,
        histogram_func: Callable,
        entity: EntityT,
        axis: AxisT
    ) -> None:
        """
        Updates a 1D histogram in cylindrical coordinates.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction (`r`, `z`, or `theta`).
        histogram_func : Callable
            Function to compute the histogram.
        entity : EntityT
            Particle entity (e.g., `Mon` or `Crd`).
        axis : AxisT
            Spatial axis corresponding to the direction.
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f"{direction}Hist{entity}"]['bin_edges'],
            self._hist1d_props[f"{direction}Hist{entity}"]['range'],
            axis
        )
        self._collectors[f"{direction}Hist{entity}"] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)


class SphericalHistogramMixin(ABC):
    """
    Mixin for 1D histograms in spherical coordinate systems. Currently, only handles histograms for `r`.
    """
    _hist1d_bin_types: ClassVar[Dict[DirectionT, BinT]] = {'r': 'nonnegative'}

    @abstractmethod
    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations for calculating one dimensional histograms based on simulation parameters.

        Must be implemented in subclasses to supply spherical-specific bin configurations.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """

    def _initialize_hist_collectors_1d(
        self,
        entity: EntityT,
        bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes collectors for 1D histograms in spherical coordinates.

        Parameters
        ----------
        entity : EntityT
            Particle entity (e.g., `Mon` or `Crd`).
        bin_edges : Dict[str, Dict[str, Any]]
            Bin edge configurations for spherical directions.
        """
        for direction, bin_type in self._hist1d_bin_types.items():
            edge_key = f"{direction}Edge"
            output = f"{self.save_to}{self.sim_name}-{edge_key}{entity}"
            hist_data = fixedsize_bins(
                bin_edges[edge_key]['bin_size'],
                bin_edges[edge_key]['lmin'],
                bin_edges[edge_key]['lmax'],
                bin_type=bin_type,
                save_bin_edges=output
            )
            self._hist1d_props[f"{direction}Hist{entity}"] = {
                'n_bins': hist_data['n_bins'],
                'bin_edges': hist_data['bin_edges'],
                'range': hist_data['range'],
            }
            self._collectors[f"{direction}Hist{entity}"] = \
                hist_data['collector']
            self._collectors[f'{direction}HistStd{entity}'] = \
                hist_data['collector_std']

    def _update_histogram_1d(
        self,
        direction: DirectionT,
        histogram_func: Callable,
        entity: EntityT,
    ) -> None:
        """
        Updates a 1D histogram in spherical coordinates.

        Parameters
        ----------
        direction : DirectionT
            The spatial direction (`r`).
        histogram_func : Callable
            Function to compute the histogram.
        entity : EntityT
            Particle entity (e.g., `Mon` or `Crd`).
        """
        pos_hist, _ = histogram_func(
            self._atom_groups[entity].positions,
            self._hist1d_props[f"{direction}Hist{entity}"]['bin_edges'],
            self._hist1d_props[f"{direction}Hist{entity}"]['range'],
        )
        self._collectors[f"{direction}Hist{entity}"] += pos_hist
        self._collectors[f'{direction}HistStd{entity}'] += np.square(pos_hist)


class CartesianHistogram2DMixin(ABC):
    """
    Mixin for 2D histograms in Cartesian coordinate systems. Handles histograms
    for `xy`, `yz`, and `zx` planes.
    """
    _hist2d_planes_dirs: ClassVar[Dict[PlaneT, DirectionT]] = \
        {'xy': 'z', 'yz': 'x', 'zx': 'y'}
    _hist2d_planes_axes: ClassVar[Dict[PlaneT, AxisT]] = {'xy': 2, 'yz': 0, 'zx': 1}

    @abstractmethod
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations based on simulation parameters.

        Must be implemented in subclasses to supply Cartesian-specific bin configurations.

        Returns
        -------
        Dict[str, Dict[str, Any]]
            A dictionary of bin edge settings, including bin size, min, and max
            for each spatial dimension.
        """

    def _initialize_hist_collectors_2d(
        self,
        entity: EntityT,
        bin_edges: Dict[str, Dict[str, Any]]
    ) -> None:
        """
        Initializes collectors for 2D histograms in Cartesian coordinates.

        Parameters
        ----------
        entity : EntityT
            Particle entity (e.g., `Mon` or `Crd`).
        bin_edges : Dict[str, Dict[str, Any]]
            Bin edge configurations for Cartesian planes.
        """
        for plane in self._hist2d_planes_dirs:
            self._hist2d_props[f'{plane}Hist{entity}'] = \
                    {'n_bins': [], 'bin_edges': [], 'range': []}
            for axis in plane:
                edge_key = f"{axis}Edge"
                output = f"{self.save_to}{self.sim_name}-{edge_key}{entity}"
                hist_data = fixedsize_bins(
                    bin_edges[edge_key]['bin_size'],
                    bin_edges[edge_key]['lmin'],
                    bin_edges[edge_key]['lmax'],
                    bin_type='ordinary',
                    save_bin_edges=output,
                )
                self._hist2d_props[f'{plane}Hist{entity}']['n_bins']\
                    .append(hist_data['n_bins'])
                self._hist2d_props[f'{plane}Hist{entity}']['bin_edges']\
                    .append(hist_data['bin_edges'])
                self._hist2d_props[f'{plane}Hist{entity}']['range']\
                    .append(hist_data['range'])

            self._collectors[f'{plane}Hist{entity}'] = np.zeros(
                self._hist2d_props[f'{plane}Hist{entity}']['n_bins']
            )

    def _update_histogram_2d(self, plane: PlaneT, entity: EntityT) -> None:
        """
        Updates a 2D histogram in Cartesian coordinates.

        Parameters
        ----------
        plane : PlaneT
            A Cartesian plane ('xy', 'yz', or 'zx') for the histogram.
        entity : EntityT
            The type of entity ('Mon' for monomers, 'Crd' for crowders).
        """
        pos_hist, _ = planar_cartesian_histogram(
            self._atom_groups[entity].positions,
            self._hist2d_props[f"{plane}Hist{entity}"]['bin_edges'],
            self._hist2d_props[f"{plane}Hist{entity}"]['range'],
            self._hist2d_planes_axes[plane],
        )
        self._collectors[f"{plane}Hist{entity}"] += pos_hist


class GyrationDataMixin:
    """
    Mixin class for collecting gyration-related data for the 'Mon' entity 
    in molecular dynamics simulations.

    This mixin provides functionality to compute and store physical properties 
    related to the gyration of a polymer, including radius of gyration, 
    asphericity, shape parameter, and principal axes.

    Methods
    -------
    _collect_gyration_data(idx: int) -> None
        Collects and stores gyration-related data for the current frame.

    Requirements
    ------------
    - The class that uses this mixin must have the following attributes:
        - `_atom_groups`: Dict[str, MDAnalysis.AtomGroup]
            Dictionary mapping particle entities to MDAnalysis AtomGroups.
        - `_collectors`: Dict[str, Any]
            Dictionary for storing computed physical properties.
    """
    def _collect_gyration_data(self, idx: int, unwrap: bool = False) -> None:
        """
        Collects gyration-related data for the 'Mon' entity.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        unwrap: bool, optional
             If True, compounds will be unwrapped before computing their centers. default to False

        Updates
        -------
        - `gyrTMon`: Radius of gyration.
        - `asphericityTMon`: Asphericity of the polymer.
        - `shapeTMon`: Shape parameter of the polymer.
        - `principalTMon`: Principal axes of the polymer.
        """
        self._collectors['gyrTMon'][idx] = \
            self._atom_groups['Mon'].radius_of_gyration()
        self._collectors['asphericityTMon'][idx] = \
            self._atom_groups['Mon'].asphericity(wrap=False, unwrap=unwrap)
        self._collectors['shapeTMon'][idx] = \
            self._atom_groups['Mon'].shape_parameter(wrap=False)
        self._collectors['principalTMon'][idx] = \
            self._atom_groups['Mon'].principal_axes(wrap=False)


class BiaxialConfinementSizeMixin:
    """
    Mixin class for collecting biaxial confinement size data in molecular 
    dynamics simulations along a specified longitudinal axis.

    This mixin provides functionality to compute and store transverse size 
    and farthermost distance (Feret's statistical diameter) for a particle
    entity perpendicular to the specified longitudinal axis.

    Methods
    -------
    _collect_biaxial_confinement_size(idx: int, axis: int) -> None
        Collects and stores biaxial confinement size data for the current 
        frame along the specified axis.

    Requirements
    ------------
    - The class that uses this mixin must have the following attributes:
        - `_atom_groups`: Dict[str, Any]
            Dictionary mapping particle entities to atom groups.
        - `_collectors`: Dict[str, Any]
            Dictionary for storing computed physical properties.
    """
    def _collect_biaxial_confinement_size(self, idx: int, axis: int) -> None:
        """
        Collects biaxial confinement size data for the `Mon` entity along the 
        specified longitudinal axis.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        axis : int
            The longitudinal axis (0 for x, 1 for y, 2 for z) along which the 
            cylinder or cylinder-like confinement is located.

        Updates
        -------
        - `transSizeTMon`: Transverse size in the plan perpendicular to `axis` 
        - `fsdTMon`: Farthermost distance of positions along the `axis`.
        """
        self._collectors['transSizeTMon'][idx] = \
            transverse_size(self._atom_groups['Mon'].positions, axis=axis)
        self._collectors['fsdTMon'][idx] = \
            fsd(self._atom_groups['Mon'].positions, axis=axis)


class FociContactDataMixin:
    """
    Mixin class for collecting contact-related data for the 'Foci' entity 
    in molecular dynamics simulations.

    This mixin provides functionality to compute and store contact-related 
    properties, including pairwise distances, direct contacts, bond statistics, 
    and cluster sizes for the 'Foci' entity.

    Methods
    -------
    _collect_foci_contact_data(idx: int) -> None
        Collects and stores contact-related data for the current frame.

    Requirements
    ------------
    - The class that uses this mixin must have the following attributes:
        - `_atom_groups`: Dict[str, MDAnalysis.AtomGroup]
            Dictionary mapping particle entities to MDAnalysis AtomGroups.
        - `_collectors`: Dict[str, Any]
            Dictionary for storing computed physical properties.
        - `cluster_cutoff`: float
            Distance cutoff for defining direct contacts between large monomers.
    """
    def _collect_foci_contact_data(self, idx: int) -> None:
        """
        Collects contact-related data for the 'Foci' entity.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.

        Updates
        -------
        - `distMatTFoci`: Pairwise distance matrix for large monomers.
        - `directContactsMatTFoci`: Direct contact matrix.
        - `bondsHistTFoci`: Histogram of direct contacts per monomer.
        - `clustersHistTFoci`: Histogram of cluster sizes.
        """
        dist_mat = \
            clusters.self_dist_array(self._atom_groups['Foci'].positions)
        foci_pair_dist = np.triu(dist_mat, 1)
        np.fill_diagonal(foci_pair_dist, self._atom_groups['Foci'].atoms.ids)
        self._collectors['distMatTFoci'][idx] = foci_pair_dist
        dir_contacts = clusters.find_direct_contacts(dist_mat, self.cluster_cutoff)
        self._collectors['directContactsMatTFoci'][idx] = dir_contacts
        bonds_stat = clusters.count_foci_bonds(dir_contacts)
        self._collectors['bondsHistTFoci'][idx] = bonds_stat
        contacts = clusters.generate_contact_matrix(dir_contacts)
        clusters_stat = clusters.count_foci_clusters(contacts)
        self._collectors['clustersHistTFoci'][idx] = clusters_stat


class HnsBindingMixin:
    """
    Mixin for calculating the binding statistics of H-NS proteins in a polymer system.

    This mixin enables calculation of H-NS binding behavior, handling
    statistics for free, bridged, and dangling H-NS proteins. It also computes
    genomic distances between monomers bridged by H-NS cores and classifies
    them as cis or trans based on a distance threshold.

    Requirements
    ------------
    - The class using this mixin must define the following attributes:
        - `_atom_groups`: Dict[str, MDAnalysis.AtomGroup]
            Atom groups for 'Mon' (monomers) and 'Hns' (H-NS proteins).
        - `collectors`: Dict[str, List or np.ndarray]
            Dictionary to store collected data.
        - `parser`: Instance of a parser providing simulation metadata.
        - `cis_threshold`: int
            Distance threshold for cis-binding classification.
        - `r_cutoff`: float
            Distance cutoff for direct monomer-patch contact.

    Methods
    -------
    _collect_hns_binding(idx: int) -> None
        Collects H-NS binding statistics for a single trajectory frame.
    """
    
    def _collect_hns_binding(self, idx: int) -> None:
        """
        Collects H-NS binding statistics for a single trajectory frame.

        Parameters
        ----------
        idx : int
            Index of the trajectory frame to analyze.

        Updates
        -------
        Collectors:
        - `distMatTMonHnsPatch`: Distance matrix between monomers and H-NS patches.
        - `nBoundTHnsPatch`, `nEngagedTHnsPatch`, `nFreeTHnsPatch`: Patch-level stats.
        - `nFreeTHnsCore`, `nBridgeTHnsCore`, `nDangleTHnsCore`: Core-level stats.
        - `nCisTHnsCore`, `nTransTHnsCore`: Cis/trans-binding stats.
        - Histogram of genomic distances: Updated with loop sizes.
        """
        # Calculate distance matrix and direct contacts
        pair_dist_mat = mda_dist.distance_array(
            self._atom_groups['Mon'], 
            self._atom_groups['Hns']
        )
        self.collectors['distMatTMonHnsPatch'].append(pair_dist_mat)
        dir_cont_m_hpatch = clusters.find_direct_contacts(
            pair_dist_mat, self.r_cutoff, inclusive=False
        )
        
        # Dimensions and histogram setup
        n_mon, n_hpatch = dir_cont_m_hpatch.shape
        n_hcore = n_hpatch // 2
        max_gen_distance = {'linear': n_mon, 'ring': n_mon // 2 + 1}
        hist = np.zeros(max_gen_distance[self.parser.topology], dtype=int)

        # Patch-level statistics
        d_cont_per_hpatch = np.sum(pair_dist_mat, axis=0)
        self.collectors['nBoundTHnsPatch'][idx] = np.sum(pair_dist_mat)
        n_engaged_hpatch = np.count_nonzero(d_cont_per_hpatch)
        self.collectors['nEngagedTHnsPatch'][idx] = n_engaged_hpatch
        self.collectors['nFreeTHnsPatch'][idx] = n_hpatch - n_engaged_hpatch

        # Core-level statistics
        free_cores = np.sum(
            (d_cont_per_hpatch[0::2] == 0) & (d_cont_per_hpatch[1::2] == 0)
        )
        bridged_cores = np.sum(
            (d_cont_per_hpatch[0::2] > 0) & (d_cont_per_hpatch[1::2] > 0)
        )
        dangling_cores = n_hcore - free_cores - bridged_cores
        self.collectors['nFreeTHnsCore'][idx] = free_cores
        self.collectors['nBridgeTHnsCore'][idx] = bridged_cores
        self.collectors['nDangleTHnsCore'][idx] = dangling_cores

        # Enforce single-patch binding and calculate genomic distances
        single_patch_dir_contact = \
            clusters.enforce_single_patch_dir_contact(pair_dist_mat)
        cont_m_hpatch = clusters.generate_contact_matrix(single_patch_dir_contact)
        cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
        cont_m_m = np.matmul(cont_m_hcore, cont_m_hcore.T)
        cont_m_m_nonzeros = np.array(np.where(np.triu(cont_m_m, 1) > 0)).T

        genomic_distances = clusters.genomic_distance(
            cont_m_m_nonzeros[:, 0],
            cont_m_m_nonzeros[:, 1],
            self.parser.topology,
            n_mon
        )
        m_m_gen_dist = np.column_stack((cont_m_m_nonzeros, genomic_distances))

        # Cis/trans classification
        cis_bindings = np.count_nonzero(
            (m_m_gen_dist[:, 2] > 0) & (m_m_gen_dist[:, 2] <= self.cis_threshold)
        )
        trans_bindings = np.count_nonzero(m_m_gen_dist[:, 2] > self.cis_threshold)
        self.collectors['nCisTHnsCore'][idx] = cis_bindings
        self.collectors['nTransTHnsCore'][idx] = trans_bindings
        np.add.at(hist, m_m_gen_dist[:, 2], 1)


class HnsBindingMixin:
    """
    Calculate the binding statistics of H-NS proteins in a polymer system.

    This function evaluates the binding behavior of H-NS proteins, which
    consist of a core and two patches capable of binding monomers in the
    polymer. It computes statistics for free, bridged, and dangling H-NS
    proteins, and determines genomic distances between monomers bridged by H-NS
    cores.

    The genomic distance is classified as 'cis' or 'trans' based on the
    `cis_threshold`, and a histogram of genomic distances (loop sizes) is
    updated accordingly.

    Parameters
    ----------
    direct_contact : np.ndarray
        A binary matrix of shape `(n_mon, n_hpatch)` where:
        - `n_mon` is the number of monomers.
        - `n_hpatch` is the number of H-NS patches (twice the number of H-NS cores).
        Each element `(i, j)` is `1` if monomer `i` is in contact with patch
        `j`, and `0` otherwise.

    topology : TopologyT
        The topology of the polymer, either 'linear' or 'ring'.

    cis_threshold : int, optional
        The genomic distance threshold (in bonds) below which an H-NS binding
        is classified as 'cis'. Default is `4`.

    Returns
    -------
    Tuple[Dict[str, int], np.ndarray]
        - `stats`: A dictionary containing various binding statistics:
          - 'nBoundHnsPatch': Total number of monomer-patch contacts.
          - 'nEngagedHnsPatch': Total number of engaged patches.
          - 'nFreeHnsPatch': Total number of free patches.
          - 'nFreeHnsCore': Total number of free H-NS cores.
          - 'nBridgeHnsCore': Total number of bridging H-NS cores.
          - 'nDangleHnsCore': Total number of dangling H-NS cores.
          - 'nCisHnsCore': Total number of cis-binding H-NS cores.
          - 'nTransHnsCore': Total number of trans-binding H-NS cores.

        - `hist`: A 1D array representing a histogram of genomic distances
          (loop sizes).

    Raises
    ------
    ValueError
        If the provided topology is invalid ('linear' or 'ring' are the only
        valid options).

    Notes
    -----
    - An H-NS core is considered free if neither of its two patches is bound.
    - A core is considered bridging if both patches are bound to different
      monomers.
    - A core is considered dangling if only one of its patches is bound.
    - Genomic distances are computed between monomer pairs bridged by H-NS cores.
    - For a 'linear' topology, the genomic distance is the absolute difference
      between monomer indices.
    - For a 'ring' topology, the genomic distance is the minimum of the direct
      distance and the distance around the circle.
    """
    def _collect_hns_binding(self, idx: int) -> None:
        pair_dist_mat = mda_dist.distance_array(
            self._atom_groups['Mon'], 
            self._atom_groups['Hns']
        )
        self.collectors['distMatTMonHnsPatch'].append(pair_dist_mat)
        dir_cont_m_hpatch = clusters.find_direct_contacts(
            pair_dist_mat, self.r_cutoff, inclusive=False
        )      
        n_mon, n_hpatch = dir_cont_m_hpatch.shape
        n_hcore = n_hpatch // 2

        max_gen_distance = {'linear': n_mon, 'ring': n_mon // 2 + 1}
        hist = np.zeros(max_gen_distance[self.parser.topology], dtype=int)

        # Compute patch-level and core-level stats
        d_cont_per_hpatch = np.sum(pair_dist_mat, axis=0)
        n_engaged_hpatch = np.count_nonzero(d_cont_per_hpatch)
        self.collectors['nBoundTHnsPatch'][idx] = np.sum(pair_dist_mat)
        self.collectors['nEngagedTHnsPatch'][idx] = n_engaged_hpatch 
        self.collectors['nFreeTHnsPatch'][idx] = n_hpatch - n_engaged_hpatch

        free_cores = \
            np.sum(
                (d_cont_per_hpatch[0::2] == 0) & (d_cont_per_hpatch[1::2] == 0)
                )
        bridged_cores = \
            np.sum(
                (d_cont_per_hpatch[0::2] > 0) & (d_cont_per_hpatch[1::2] > 0)
                )
        dangling_cores = n_hcore - free_cores - bridged_cores

        self.collectors['nFreeTHnsCore'][idx] = free_cores
        self.collectors['nBridgeTHnsCore'][idx] = bridged_cores
        self.collectors['nDangleTHnsCore'][idx] = dangling_cores

        # Enforce single-patch binding and calculate genomic distances
        single_patch_dir_contact = \
            clusters.enforce_single_patch_dir_contact(pair_dist_mat)
        cont_m_hpatch = \
            clusters.generate_contact_matrix(single_patch_dir_contact)
        cont_m_hcore = np.logical_or(cont_m_hpatch[:, ::2], cont_m_hpatch[:, 1::2])
        cont_m_m = np.matmul(cont_m_hcore, cont_m_hcore.T)
        cont_m_m_nonzeros = np.array(np.where(np.triu(cont_m_m, 1) > 0)).T
        # Calculate genomic distances for monomer pairs bridged by an H-NS protein:
        genomic_distances = clusters.genomic_distance(
            cont_m_m_nonzeros[:, 0],
            cont_m_m_nonzeros[:, 1],
            self.parser.topology,
            n_mon
        )
        # A 2D array with shape `(n_pairs, 3)`: Column 0 is the index of the first
        # monomer in the pair. Column 1 is the index of the second monomer in the
        # pair. And, column 3 is the genomic distance between the two monomers:
        m_m_gen_dist = np.column_stack((cont_m_m_nonzeros, genomic_distances))

        # Classify bindings and update histogram
        cis_bindings = np.count_nonzero(
            (m_m_gen_dist[:, 2] > 0) & 
            (m_m_gen_dist[:, 2] <= self.cis_threshold)
         )
        trans_bindings = np.count_nonzero(
            m_m_gen_dist[:, 2] > self.cis_threshold)
        self.collectors['nCisTHnsCore'][idx] = cis_bindings
        self.collectors['nTransTHnsCore'][idx] = trans_bindings
        np.add.at(hist, m_m_gen_dist[:, 2], 1)


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
    _entities: ClassVar[Optional[List[EntityT]]] = None

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        if lineage == "segment" and not continuous:
            warnings.warn(
                f"Lineage is '{lineage}' and 'continuous' is False. "
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
        self._atom_groups: Optional[Dict[EntityT, mda.AtomGroup]] = None
        self._define_atom_groups()
        self._collectors: Dict[str, Any] = {}
        self._prepare()

    @property
    def entities(self) -> List[EntityT]:
        """
        Returns the list particle entities in  the molecular dynamics system.
        """
        if self._entities is None:
            raise AttributeError("'_entities' has not been initialized.")
        return self._entities

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
    def lineage(self) -> PrimitiveLineageT:
        """
        Returns the lineage of the trajectory file associated with the
        molecular dynamics system.
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
    def atom_groups(self) -> Dict[EntityT, mda.AtomGroup]:
        """
        Returns a dictionary of atom groups and their associated MDAnalysis
        AtomGroup instances.
        """
        if self._atom_groups is None:
            raise AttributeError("'_atom_groups' has not been initialized.")
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
        lineage: PrimitiveLineageT,
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
            'gyrTMon': np.zeros(self.n_frames),
            'dxTMon': np.zeros(self.n_frames),
            'dyTMon': np.zeros(self.n_frames),
            'dzTMon': np.zeros(self.n_frames)
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


class TwoMonDepCubAllProber(
    SphericalHistogramMixin,
    CartesianHistogram2DMixin,
    ProberBase
):
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

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
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
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
        }

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
            self._update_histogram_1d('r', radial_histogram, ent)
            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)

class SumRuleCylBugProber(
    BiaxialConfinementSizeMixin,
    GyrationDataMixin,
    ProberBase):
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
        lineage: PrimitiveLineageT,
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
            'transSizeTMon': np.zeros(self.n_frames),
            'fsdTMon': np.zeros(self.n_frames),
            'gyrTMon': np.zeros(self.n_frames),
            'rfloryTMon': np.zeros(self.n_frames),
            'asphericityTMon': np.zeros(self.n_frames),
            'shapeTMon': np.zeros(self.n_frames),
            'principalTMon': np.zeros([self.n_frames, 3, 3])
        }

    def _probe_frame(self, idx) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collect_biaxial_confinement_size(idx, axis=2)
        self._collectors['rfloryTMon'][idx] = \
            end_to_end(self._atom_groups['Mon'].positions)
        self._collect_gyration_data(idx, unwrap=False)


class SumRuleCylAllProber(
    CylindricalHistogramMixIn,
    CartesianHistogram2DMixin,
    ProberBase
):
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

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
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
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations for calculating one dimentional histograms based on simulation parameters.

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
            'rEdge': {'bin_size': 0.1 * min(dmon, dcrowd),
                      'lmin': 0,
                      'lmax': 0.5 * dcyl},
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'thetaEdge': {'bin_size': np.pi / 36,
                          'lmin': -np.pi,
                          'lmax': np.pi},
        }
    
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
        """
        Computes bin edge configurations for calculating two dimentional histograms based on simulation parameters.

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
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'xEdge': {'bin_size': 0.1 * min(dmon, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
            'yEdge': {'bin_size': 0.1 * min(dmon, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
        }

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
            self._update_histogram_1d('r', radial_cyl_histogram, ent, 2)
            self._update_histogram_1d('z', axial_histogram, ent, 2)
            self._update_histogram_1d('theta', azimuth_cyl_histogram, ent, 2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class TransFociCylBugProber(
    BiaxialConfinementSizeMixin,
    GyrationDataMixin,
    FociContactDataMixin,
    ProberBase
):
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
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        self._parser: TransFociCylInstance = TransFociCyl(trajectory, lineage, 'bug')
        self._damping_time = (
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
        )
        self._cluster_cutoff: float = (
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
            'transSizeTMon': np.zeros(self.n_frames),
            'fsdTMon': np.zeros(self.n_frames),
            'gyrTMon': np.zeros(self.n_frames),
            'asphericityTMon': np.zeros(self.n_frames),
            'shapeTMon': np.zeros(self.n_frames),
            'principalTMon': np.zeros([self.n_frames, 3, 3]),
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
        self._collect_biaxial_confinement_size(idx, axis=2)
        self._collect_gyration_data(idx, unwrap=False)
        self._collect_foci_contact_data(idx)
        

class TransFociCylAllProber(
    CylindricalHistogramMixIn,
    CartesianHistogram2DMixin,
    ProberBase
):
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

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
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
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
                      'lmax': 0.5 * dcyl},
            'zEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'thetaEdge': {'bin_size': np.pi / 36,
                          'lmin': -np.pi,
                          'lmax': np.pi}
        }
    
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
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
            'zEdge': {'bin_size': 0.5 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'xEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
            'yEdge': {'bin_size': 0.1 * min(dmon_small, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl}
        }

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
            self._update_histogram_1d('r', radial_cyl_histogram, ent, 2)
            self._update_histogram_1d('z', axial_histogram, ent, 2)
            self._update_histogram_1d('theta', azimuth_cyl_histogram, ent, 2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class SumRuleCubHeteroBugProberBase(
    GyrationDataMixin,
    FociContactDataMixin,
    ProberBase
):
    """
    Base class for probing simulations of the LAMMPS 'bug' atom group in the *TransFociCub*, *SumRuleCubHeteroRing*, and *SumRuleCubHeteroLinear* molecular dynamics projects to extract specific physical properties.

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
    """
    parser_cls: Optional[ParserType] = None  # Define in subclasses

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        if self.parser_cls is None:
            raise NotImplementedError(
                f"{self.__class__.__name__} must define the `parser_cls` attribute."
            )
        self._parser = self.parser_cls(trajectory, lineage, 'bug')
        self._damping_time = (
            getattr(self.parser, 'bdump') * getattr(self.parser, 'dt')
        )
        self._cluster_cutoff: float = (
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
            'gyrTMon': np.zeros(self.n_frames),
            'asphericityTMon': np.zeros(self.n_frames),
            'shapeTMon': np.zeros(self.n_frames),
            'principalTMon': np.zeros([self.n_frames, 3, 3]),
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
        self._collect_gyration_data(idx, unwrap=True)
        self._collect_foci_contact_data(idx)

class SumRuleCubHeteroAllProberBase(
    SphericalHistogramMixin,
    CartesianHistogram2DMixin,
    ProberBase
):
    """
    Base class for probing simulations of the LAMMPS 'all' atom group in the in the *TransFociCub*, *SumRuleCubHeteroRing*, and *SumRuleCubHeteroLinear* molecular dynamics projects to extract spatial distributions of particles.

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
    parser_cls: Optional[ParserType] = None  # Define in subclasses
    
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        if self.parser_cls is None:
            raise NotImplementedError(
                f"{self.__class__.__name__} must define the `parser_cls` attribute."
            )
        self._parser = self.parser_cls(trajectory, lineage, 'all')
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
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
        }
    
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
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
            self._update_histogram_1d('r', radial_histogram, ent)
            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class TransFociCubBugProber(SumRuleCubHeteroBugProberBase):
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
    _entities = ['Mon', 'Foci']
    parser_cls: TransFociCubT = TransFociCub 

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)


class TransFociCubAllProber(SumRuleCubHeteroAllProberBase):
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
    _entities = ['Mon', 'Dna', 'Foci', 'Crd']
    parser_cls: TransFociCubT = TransFociCub
    
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)


class SumRuleCubHeteroLinearBugProber(SumRuleCubHeteroBugProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the
    *SumRuleCubHeteroLinear* molecular dynamics project to extract specific
    physical properties.

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
    Creating an instance of `SumRuleCubHeteroLinearBugProber` for a specific
    simulation:

    >>> prober = SumRuleCubHeteroLinearBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities= ['Mon', 'Foci']
    parser_cls: SumRuleCubHeteroLinearT = SumRuleCubHeteroLinear     

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)


class SumRuleCubHeteroLinearAllProber(SumRuleCubHeteroAllProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the
    *SumRuleCubHeteroLinear* molecular dynamics project to extract spatial
    distributions of particles.

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
    Creating an instance of `SumRuleCubHeteroLinearAllProber` for a specific
    simulation:

    >>> prober = SumRuleCubHeteroLinearAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Dna', 'Foci', 'Crd']
    parser_cls: SumRuleCubHeteroLinearT = SumRuleCubHeteroLinear
    
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)

class SumRuleHeteroCubRingBugProber(SumRuleCubHeteroBugProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *SumRuleCubHeteroRing* molecular dynamics project to extract specific
    physical properties.

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
    Creating an instance of `SumRuleHeteroCubRingBugProber` for a specific
    simulation:

    >>> prober = SumRuleHeteroRingCubBugProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Foci']
    parser_cls: SumRuleCubHeteroRingT = SumRuleCubHeteroRing 

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)


class SumRuleCubHeteroRingAllProber(ProberBase):
    """
    Probes simulations of the LAMMPS 'all' atom group in the
    *SumRuleCubHeteroRing* molecular dynamics project to extract spatial
    distributions of particles.

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
    Creating an instance of `SumRuleCubHeteroRingAllProber` for a specific
    simulation:

    >>> prober = SumRuleCubHeteroRingAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Dna', 'Foci', 'Crd']
    parser_cls: SumRuleCubHeteroRingT = SumRuleCubHeteroRing
    
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)


class HnsCylNucleoidProber(
    BiaxialConfinementSizeMixin,
    GyrationDataMixin,
    ProberBase):
    """
    Probes simulations of the LAMMPS 'bug' atom group in the *HnsCyl*
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
    Creating an instance of `HnsCylNucleoidProber` for a specific simulation:

    >>> prober = HnsCylNucleoidProber(
    ...     topology="topology.bug.data",
    ...     trajectory="trajectory.bug.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities: List[EntityT] = ['Mon', 'Hns']

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        self._parser: HnsCylInstance = HnsCyl(trajectory, lineage, 'nucleoid')
        self._damping_time = \
            getattr(self.parser, 'ndump') * getattr(self.parser, 'dt')
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
        self._atom_groups['Hns'] = \
            self.universe.select_atoms('type 2')  # monomers

    @property
    def n_bonds(self) -> int:
        """
        Returns the number of bonds between monomers in a polymer.
        """
        return self._n_bonds

    def _prepare(self) -> None:
        self.cis_threshold = 4
        dist_m_hpatch = []
        lj_cut = 2**(1/6)
        dmon = getattr(self.parser, 'dmon')
        self.dhns_patch = getattr(self.parser, 'dhns_patch') 
        self.r_cutoff = np.round(0.5 * lj_cut * (dmon + dhns_patch), 3)
        self._n_bonds = len(self._atom_groups['Mon'].bonds.indices)
        self._collectors = {
            'transSizeTMon': np.zeros(self.n_frames),
            'fsdTMon': np.zeros(self.n_frames),
            'gyrTMon': np.zeros(self.n_frames),
            'asphericityTMon': np.zeros(self.n_frames),
            'shapeTMon': np.zeros(self.n_frames),
            'principalTMon': np.zeros([self.n_frames, 3, 3]),
            'bondLengthVecMon': np.zeros((self.n_bonds, 1), dtype=np.floating),
            'bondCosineCorrVecMon': np.zeros(self.n_bonds, dtype=np.floating),
            'distMatTMonHnsPatch': [],
            'nBoundTHnsPatch': np.zeros(self.n_frames, dtype=np.integer),
            'nEngagedTHnsPatch': np.zeros(self.n_frames, dtype=np.integer),
            'nFreeTHnsPatch': np.zeros(self.n_frames, dtype=np.integer),
            'nFreeTHnsCore': np.zeros(self.n_frames, dtype=np.integer),
            'nBridgeTHnsCore': np.zeros(self.n_frames, dtype=np.integer),
            'nDangleTHnsCore': np.zeros(self.n_frames, dtype=np.integer),
            'nCisTHnsCore': np.zeros(self.n_frames, dtype=np.integer),
            'nTransTHnsCore': np.zeros(self.n_frames, dtype=np.integer)
        }
        nmon = getattr(self.parser, 'nmon')
        if self.parser.topology == 'linear':
            self.loop_length_hist_t = np.zeros(nmon, dtype=int)
        elif self.parser.topology == 'ring':
            self.loop_length_hist_t = np.zeros((nmon//2)+1, dtype=int)
        else:
            raise ValueError(
                "The genomic distance is not defined for"
                f" '{self.parser.topology}' topology"
            )

    def _probe_frame(self, idx) -> None:
        """
        Probes a single trajectory frame and updates histogram collectors.

        Parameters
        ----------
        idx : int
            The index of the trajectory frame to probe.
        """
        self._collect_biaxial_confinement_size(idx, axis=2)
        self._collect_gyration_data(idx, unwrap=False)

        bond_dummy, cosine_dummy = correlations.bond_info(
            self._atom_groups['Mon'],
            self.parser.topology,
            )
        self.collectors['bondLengthVecMon'] += bond_dummy
        self.collectors['bondCosineCorrVecMon'] += cosine_dummy
        
        
    
    def save_artifacts(self) -> None:
        """
        Saves all collected analysis data to the specified output directory.
        """
        self.collectors['bondLengthVecMon'] = \
            self.collectors['bondLengthVecMon'] / self.n_frames
        self.collectors['bondLengthVecMon'] = \
            self.collectors['bondLengthVecMon'].reshape(self.n_bonds,)
        bonds_per_lag = np.arange(self.n_bonds, 0, -1)
        self.collectors['bondCosineCorrVecMon'] = \
           self.collectors['bondCosineCorrVecMon'] / (self.n_frames *
                                                      bonds_per_lag)
        for prop, artifact in self.collectors.items():
            filename = f"{self.save_to}{self.sim_name}-{prop}.npy"
            np.save(filename, artifact)
        print("All artifacts saved.")


class HnsCylAllProber(
    CylindricalHistogramMixIn,
    CartesianHistogram2DMixin,
    ProberBase
):
    """
    Probes simulations of the LAMMPS 'all' atom group in the *HnsCylAll*
    molecular dynamics project to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> (`Mon` for monomers, `Hns` for H-NS proteins , and `Crd`
    for crowders), the following spatial distributions are extracted:

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
    (crowders), atoms with `resid 1` represent `Mon` (monomers), and atoms with
    `type 3` represent `Hns` (H-NS proteins)
    - Coordinates are wrapped and unscaled. The center of mass of the `bug`
    group is recentered to the simulation box's center.

    Examples
    --------
    Creating an instance of `HnsCylAllProber` for a specific simulation:

    >>> prober = HnsCylAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Hns', 'Crd']

    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)

        self._parser: HnsCylInstance = HnsCyl(trajectory, lineage, 'all')
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
            self.universe.select_atoms("resid 0")  # Crowders
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # Monomers
        self._atom_groups['Hns'] = \
            self.universe.select_atoms("type 3")  # H-NS proteins

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
            'rEdge': {'bin_size': 0.1 * min(dmon, dcrowd), 
                      'lmin': 0,
                      'lmax': 0.5 * dcyl},
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'thetaEdge': {'bin_size': np.pi / 36,
                          'lmin': -np.pi,
                          'lmax': np.pi}
        }
    
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
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
            'zEdge': {'bin_size': 0.5 * min(dmon, dcrowd),
                      'lmin': -0.5 * lcyl,
                      'lmax': 0.5 * lcyl},
            'xEdge': {'bin_size': 0.1 * min(dmon, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl},
            'yEdge': {'bin_size': 0.1 * min(dmon, dcrowd),
                      'lmin': -0.5 * dcyl,
                      'lmax': 0.5 * dcyl}
        }

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
            self._update_histogram_1d('r', radial_cyl_histogram, ent, 2)
            self._update_histogram_1d('z', axial_histogram, ent, 2)
            self._update_histogram_1d('theta', azimuth_cyl_histogram, ent, 2)

            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)


class HnsCubNucleoidProber(ProberBase):


class HnsCubAllProber(
    SphericalHistogramMixin,
    CartesianHistogram2DMixin,
    ProberBase
):
    """
    Base class for probing simulations of the LAMMPS 'all' atom group in the in the *TransFociCub*, *SumRuleCubHeteroRing*, and *HnsCubAll* molecular dynamics projects to extract spatial distributions of particles.

    Physical Properties Extracted
    -----------------------------
    For each <entity> (`Mon` for monomers, `Hns` for H-NS proteins , and `Crd`
    for crowders), the following spatial distributions are extracted:

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
    Creating an instance of `HnsCubAllProber` for a specific simulation:

    >>> prober = HnsCubAllProber(
    ...     topology="topology.all.data",
    ...     trajectory="trajectory.all.lammpstrj",
    ...     lineage="whole",
    ...     save_to="./results/",
    ...     continuous=False
    ... )
    >>> prober.run()
    >>> prober.save_artifacts()
    """
    _entities = ['Mon', 'Hns', 'Crd']
    
    def __init__(
        self,
        topology: str,
        trajectory: str,
        lineage: PrimitiveLineageT,
        save_to: str,
        continuous: bool = False,
    ) -> None:
        super().__init__(topology, trajectory, lineage, save_to,
                         continuous=continuous)
        self._parser: HnsCubInstance = HnsCub(trajectory, lineage, 'all')
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
            self.universe.select_atoms("resid 0")  # Crowders
        self._atom_groups['Mon'] = \
            self.universe.select_atoms("resid 1")  # Monomers
        self._atom_groups['Hns'] = \
            self.universe.select_atoms("type 3")  # H-NS proteins

    def _prepare(self) -> None:
        self._hist1d_props: Dict[str, Dict[str, Any]] = {}
        self._hist2d_props: Dict[str, Dict[str, Any]] = {}
        bin_edges_1d = self._initialize_bin_edges_1d()
        bin_edges_2d = self._initialize_bin_edges_2d()
        for entity in self._entities:
            self._initialize_hist_collectors_1d(entity, bin_edges_1d)
            self._initialize_hist_collectors_2d(entity, bin_edges_2d)
        print("bin edges data ('xEdge', 'yEdge', 'zEdge', 'thetaEdge',"
                  " 'rEdge') saved to storage for each particle entity.")

    def _initialize_bin_edges_1d(self) -> Dict[str, Dict[str, Any]]:
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
        }
    
    def _initialize_bin_edges_2d(self) -> Dict[str, Dict[str, Any]]:
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
            self._update_histogram_1d('r', radial_histogram, ent)
            # 2D histograms
            for plane in self._hist2d_planes_dirs:
                self._update_histogram_2d(plane, ent)