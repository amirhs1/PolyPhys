"""
A submodule to calculate various 1D spatial distributions of finite (not
point-like) particles (i.e. beads) in spherical, cylindrical, cartesian, or
polar coordinate system in a box, slit, or cylindrical geometry.

Finite particle can be hard or soft. The volume of a hard particle is set by
its geometrical shape. However, the volume of a soft particle is given by the
probability distribution function (PDF) of its mass around its center of mass
(COM) or its center of geometry (COG). Such a PDF can have some connection
with the soft potential defined for the interaction of this particle with
other particles in a system.

Currently, the following distribution are implemented:
    1. The local number density of hard beads.
    2. The local volume fraction of hard beads.
"""
from glob import glob
import warnings
from abc import ABC, abstractmethod
from typing import Dict, Tuple, Optional, Callable, Literal, List
import numpy as np
import pandas as pd
from scipy import integrate
from ..manage.typer import ParserType, FreqDataT, EdgeDataT, EntityT
from ..manage.organizer import make_database, sort_filenames
from ..manage.utilizer import invalid_keyword, round_up_nearest
from ..analyze.measurer import spherical_segment, sphere_sphere_intersection


class ConsecutiveDiskBins:
    """
    Handle volume distribution for a bead in consecutive disk-like bins along a
    longitudinal axis.

    This class manages bin boundaries and volume share calculations for beads
    within consecutive disk-like bins along a specified axis (e.g., z-axis) in
    a system with periodic boundary conditions (PBC). Assuming the bead resides
    at the center of its bin, it may intersect several neighboring bins. This
    class calculates the shared volume between the bead and each intersected
    bin as the bead moves from the leftmost to the rightmost bin along the
    longitudinal axis.

    Parameters
    ----------
    edges : np.ndarray
        Array of bin edges along the longitudinal direction.
    r_bead : float
        Radius of the bead (particle) for volume calculations.

    Attributes
    ----------
    edges : np.ndarray
        Array of bin edges along the longitudinal direction.
    r_bead : float
        Radius of the bead.
    centers : np.ndarray
        Array of bin centers calculated as midpoints between edges.
    box_length : float
        Total length of the simulation box (used for PBC adjustments).
    bounds : np.ndarray
        Array marking the leftmost and rightmost bins intersected by each bead.
    volume_shares : Dict[int, Dict[int, float]]
        Nested dictionary where each key is a bin index containing a dictionary
        of volume shares for each intersected bin.

    Methods
    -------
    _find_bounds() -> None
        Determine the left and right bounds of bins intersected by each bead.
    _calculate_vol_shares() -> None
        Calculate the volume share for each bin intersected by each bead.

    Raises
    ------
    TypeError
        If `edges` is not an instance of `np.ndarray`.

    Notes
    -----
    This class accounts for periodic boundary conditions, allowing particles
    near the edges of the system to contribute volume shares to bins on both
    ends of the axis.
    """
    def __init__(self, edges: np.ndarray, r_bead: float) -> None:
        if not isinstance(edges, np.ndarray):
            raise TypeError("'edges' must be an 'np.ndarray',"
                            f" got {type(edges)}.")
        self.edges = edges
        self.r_bead = r_bead
        self.box_length = self.edges[-1] - self.edges[0]
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        n_centers = len(self.centers)
        self.bounds = np.empty([n_centers, 2])
        self.volume_shares: Dict[int, Dict[int, float]] = {}
        self._find_bounds()
        self._calculate_vol_shares()

    def _find_bounds(self) -> None:
        """
        Find left- and right-most bin edges, a.k.a. bounds, by which a bead is
        limited, adjusting bin bounds to account for periodic boundary
        conditions.
        """
        # The left- and right-most bounds by which a bead is limited:
        leftmost = self.centers - self.r_bead
        rightmost = self.centers + self.r_bead

        # Apply PBC to conserve volume
        leftmost = np.where(
            leftmost < self.edges[0],
            leftmost + self.box_length,
            leftmost
        )
        rightmost = np.where(
            rightmost > self.edges[-1],
            rightmost - self.box_length,
            rightmost
        )
        # Initiate the leftmost with the lowest possible center index:
        left_bound = np.zeros(len(leftmost), dtype=int)
        # Initiate the rightmost with the highest possible center index:
        right_bound = (len(rightmost) - 1) * np.ones(len(rightmost), dtype=int)

        for i, left_val in enumerate(leftmost):
            for edge_i in range(len(self.edges) - 1):
                if self.edges[edge_i] <= left_val < self.edges[edge_i + 1]:
                    # left_bound is the index of the smaller edges of
                    # the leftmost bin
                    left_bound[i] = edge_i
                if self.edges[edge_i] <= rightmost[i] < self.edges[edge_i + 1]:
                    # right_bound is the index of the smaller edges of
                    # the rightmost bin
                    right_bound[i] = edge_i
        self.bounds = np.column_stack((left_bound, right_bound))

    def _calculate_vol_shares(self) -> None:
        """
        Calculate the intersection volume between a bead and each intersected
        bin, adjusting for periodic boundaries when necessary.
        """
        for center_idx, bounds_minmax in enumerate(self.bounds):
            center = self.centers[center_idx]
            self.volume_shares[center_idx] = {}
            # A bead near the end of a system along a periodic direction
            # contributes to bins at both ends of that direction. If
            # bounds_minmax[0] > bounds_minmax[1], the BPC should be imposed:
            if bounds_minmax[0] > bounds_minmax[1]:
                # If the center is smaller than the half the box length, then
                # the bead's center is moved to the right, so the volume of
                # intersection can be calculated correctly:
                if center_idx <= len(self.centers) // 2:
                    center += self.box_length
                # The first element of the bounds_minmax is larger:
                for edge_idx in range(bounds_minmax[0], len(self.edges) - 1):
                    # Distances of a bin's left and right edges of from the
                    # center of the bead (in reality, the center of another
                    # bin):
                    left_dist = self.edges[edge_idx] - center
                    right_dist = self.edges[edge_idx + 1] - center
                    # the most right bound can be a spherical cap or
                    # segment; the `spherical_segment` is used to find
                    # the volume of intersection:
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(self.r_bead, left_dist, right_dist)
                # Reset center to its initial value to calculate the rest
                # of shares:
                center = self.centers[center_idx]
                # If the center is larger the half the box length, then
                # the bead's center is moved to the left, so the volume of
                # intersection can be calculated correctly.
                if center_idx > len(self.centers) // 2:
                    center -= self.box_length
                for edge_idx in range(bounds_minmax[1] + 1):
                    left_dist = self.edges[edge_idx] - center
                    right_dist = self.edges[edge_idx + 1] - center
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(self.r_bead, left_dist, right_dist)
            # When bounds_minmax[0] <= bounds_minmax[1], everything is normal
            else:
                for edge_idx in range(bounds_minmax[0], bounds_minmax[1] + 1):
                    left_dist = self.edges[edge_idx] - center
                    right_dist = self.edges[edge_idx + 1] - center
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(self.r_bead, left_dist, right_dist)


class ConcentricSphericalBins:
    """
    Handle volume distribution for a bead in concentric spherical bins along a
    longitudinal axis.

    This class manages bin boundaries and volume share calculations for beads
    within concentric spherical bins along the radial axis. Assuming the bead
    resides at the center of its bin, it may intersect several neighboring
    bins. This class calculates the shared volume between the bead and each
    intersected bin as the bead moves from the innermost to the outermost bin
    along the radial axis.

    Parameters
    ----------
    edges : np.ndarray
        Array of bin edges along the radial direction.
    r_bead : float
        Radius of the bead (particle) for volume calculations.

    Attributes
    ----------
    edges : np.ndarray
        Array of bin edges along the radial direction.
    r_bead : float
        Radius of the bead.
    centers : np.ndarray
        Array of bin centers calculated as midpoints between edges.
    bounds : np.ndarray
        Array marking the inner and outer bin boundaries for each bead.
    volume_shares : dict of dict
        Nested dictionary where each key is a bin index containing a dictionary
        of volume shares for each intersected bin.

    Methods
    -------
    _find_bounds() -> None
        Determine the inner and outer bounds of bins intersected by each bead.
    _calculate_vol_shares() -> None
        Calculate the volume share for each concentric bin intersected by each
        bead.

    Raises
    ------
    TypeError
        If `edges` is not an instance of `np.ndarray`.
    """
    def __init__(self, edges: np.ndarray, r_bead: float) -> None:
        if not isinstance(edges, np.ndarray):
            raise TypeError("'edges' must be an 'np.ndarray',"
                            f" got {type(edges)}.")
        self.edges = edges
        self.r_bead = r_bead
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        n_centers = len(self.centers)
        self.bounds = np.empty([n_centers, 2])
        self.volume_shares: Dict[int, Dict[int, float]] = {}
        self._find_bounds()
        self._calculate_vol_shares()

    def _find_bounds(self) -> None:
        """
        Find inner- and outer-most bin edges, a.k.a. bounds, by which a bead is
        limited, constraining bounds within the edges range.
        """
        innermost = self.centers - self.r_bead
        outermost = self.centers + self.r_bead

        # Constrain bounds within edges range
        innermost = \
            np.where(innermost < self.edges[0], self.edges[0], innermost)
        outermost = \
            np.where(outermost > self.edges[-1], self.edges[-1], outermost)

        # Initiate the innermost bound with the lowest possible bound
        inner_bound = np.zeros(len(innermost), dtype=int)
        # Initiate the outermost bound with the highest possible bound
        outer_bound = (len(outermost) - 1) * np.ones(len(outermost), dtype=int)

        for i, inner_val in enumerate(innermost):
            for edge_i in range(len(self.edges) - 1):
                if self.edges[edge_i] <= inner_val < self.edges[edge_i + 1]:
                    # inner_bound is the index of the larger edge of the
                    # innermost bin, since the intersection of the bead
                    # and that edge is important:
                    inner_bound[i] = edge_i + 1
                if self.edges[edge_i] <= outermost[i] < self.edges[edge_i + 1]:
                    # outer_bound is the index of the smaller edge of the
                    # outermost bin, since the intersection of the bead
                    # and that edge is important:
                    outer_bound[i] = edge_i
        self.bounds = np.column_stack((inner_bound, outer_bound))

    def _calculate_vol_shares(self) -> None:
        """
        Calculate the intersection volume between a bead and each intersected
        bin, adjusting for periodic boundaries when necessary.
        """
        for center_idx, bound_minmax in enumerate(self.bounds):
            self.volume_shares[center_idx] = {}
            # The intersection volume of the previous bin; for the innermost
            # bin, there is no previous bin, so this quantity is initiated
            # by 0:
            intersect_vol_previous = 0.0
            for edge_idx in range(bound_minmax[0], bound_minmax[1] + 2):
                # The difference between the intersection volume of a bead
                # with a bin edge with index idx and the previous bin is the
                # volume share of the bin (or center) with index idx-1:
                intersect_vol = sphere_sphere_intersection(
                    self.r_bead, self.edges[edge_idx], self.centers[center_idx]
                )
                self.volume_shares[center_idx][edge_idx - 1] = (
                    intersect_vol - intersect_vol_previous
                )
                intersect_vol_previous = intersect_vol


class DistributionBase(ABC):
    """
    Base class for managing distributions in different geometries.

    This class provides core functionality for calculating local number density
    (`rho`) and local volume fraction (`phi`) for various geometries. It
    defines shared attributes and methods that are inherited by
    geometry-specific subclasses, each of which implements geometry-specific
    behavior via the `_get_geometry_args` method.

    Parameters
    ----------
    frequencies : np.ndarray
        Array of frequencies in each bin, representing particle counts or
        measurements within each bin.
    edges : np.ndarray
        Array of bin edges, where `len(edges) = len(frequencies) + 1`.
    integrand : Callable
        Integrand function specific to the geometry and direction for volume
        integration.
    volume_shares : Dict[int, Dict[int, float]]
        Nested dictionary mapping each bin index to another dictionary that
        contains the volume shares of particles within intersected bins.
    normalize : bool, optional
        If True, normalizes `rho` and `phi` to unit area under the distribution
        (default is False).

    Attributes
    ----------
    frequencies : np.ndarray
        Frequency values for each bin, provided during initialization.
    n_bins : int
        The total number of bins, derived from `frequencies`.
    edges : np.ndarray
        Array of bin edges along the chosen direction.
    bin_size : np.ndarray
        Array of bin widths, calculated from consecutive differences in
        `edges`.
    integrand : Callable
        The integrand function used for volume integration within bins.
    args : tuple
        Tuple of arguments specific to the geometry and direction, used by the
        integrand.
    volume_shares : Dict[int, Dict[int, float]]
        Nested dictionary where each key is a bead index containing a
        dictionary of volume shares for each intersected bin.
    normalize : bool
        Flag indicating whether the distributions should be normalized.
    rho : np.ndarray
        Calculated local number density values for each bin, obtained by
        dividing `frequencies` by the volume of each bin.
    phi : np.ndarray
        Calculated local volume fraction values for each bin, representing the
        fraction of bin volume occupied by particles.

    Methods
    -------
    _set_args() -> Tuple[float, ...]
        Sets up and returns geometry-specific arguments for the integrand by
        calling the subclass-defined `_get_geometry_args`.
    _get_geometry_args() -> Tuple[float, ...]
        Abstract method; must be implemented by subclasses to return geometry-
        specific arguments for the integrand.
    _number_density() -> None
        Calculates the local number density (`rho`) for each bin by dividing
        the frequency values by bin volumes.
    _volume_fraction() -> None
        Calculates the local volume fraction (`phi`) for each bin by summing
        the product of `rho` and volume shares across intersected bins.
    _normalize_distributions() -> None
        Normalizes the `rho` and `phi` arrays if `normalize` is set to True.
    save_number_density(filename: str) -> None
        Saves the `rho` array (local number density) to a `.npy` file.
    save_volume_fraction(filename: str) -> None
        Saves the `phi` array (local volume fraction) to a `.npy` file.

    Raises
    ------
    TypeError
        If `frequencies` or `edges` is not an instance of `np.ndarray`.
    """

    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        integrand: Callable,
        volume_shares: Dict[int, Dict[int, float]],
        normalize: bool = False
    ) -> None:
        if not isinstance(frequencies, np.ndarray):
            raise TypeError("'frequencies' must be an 'np.ndarray',"
                            f" got {type(frequencies)}.")
        if not isinstance(edges, np.ndarray):
            raise TypeError("'edges' must be an 'np.ndarray',"
                            f" got {type(edges)}.")

        self.frequencies = frequencies
        self.n_bins = len(self.frequencies)
        self.edges = edges
        self.bin_size = self.edges[1:] - self.edges[:-1]
        self.integrand = integrand
        self.volume_shares = volume_shares
        self.normalize = normalize
        self.rho = np.zeros(self.n_bins)
        self.phi = np.zeros(self.n_bins)

        # Set up geometry-specific arguments
        self.args = self._set_args()

        # Calculate number density and volume fraction
        self._number_density()
        self._volume_fraction()

        if self.normalize:
            self._normalize_distributions()

    def _set_args(self) -> Tuple[float, ...]:
        """
        Set up geometry-specific arguments for the integrand.

        This method calls the subclass-specific `_get_geometry_args` method to
        retrieve arguments that are unique to the geometry and direction.

        Returns
        -------
        tuple
            Geometry-specific integration arguments.
        """
        return self._get_geometry_args()

    @abstractmethod
    def _get_geometry_args(self) -> Tuple[float, ...]:
        """
        Abstract method to provide geometry-specific arguments.

        This method must be implemented by subclasses to supply arguments
        specific to the geometry being represented.

        Returns
        -------
        tuple
            Geometry-specific integration arguments.
        """

    def _number_density(self) -> None:
        """
        Calculate the local number density (`rho`) for each bin.

        Integrates the `integrand` over each bin to obtain bin volumes, and
        then divides `frequencies` by the volumes to get `rho`.
        """
        bin_vols = np.array([
            integrate.quad(
                self.integrand,
                self.edges[idx],
                self.edges[idx + 1],
                args=self.args
            )[0] for idx in range(self.n_bins)
        ])
        self.rho = self.frequencies / bin_vols

    def _volume_fraction(self) -> None:
        """
        Calculate the local volume fraction (`phi`) for each bin.

        Uses `volume_shares` to accumulate the fraction of each intersected bin
        that the particles occupy, based on `rho` values.
        """
        for c_idx in range(self.n_bins):
            for idx, vol in self.volume_shares[c_idx].items():
                self.phi[c_idx] += self.rho[idx] * vol

    def _normalize_distributions(self) -> None:
        """
        Normalize distributions.
        """
        self.rho /= self.rho.sum()
        self.phi /= self.phi.sum()

    def save_number_density(self, filename: str) -> None:
        """
        Save the local number density (`rho`) to a `.npy` binary file.

        Parameters
        ----------
        filename : str
            Base filename or path to save the density data. A `.npy` extension
            will be appended to the filename if it does not already have one.

        Notes
        -----
        The `.npy` extension is handled internally by `np.save`.
        """
        np.save(filename, self.rho)

    def save_volume_fraction(self, filename: str) -> None:
        """
        Save the local volume fraction density (`phi`) to a `.npy` binary file.

        Parameters
        ----------
        filename : str
            Base filename or path to save the volume fraction data. A `.npy`
            extension will be appended to the filename if it does not already
            have one.

        Notes
        -----
        The `.npy` extension is handled internally by `np.save`.
        """
        np.save(filename, self.phi)


class DistributionCylinder(DistributionBase):
    """
    Calculate local distributions (local number density and volume fraction)
    within a cylindrical geometry.

    This class handles calculations for `rho` (local number density) and `phi`
    (local volume fraction) in a cylindrical system:

        - Along the longitudinal axis, it utilizes consecutive disk-like bins
        and handles periodic boundary conditions as necessary.

    Parameters
    ----------
    frequencies : np.ndarray
        Array of frequencies in each bin, representing particle counts or
        measurements within each bin.
    edges : np.ndarray
        Array of bin edges along the longitudinal axis, where
        `len(edges) = len(frequencies) + 1`.
    rbead : float
        Radius of the particle species being measured.
    lcyl : float
        Length of the cylinder along the longitudinal axis.
    dcyl : float
        Diameter of the cylinder's cross-section.
    direction : str
        Direction of interest for the distribution calculation, which must be
        one of ['r', 'phi', 'z'].
    normalize : bool, optional
        If True, normalizes `rho` and `phi` to a unit area under the
        distribution (default is False).

    Attributes
    ----------
    lcyl : float
        Length of the cylinder, used as an argument in the integration.
    dcyl : float
        Diameter of the cylinder's cross-section, used for integration.
    direction : {'r', 'phi', 'z'}
        Direction along which the distribution is computed within the
        cylindrical geometry
    integrand : Callable
        Integrand function used for volume integration, chosen based on
        `direction`.
    args : Tuple[float, ...]
        Arguments specific to a cylindrical direction, set by
        `_get_geometry_args`.
    consecutive_bins : ConsecutiveDiskBins
        Instance of `ConsecutiveDiskBins` used to manage volume shares of bins.

    Methods
    -------
    _get_geometry_args() -> Tuple[float, ...]
        Provides geometry-specific arguments for cylindrical integration,
        based on the specified `direction`.

    Raises
    ------
    ValueError
        If `direction` has invalid values.
    NotImplementedError
        If `direction` is 'r' or 'phi'.
    """
    _directions = ['r', 'phi', 'z']
    _integrands: Dict[str, Callable] = {
        'r': lambda r, lcyl: 2 * np.pi * lcyl * r,
        'phi': lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2,
        'z': lambda z, dcyl: 0.25 * np.pi * dcyl**2
    }

    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        rbead: float,
        lcyl: float,
        dcyl: float,
        direction: Literal['r', 'phi', 'z'],
        normalize: bool = False
    ) -> None:
        self.lcyl = lcyl
        self.dcyl = dcyl
        invalid_keyword(direction, self._directions)
        if direction in ['r', 'phi']:
            raise NotImplementedError(
                "This class has not yet implemented for 'r' and 'phi'"
                "directions"
            )
        self.direction = direction
        self.integrand = self._integrands[direction]
        self.consecutive_bins = ConsecutiveDiskBins(edges, rbead)

        super().__init__(
            frequencies,
            edges,
            self.integrand,
            self.consecutive_bins.volume_shares,
            normalize=normalize
        )

    def _get_geometry_args(self) -> Tuple[float, ...]:
        """Provide arguments specific to the cylindrical geometry."""
        direction_args = {
            'r': (self.lcyl,),
            'phi': (self.lcyl, self.dcyl),
            'z': (self.lcyl,)
        }
        return direction_args[self.direction]

    #np.save(filaname+self.direction+, self.phi)
    #    np.save(filaname, self.rho)


class DistributionSphere(DistributionBase):
    """
    Calculate local distributions (local number density and volume fraction)
    within a spherical geometry.

    This class handles calculations for `rho` (local number density) and `phi`
    (local volume fraction) in a spherical system, using concentric spherical
    bins. The spherical geometry can accommodate different directions (radial,
    azimuthal, and polar).

    Parameters
    ----------
    frequencies : np.ndarray
        Array of frequencies in each bin, representing particle counts or
        measurements within each bin.
    edges : np.ndarray
        Array of bin edges along the radial direction, where
        `len(edges) = len(frequencies) + 1`.
    rbead : float
        Radius of the particle species being measured.
    lcube : float
        Length of the cube's side, used as an effective diameter for spherical
        symmetry.
    direction : str
        Direction of interest for the distribution calculation, which must be
        one of ['r', 'theta', 'phi'].
    normalize : bool, optional
        If True, normalizes `rho` and `phi` to a unit area under the
        distribution (default is False).

    Attributes
    ----------
    lcube : float
        Effective radius or side length of the cube in which the spherical
        distribution is computed.
    direction : {'r', 'theta', 'phi'}
        Direction along which the distribution is computed within the
        spherical geometry.
    integrand : Callable
        Integrand function used for volume integration, chosen based on
        `direction`.
    args : Tuple[float, ...]
        Arguments specific to the spherical geometry and direction, set by
        `_get_geometry_args`.
    concentric_bins : ConcentricSphericalBins
        Instance of `ConcentricSphericalBins` used to manage volume shares of
        bins.

    Methods
    -------
    _get_geometry_args() -> Tuple[float, ...]
        Provides geometry-specific arguments for spherical integration,
        based on the specified `direction`.

    Raises
    ------
    ValueError
        If `direction` has invalid values.
    NotImplementedError
        If `direction` is 'theta' or 'phi'.
    """
    _directions = ['r', 'theta', 'phi']
    _integrands: Dict[str, Callable] = {
        'r': lambda r, const: 4 * const * np.pi * r**2,
        # constant is redundant and merely defined to consistently use
        # of 'args' parameter of "scipi.integral.quad" among integrands.
        'theta': lambda theta, lcube: np.pi * lcube**3 * np.sin(theta) / 12,
        'phi': lambda lcube: lcube**3 / 12
    }

    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        rbead: float,
        lcube: float,
        direction: Literal['r', 'theta', 'phi'],
        normalize: bool = False
    ) -> None:
        self.lcube = lcube
        invalid_keyword(direction, self._directions)
        if direction in ['theta', 'phi']:
            raise NotImplementedError(
                "This class has not yet implemented for 'r' and 'phi'"
                "directions"
            )
        self.direction = direction
        self.integrand = self._integrands[direction]
        self.concentric_bins = ConcentricSphericalBins(edges, rbead)

        super().__init__(
            frequencies,
            edges,
            self.integrand,
            self.concentric_bins.volume_shares,
            normalize=normalize
        )

    def _get_geometry_args(self) -> Tuple[float, ...]:
        """Provide arguments specific to the cylindrical geometry."""
        direction_args: Dict[str, Tuple[float, ...]] = {
            'r': (1,),  # set const to 1, see above for more discussion
            'theta': (self.lcube,),
            'phi': (self.lcube,)
        }
        return direction_args[self.direction]


def distributions_generator(
    freqs: FreqDataT,
    bin_edges: EdgeDataT,
    group: str,
    species: str,
    geometry: str,
    topology: str,
    direction: str,
    parser: ParserType,
    save_to: Optional[str] = None,
    normalized: bool = False
) -> Tuple[Dict, Dict]:
    """
    Generates the local number densities (rho) and volume fractions (phi).

    Parameters
    ----------
    freqs: dict
        A dictionary in which the keys are the filenames and values are 1D
        frequencies.
    bin_edges:
        A dictionary in which the keys are the filenames and values are 1D
        bin edges. A 1D array of the bin edge where `len(edge)=len(freq)`.
    group: str in {'bug', 'all'}
        The type of the particle group.
    species: str in {'Mon', 'Crd', 'Foci'}
        The species of particles.
    geometry: str in {'cylindrical', 'slit', 'cubic'}
        The shape of the simulation box.
    topology: str in {'ring', 'linear'}
        The topology of the polymer.
    direction: str in {'azimuthal', 'polar', 'radial', 'longitudinal'}
        The direction of interest in the `geometry` of interest.
    parser: ParserT
        A class from 'PolyPhys.manage.parser' module that parses filenames
        or filepaths to infer information about a file.
    save_to: str
        path to which the output saved.
    normalized: bool
        whether normalized the distributions or not.

    Return
    ------
    densities: dict
        A dictionary in which the keys are 'whole' names and the values are the
        local number densities.
    vol_fractions: dict
        A dictionary in which the keys are 'whole' names and the values are the
        local volume fractions.
    """
    densities = {}
    vol_fractions = {}
    parser_name = parser.__name__
    radius_attrs = {
        'SumRuleCyl': {
            'Mon': 'dmon',
            'Crd': 'dcrowd'
        },
        'TransFociCyl': {
            'Mon': 'dmon_small',
            'Crd': 'dcrowd',
            'Foci': 'dmon_large'
        },
        'TransFociCub': {
            'Mon': 'dmon_small',
            'Crd': 'dcrowd',
            'Foci': 'dmon_large'
        },
        'SumRuleCubHeteroRing': {
            'Mon': 'dmon_small',
            'Crd': 'dcrowd',
            'Foci': 'dmon_large'
        },
        'SumRuleCubHeteroLinear': {
            'Mon': 'dmon_small',
            'Crd': 'dcrowd',
            'Foci': 'dmon_large'
        },
        'HnsCub': {
            'Mon': 'dmon',
            'Hns': 'dhns',
            'Crd': 'dcrowd'
        },
        'HnsCyl': {
            'Mon': 'dmon',
            'Hns': 'dhns',
            'Crd': 'dcrowd'
        }
    }
    lbox_attrs = {
        'cubic': {
            'TransFociCub': 'lcube',
            'SumRuleCubHeteroLinear': 'lcube',
            'SumRuleCubHeteroRing': 'lcube',
            'HnsCub': 'lcube'
        },
        'cylindrical': {
            'SumRuleCyl': 'lcyl',
            'TransFociCyl': 'lcyl',
            'HnsCyl': 'lcyl'
        }
    }
    dbox_attrs = {
        'cubic': {
            'TransFociCub': 'N/A',
            'SumRuleCubHeteroLinear': 'N/A',
            'SumRuleCubHeteroRing': 'N/A',
            'HnsCub': 'N/A'
        },
        'cylindrical': {
            'SumRuleCyl': 'dcyl',
            'TransFociCyl': 'dcyl',
            'HnsCyl': 'dcyl'
        }
    }
    geometries = ['cubic', 'cylindrical']
    if geometry not in geometries:
        raise NotImplementedError("The number density and volume fraction"
                                  f"are currently available in '{geometries}")
    for whole, freq in freqs.items():
        whole_info = parser(
            whole,
            'whole',
            geometry,
            group,
            topology,
            ispath=False,
        )
        if geometry == 'cubic':
            distributions = DistributionSphere(
                freq,
                bin_edges[whole],
                radius_attrs[parser_name][species],
                lbox_attrs[geometry][parser_name],
                dbox_attrs[geometry][parser_name],
                direction,
                normalized=normalized
            )
        else:
            densities[whole] = distributions.rho
            vol_fractions[whole] = distributions.phi

        if save_to is not None:
            filename = f"{save_to}{whole}-{group}-{direction}"
            distributions.save_number_density(f"{filename}Rho{species}")
            distributions.save_volume_fraction(f"{filename}Phi{species}")

    return densities, vol_fractions


def entropic_energy(
    histo_collections: np.ndarray,
    beta: Optional[float] = 1.0
) -> np.ndarray:
    """
    this is the effective free energy potential as the function of the
    end-to-end distance. Here, we use the end-to-end distance as the
    action parameters and write a Langevin equation for that.

    Parameters
    ----------
    histo_collections: ArrayLike
        radial distribution function or the probability distribution
        function of *the magnitude of the end-to-end vector *
        (or *the end-to-end distance*) between R and R+dR.
    beta: float
        k_B*T is the inverse of the thermal energy.

    Returns
    -------
    free_energy: np.ndarray
        The dimensionless free energy potential of the end-to-end distance

        -beta * Ln(P(r)*4pi*r^2) = -beta * Ln(histo_collections) where
        beta=k_BT is set 1.0.
    """
    histo_collections = histo_collections / np.sum(histo_collections)
    free_energy = -1.0 * (
        np.log(histo_collections) - np.log(np.sum(histo_collections))) / beta
    return free_energy


def looping_prop(
    histo_collections: np.ndarray,
    bin_edges: np.ndarray,
    r_max: float,
    r_min: float
) -> float:
    """
    Calculate the looping probability, defined as the probability that two
    chain ends are within a specific distance range.

    The looping entropy, :math:`P_L`, is calculated as the integral from
    `r_min` to `r_max` over the distribution :math:`P(R) \\cdot 4 \\pi R^2 dR`.
    Since :math:`P(R) \\cdot 4 \\pi R^2` is approximated by
    `histo_collections[i]` for each bin between `bin_edges[i]` and
    `bin_edges[i+1]`, the looping probability is given by:

    .. math::
        P_L \\approx \\sum_{i} P(\\text{bin_center}_i) \\cdot 4\\pi
        \\cdot \\text{bin_center}_i^2 \\cdot (\\text{bin_edges}_{i+1} -
        \\text{bin_edges}_i).

    Parameters
    ----------
    histo_collections : np.ndarray
        Array containing histogram values of the Flory radius or end-to-end
        distance over the simulation.
    bin_edges : np.ndarray
        Array of bin edges corresponding to `histo_collections`.
    r_max : float
        Maximum center-to-center distance between two chain ends for looping
        to occur.
    r_min : float
        Minimum center-to-center distance between two chain ends where crowders
        cannot fit, allowing non-zero depletion forces.

    Returns
    -------
    float
        The probability of looping.

    Notes
    -----
    - The function normalizes the histogram to ensure total probability equals
    one.
    - It is assumed that bins are eqaul in width, so they do not change the
    probability after normalization.

    Examples
    --------
    >>> histo_collections = np.array([0.1, 0.2, 0.3, 0.4])
    >>> bin_edges = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    >>> r_max = 2.5
    >>> r_min = 0.5
    >>> looping_p(histo_collections, bin_edges, r_max, r_min)
    0.6
    """
    # Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Normalize the histogram by bin width and total count
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    normalized_hists = histo_collections * bin_widths
    normalized_hists /= np.sum(normalized_hists)

    # Vectorized computation of looping probability
    within_range = bin_centers <= (r_max + r_min)
    prob = np.sum(normalized_hists[within_range])

    return prob


def normalize_z(
    prop: str,
    ens_avg: pd.DataFrame,
    trans_symmetry: bool = True
) -> pd.DataFrame:
    """
    Normalizes the ensemble-average local distribution `ens_avg` of the
    specified property `prop` along the z-axis in cylindrical geometry.
    Normalization is performed by dividing by the maximum value of `ens_avg`.
    If `trans_symmetry` is `True`, it averages over the absolute values of bin
    centers along the z-axis.


    Parameters
    ----------
    prop : str
        Name of the physical property to normalize.
    ens_avg : pd.DataFrame
        DataFrame containing the ensemble-average local distribution with a
        column named `<prop>-scale`.
    trans_symmetry : bool, default True
        If `True`, the system has translational symmetry along the z axis, so
        averaging is performed over absolute values of bin centers along this
        axis. See *Notes* below.

    Returns
    -------
    pd.DataFrame
        Normalized ensemble-average local distribution. If `trans_symmetry` is
        `True`, returns the averaged distribution over the absolute values of
        bin centers; otherwise, returns the original DataFrame with
        normalization applied.

    Examples
    --------
    >>> data = pd.DataFrame({"bin_center": [-1, 0, 1],
                             "density-scale": [0.5, 1.0, 0.5]})
    >>> normalize_z("density", data)
       bin_center  density-scale  density-normalizer  density-norm
    0          -1            0.5                 1.0           0.5
    1           0            1.0                 1.0           1.0
    2           1            0.5                 1.0           0.5

    Notes
    -----
    - This function assumes the z values range from negative to positive. If
      this condition is not met, consider using the `normalize_r` function.

    - **Bin Center Behavior After Averaging**:
      If `trans_symmetry` is `True`, the function averages the data based on
      the absolute values of bin centers along the z-axis. When the bin
      centers are symmetrical (equal positive and negative values around
      (:math:`z = 0`), the absolute values of bin centers remain unchanged in
      the averaged output.

      However, if bin centers are asymmetrical (unequal positive and negative
      values), the output bin centers represent the mean of the absolute
      values of each pair of positive and negative bin centers. For odd-length
      DataFrames, the middle bin center is retained as is, as it represents
      (:math:`z = 0`) or the closest value to zero.

    - **Indexing Strategy**:
      The function uses indexing to handle bin center pairs:
        - For an odd-length `ens_avg`, the middle bin center value represents
          (:math:`z = 0`) (or the closest to zero) and remains as is in the
          final output while the average between each positive and negative bin
          center pair is computed.
        - For an even-length `ens_avg`, the average between each positive
          and negative bin center pair is computed, with no (:math:`z = 0`)
          center etained.

    """
    ens_avg_max = ens_avg[f'{prop}-scale'].max()
    ens_avg[f'{prop}-normalizer'] = ens_avg_max

    if ens_avg_max != 0:
        ens_avg[f'{prop}-norm'] = ens_avg[f'{prop}-scale'] / ens_avg_max
    else:
        warnings.warn(
            "All values are zero; normalized values set to zero.", UserWarning)
        ens_avg[f'{prop}-norm'] = 0

    if not trans_symmetry:
        return ens_avg

    df_len = len(ens_avg)
    mid_point = df_len // 2
    odd_length = bool(df_len % 2 == 1)

    # Separate data into positive and negative z regions:
    ens_neg = ens_avg.iloc[:mid_point].copy()
    ens_neg.index = -1 * ens_neg.index
    ens_neg.sort_index(inplace=True)
    ens_neg.reset_index(inplace=True, drop=True)
    ens_neg['bin_center'] *= -1

    ens_pos = ens_avg.iloc[mid_point + odd_length:].copy()
    ens_pos.reset_index(inplace=True, drop=True)

    # Averaging over |z|>0:
    ens_sum = 0.5 * (ens_pos + ens_neg)

    # Handling z=0 or z close to 0:
    if odd_length:
        ens_sum.set_index('bin_center', inplace=True)
        ens_sum.loc[0, :] = ens_avg.loc[mid_point, :]

    ens_sum.sort_index(inplace=True)
    ens_sum.reset_index(drop=True)

    return ens_sum


def normalize_r(
    prop: str,
    ens_avg: pd.DataFrame,
    method: Literal['first', 'max', 'range_mean', 'range_max'] = 'first',
    index_range: Optional[slice] = None
) -> pd.DataFrame:
    """
    Normalizes the ensemble-average local distribution `ens_avg` of `prop`
    along the radial (r) direction by a specified `method`.

    Parameters
    ----------
    prop : str
        Name of the physical property to normalize.
    ens_avg : pd.DataFrame
        DataFrame containing the ensemble-average local distribution, with a
        column named `<prop>-scale`.
    method : {'first', 'max', 'range_mean', 'range_max'}, default 'first'
        Normalization method:

        - 'first': Normalizes by the first value in `ens_avg`.
        - 'max': Normalizes by the maximum value in `ens_avg`.
        - 'range_mean': Normalizes by the mean of values within the range
          specified by `index_range`.
        - 'range_max': Normalizes by the maximum value within the range
          specified by `index_range`.

    index_range : Optional[slice], default None
        A slice specifying the index range over which the mean or maximum of
        `prop` in `ens_avg` is computed. Required if `method` is 'range_mean'
        or 'range_max'.

    Returns
    -------
    pd.DataFrame
        DataFrame with the normalized ensemble-average local distribution,
        where the new column `<prop>-norm` holds the normalized values.

    Raises
    ------
    ValueError
        If `index_range` is not specified when `method` is 'range_mean' or
        'range_max'.

    Examples
    --------
    >>> data = pd.DataFrame({"bin_center": [0, 1, 2],
                             "density-scale": [1.0, 0.8, 0.6]})
    >>> normalize_r("density", data, method="max")
       bin_center  density-scale  density-normalizer  density-norm
    0           0            1.0                 1.0           1.0
    1           1            0.8                 1.0           0.8
    2           2            0.6                 1.0           0.6

    Notes
    -----
    - Ensure `index_range` is specified if using the 'range_mean' or
    'range_max' methods. For instance, `index_range=slice(start, end)` allows
    indexing over a defined section of `ens_avg`.
    """
    # Determine the normalizing value
    normalizer: float = 0
    if method == 'first':
        normalizer = ens_avg[prop + '-scale'].iloc[0]
    elif method == 'max':
        normalizer = ens_avg[prop + '-scale'].max()
    elif method == 'range_mean' and index_range is not None:
        normalizer = ens_avg.loc[index_range, prop + '-scale'].mean()
    elif method == 'range_max' and index_range is not None:
        normalizer = ens_avg.loc[index_range, prop + '-scale'].max()
    else:
        invalid_keyword(method, ['first', 'max', 'range_mean', 'range_max'])

    # Apply normalization
    ens_avg[f'{prop}-normalizer'] = normalizer
    if normalizer != 0:
        ens_avg[f'{prop}-norm'] = ens_avg[f'{prop}-scale'] / normalizer
    else:
        warnings.warn(
            "All values are zero; normalized values set to zero.",
            UserWarning
        )
        ens_avg[f'{prop}-norm'] = 0

    return ens_avg


def space_sum_rule(
    input_database: str,
    property_: str,
    parser: ParserType,
    hierarchy: str,
    physical_attrs: List[str],
    species: Literal['Mon', 'Crd', 'Foci', 'Dna'],
    size_attr: str,
    group: Literal["bug", "nucleoid", "all"],
    geometry: Literal["cylindrical", "slit", "cubic"],
    topology: Literal["ring", "linear", "branched"],
    direction: Literal["r", "z"],
    divisor: float = 0.025,
    round_to: int = 3,
    is_save: bool = False
) -> pd.DataFrame:
    """
    Processes ensemble-averaged local distributions of a given `property_`
    in an `input_database`, applies normalization and scaling, and aggregates
    attributes.

    Each "ensemble-averaged" DataFrame contains columns like:
    - `<long_ensemble>-<property_>[-<measure>]-<stat>`
    - `bin_center`

    where <measure> is a physical measurement such as the auto correlation
    function (AFC) done on the physical 'property_'. [...] means this keyword
    in the column name can be optional. The <stat> keyword is either 'mean',
    'ver', or 'sem'. If the `bin_center` presents as a column in a
    'ensemble_averaged' dataframe, then it is inferred; otherwise, it should
    be passed to the function. See `bin_center` kw argument below.

    "Scaling" changes the data range, and "normalizing" changes its
    distribution. This function creates a normalized, scaled distribution for
    a given particle species along the specified `direction`.

    Parameters
    ----------
    input_database : str
        Path to the directory containing the ensemble-averaged local
        distributions.
    property_ : str
        Physical property of interest, e.g., 'density' or 'volume fraction'.
    parser : ParserT
        Parser instance for file path and attribute information.
    hierarchy : str
        Pattern prefix for filenames to select files, e.g., "N*".
    physical_attrs : List[str]
        List of physical attributes to add to the output DataFrame.
    species : {'Mon', 'Crd', 'Foci', 'Dna'}
        Particle species within the distribution:
        - 'Mon': Monomers
        - 'Crd': Crowders
        - 'Foci': Large monomers
        - 'Dna': DNA particles

    size_attr : str
        Size attribute (diameter) of the species for scaling.
    group : {'bug', 'nucleoid', 'all'}
        Particle group type.
    geometry : {'cylindrical', 'slit', 'cubic'}
        Simulation box geometry.
    topology : {'ring', 'linear', 'branched'}
        Polymer topology.
    direction : {'r', 'z'}
        Direction along which the operation is performed.
    divisor : float, default 0.025
        Step for rounding `phi_c_bulk` attribute values.
    round_to : int, default 3
        Decimal precision for `phi_c_bulk` attribute rounding.
    is_save : bool, default False
        If True, saves the output DataFrame to a CSV file.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of normalized, scaled distributions with
        added attributes.

    Raises
    ------
    NotImplementedError
        If `property_` is neither 'Phi' nor 'Rho'.

    Notes
    -----
    This function applies scaling based on particle size and species:
    - `Phi` (volume fraction) is scaled by the particle size.
    - `Rho` (density) is scaled by the particle size squared.

    Examples
    --------
    >>> space_sum_rule(input_database="data/",
                       property_="density",
                       parser=my_parser,
                       hierarchy="N*",
                       physical_attrs=["temp", "phi_c_bulk"],
                       species="Mon",
                       size_attr="diameter",
                       group="all",
                       geometry="cylindrical",
                       topology="linear",
                       direction="r")
    """
    invalid_keyword(species, ['Mon', 'Crd', 'Foci', 'Dna'])
    invalid_keyword(group, ["bug", "nucleoid", "all"])
    invalid_keyword(geometry, ["cylindrical", "slit", "cubic"])
    invalid_keyword(topology, ["ring", "linear", "branched"])
    invalid_keyword(direction, ["r", "z "])
    # Normalizer based on direction
    normalizer: Dict[str, Callable] = {
        'r': normalize_r,
        'z': normalize_z
    }

    # File extension pattern
    property_ext = f"-{group}-{direction}{property_}{species}-ensAvg.csv"
    prop = f"{direction}{property_}{species}"  # Full property name

    # Collect ensemble-averaged CSVs
    ens_avg_csvs = sort_filenames(
        glob(f"{input_database}{hierarchy}{property_ext}"),
        fmts=[property_ext])

    property_dfs = []

    for ens_avg_csv in ens_avg_csvs:
        ens_avg = pd.read_csv(ens_avg_csv[0], header=0)
        property_info: ParserType = parser(
            ens_avg_csv[0],
            'ensemble_long',
            geometry,
            group,
            topology
        )

        # Apply scaling based on property type
        scaler = getattr(property_info, size_attr)
        if property_ == 'Phi':
            ens_avg[f"{prop}-scaler"] = scaler
            ens_avg[f"{prop}-scale"] = ens_avg[f"{prop}-mean"] / scaler
        elif property_ == 'Rho':
            ens_avg[f"{prop}-scaler"] = scaler ** 2
            ens_avg[f"{prop}-scale"] = ens_avg[f"{prop}-mean"] * scaler ** 2
        else:
            raise NotImplementedError(
                "Sum rule scaler is only defined for 'Phi' (volume fraction)"
                " or 'Rho' (density) properties."
            )

        # Normalizing the area under curve
        ens_avg[f"{prop}-scale-normalized_curve"] = \
            ens_avg[f"{prop}-scale"] / ens_avg[f"{prop}-scale"].sum()

        # Normalize data using specified normalizer function
        if direction == 'r' and geometry == 'cubic':
            ens_avg = normalizer[direction](prop, ens_avg, method='max')
        else:
            ens_avg = normalizer[direction](prop, ens_avg)

        # Add specified physical attributes
        for attr_name in physical_attrs:
            ens_avg[attr_name] = getattr(property_info, attr_name)

        # Add additional attributes
        ens_avg[f"{prop}-sumrule_constant"] = scaler
        ens_avg['bin_center-norm'] = \
            ens_avg['bin_center'] / ens_avg['bin_center'].max()
        ens_avg['phi_c_bulk_round'] = \
            ens_avg['phi_c_bulk'].apply(
                lambda x: round_up_nearest(x, divisor, round_to)
            )

        # Additional processing for cylindrical geometry
        if geometry == 'cylindrical':
            ens_avg['temp'] = (
                (ens_avg['dcyl'] % ens_avg['dcrowd']) /
                (ens_avg['dcrowd'])
            )
            ens_avg['bin_center-dcrowd-recentered'] = (
                ens_avg['bin_center-dcrowd'] - ens_avg['temp']
            )
            ens_avg['bin_center-recentered-norm'] = (
                ens_avg['bin_center'] - (ens_avg['dcyl'] % ens_avg['dcrowd'])
            )
            ens_avg['bin_center-recentered-norm'] = (
                ens_avg['bin_center-recentered-norm'] /
                ens_avg['bin_center-recentered-norm'].max()
            )
            ens_avg.drop(columns=['temp'], inplace=True)

        # Store processed DataFrame
        property_dfs.append(ens_avg)

    # Concatenate all processed data
    property_db = pd.concat(property_dfs, axis=0)
    property_db.reset_index(drop=True, inplace=True)

    # Save DataFrame if requested
    if is_save:
        save_to_space = make_database(
            input_database, "analysis", stage="space", group=group)
        space = save_to_space.split("/")[-2].split("-")[0]
        filepath = save_to_space + \
            f"{space}-{group}-{property_}-{species}" + \
            "-normalizedRescaled-space.csv"
        property_db.to_csv(filepath, index=False)
    return property_db
