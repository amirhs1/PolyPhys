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
from typing import Dict, Tuple, Optional, Callable
import numpy as np
from numpy.typing import ArrayLike
import scipy.integrate as integrate
from ..manage.typer import (ParserT, FreqDataT, EdgeDataT)
from ..manage.organizer import make_database


def spherical_segment(r: float, a: float, b: float) -> float:
    """
    Compute the volume of a spherical segment defined by two parallel planes.

    Parameters
    ----------
    r : float
        Radius of the sphere. Must be positive.
    a : float
        Distance of the first plane from the center of the sphere, along the
        axis of symmetry.
    b : float
        Distance of the second plane from the center of the sphere, along the
        axis of symmetry.

    Returns
    -------
    vol : float
        Volume of the spherical segment.

    Raises
    ------
    ValueError
        If `r` is not a positive number.

    Notes
    -----
    - `a` and `b` can be positive or negative values, as long as they fall
    within the range `[-r, r]`.
    - The function will adjust `a` and `b` to `-r` or `r` if they exceed these
    bounds.
    - If `a = r` or `b = r`, the spherical segment becomes a spherical cap.
    - If both `a` and `b` lie outside the range `[-r, r]` and share the same
    sign, the volume is zero.

    References
    ----------
    .. [1] Weisstein, Eric W. "Spherical Segment." From MathWorld--A Wolfram
    Web Resource.
       https://mathworld.wolfram.com/SphericalSegment.html

    Examples
    --------
    >>> spherical_segment(3, 1, 2)
    37.69911184307752
    >>> spherical_segment(3, -3, 3)
    113.09733552923255
    """

    if r <= 0:
        raise ValueError(f"The radius 'r' must be positive. Got {r}.")

    # Ensure the bounds are within [-r, r]
    lower = max(min(a, b), -r)
    upper = min(max(a, b), r)

    # If both bounds are outside [-r, r] with the same sign, the volume is zero
    if lower * upper >= r**2:
        return 0.0
    # Calculate the volume of the spherical segment
    vol = np.pi * (r**2 * (upper - lower) - (upper**3 - lower**3) / 3)
    return vol


def sphere_sphere_intersection(r1: float, r2: float, d: float) -> float:
    """
    Compute the volume of intersection between two spheres.

    The sphere with radius `r1` is separated from the sphere with radius `r2`
    by a distance `d` along the x-axis. Thus, the vector form of the distance
    between their centers is `(d, 0, 0)` in Cartesian coordinates.

    Parameters
    ----------
    r1 : float
        Radius of the first sphere.
    r2 : float
        Radius of the second sphere.
    d : float
        Distance between the centers of the two spheres along the x-axis.

    Returns
    -------
    vol : float
        Volume of the intersection between the two spheres.

    References
    ----------
    .. [1] Weisstein, Eric W. "Sphere-Sphere Intersection."
       From MathWorld--A Wolfram Web Resource.
       https://mathworld.wolfram.com/Sphere-SphereIntersection.html

    Examples
    --------
    >>> sphere_sphere_intersection(3, 4, 2)
    75.39822368615503
    >>> sphere_sphere_intersection(3, 4, 10)
    0.0
    """

    r_max = max(r1, r2)
    r_min = min(r1, r2)

    # Volume is zero if one sphere has a radius of zero or no intersection
    # occurs:
    if r1 == 0.0 or r2 == 0.0:
        return 0.0
    if d >= r_min + r_max:
        return 0.0
    if d <= r_max - r_min:
        # The smaller sphere is entirely contained within the larger one
        return 4 * np.pi * r_min**3 / 3

    # Calculate the intersection volume for overlapping spheres
    vol = (np.pi / (12 * d)) * (r_max + r_min - d)**2 * (
        d**2 + 2 * d * (r_max + r_min)
        - 3 * (r_min**2 + r_max**2)
        + 6 * r_min * r_max
    )

    return vol


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
    _find_bounds()
        Determine the left and right bounds of bins intersected by each bead.
    _calculate_vol_shares()
        Calculate the volume share for each bin intersected by each bead.

    Notes
    -----
    This class accounts for periodic boundary conditions, allowing particles
    near the edges of the system to contribute volume shares to bins on both
    ends of the axis.
    """
    def __init__(self, edges: np.ndarray, r_bead: float) -> None:
        self.edges = edges
        self.r_bead = r_bead
        self.box_length = self.edges[-1] - self.edges[0]
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        n_centers = len(self.centers)
        self.bounds = np.empty([n_centers, 2])
        self.volume_shares: Dict[int, Dict[int, float]] = {}
        self.find_bounds()
        self._calculate_vol_shares()

    def find_bounds(self) -> None:
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
    _find_bounds()
        Determine the inner and outer bounds of bins intersected by each bead.
    _calculate_vol_shares()
        Calculate the volume share for each concentric bin intersected by each
        bead.

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


class DistributionBase:
    """
    Base class for managing distributions in different geometries.

    Parameters
    ----------
    frequencies : np.ndarray
        Array of frequencies in each bin.
    normalized : bool, optional
        Whether to normalize distributions (default is False).

    Attributes
    ----------
    frequencies : np.ndarray
        Frequency values for bins.
    normalized : bool
        Flag indicating if normalization is applied.
    rho : np.ndarray
        Calculated local number densities.
    phi : np.ndarray
        Calculated local volume fractions.

    Raises
    ------
    ValueError
        If `frequencies` is not an ndarray.
    """
    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        integrand: Callable,
        args: Tuple[float, ...],
        volume_shares: Dict[int, Dict[int, float]],
        normalize: bool = False
    ) -> None:
        # Validate frequency and edge inputs
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
        self.args = args
        self.volume_shares = volume_shares
        self.normalize = normalize
        self.rho = np.zeros(self.n_bins)
        self.phi = np.zeros(self.n_bins)
        self._number_density()
        self._volume_fraction()
        if self.normalize is True:
            self._normalize_distributions()

    def _number_density(self) -> None:
        """
        Calculate local number density along the specified direction.

        Notes
        -----
        The local number density is averaged over time steps collected in each
        simulation. For example, in a simulation with 7*10^7 total time steps
        sampled every 5000 steps, there are 14001 measurements.
        """
        bin_vols = np.array([
            integrate.quad(
                self.integrand,
                self.edges[idx],
                self.edges[idx+1],
                args=self.args
            )[0] for idx in range(self.n_bins)
        ])
        self.rho = self.frequencies / bin_vols

    def _volume_fraction(self) -> None:
        """Calculate volume fractions for cylindrical distributions."""
        volume_shares = self.volume_shares
        for c_idx in range(self.n_bins):
            for idx, vol in volume_shares[c_idx].items():
                self.phi[c_idx] += self.rho[idx] * vol

    def _normalize_distributions(self) -> None:
        """Normalize distributions."""
        self.rho /= self.rho.sum()
        self.phi /= self.phi.sum()


class DistributionCylinder(DistributionBase):
    """
    Calculate local distributions (local number density and volume fraction)
    across a cylinder.

    Parameters
    ----------
    frequencies : np.ndarray
        Array of frequencies in each bin.
    edges : np.ndarray
        Array of bin edges along the longitudinal axis.
    radius_attr : float
        Radius of the particle species.
    normalized : bool, optional
        Whether to normalize distributions (default is False).

    Attributes
    ----------
    volume_shares : dict of dict
        Nested dictionary where each key is a bead index containing a
        dictionary of volume shares for each intersected bin.
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
        direction: str,
        normalize: bool = False
    ) -> None:
        self.lcyl = lcyl
        self.dcyl = dcyl
        if direction not in self._directions:
            raise ValueError(
                f"Invalid direction '{direction}'. "
                f"Choose from {self._directions}."
            )
        self.direction = direction
        self.args = self._set_args()
        self.consective_bins = ConsecutiveDiskBins(edges, rbead)
        super().__init__(
            frequencies,
            edges,
            self.integrand,
            self.args,
            self.consective_bins.volume_shares,
            normalize=normalize
        )

    def _set_args(self) -> Tuple[float, ...]:
        """
        Set integration arguments for different directions.

        Notes
        -----
        The radius of the simulation box should be used as the radial argument
        `r`. A pre-factor of 0.5 is applied to `lbox` and `dbox` in the
        integrands for radial directions to maintain consistency.
        """
        direction_args: Dict[str, Tuple[float, ...]] = {
            'r': (self.lcyl,),
            'theta': (self.lcyl, self.dcyl),
            'z': (self.lcyl,)
        }
        return direction_args[self.direction]


class DistributionSpherical(DistributionBase):
    """
    Calculte local distributions (number density and volume fraction) across a
    sphere.
    """
    _directions = ['r', 'theta', 'phi']
    _integrands: Dict[str, Callable] = {
        'r': lambda r, const: 4 * const * np.pi * r**2,
        # constant is redundant and merely defined to consistently use
        # of 'args' parameter of "scipi.integral.quad" among integrands.
        'theta': lambda theta, dcyl: np.pi * dcyl**3 * np.sin(theta) / 12,
        'phi': lambda theta, dcyl: dcyl**3 / 12
    }

    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        rbead: float,
        lcube: float,
        direction: str,
        normalize: bool = False
    ) -> None:
        self.lcube = lcube
        if direction not in self._directions:
            raise ValueError(
                f"Invalid direction '{direction}'. "
                f"Choose from {self._directions}."
            )
        self.direction = direction
        self.args = self._set_args()
        self.concentric_bins = ConcentricSphericalBins(edges, rbead)
        super().__init__(
            frequencies,
            edges,
            self.integrand,
            self.args,
            self.concentric_bins.volume_shares,
            normalize=normalize
        )

    def _set_args(self) -> Tuple[float, ...]:
        """
        Set integration arguments for different directions.

        Notes
        -----
        The radius of the simulation box should be used as the radial argument
        `r`. A pre-factor of 0.5 is applied to `lbox` and `dbox` in the
        integrands for radial directions to maintain consistency.
        """
        direction_args: Dict[str, Tuple[float, ...]] = {
            'r': (1,),  # set const to 1, see above for more discussion
            'theta': (self.lcube,),
            'phi': (self.lcube,)
        }
        return direction_args[self.direction]

class SpatialDistribution:
    """
    Calculate local number density and volume fraction of spherical species.

    This class computes the local number density and local volume fraction
    of species within a 1D, 2D, or 3D space, with periodic boundary conditions
    (PBC) along specific axes:

    - 1D: Cylinder with PBC along z-axis (longitudinal).
    - 2D: Slit with PBC along x and y axes (radial).
    - 3D: Cube with PBC along x, y, and z axes (spherical).

    The radius of each species is needed for computing the local volume
    fraction.

    To-do List
    ----------
    - Vectorize the for-loops.

    Parameters
    ----------
    freqs : np.ndarray
        A 1D array of frequencies in bins.
    edges : np.ndarray
        A 1D array of bin edges, where `len(edges) = len(freqs) + 1`.
    hist_info : ParserT
        Parser containing information about the system.
    radius_attr : str
        Attribute name in `hist_info` with the species radius.
    lbox : str
        Length attribute for the simulation box.
    dbox : str
        Diameter attribute for cylindrical or slit geometries.
    geometry : str
        Shape of the simulation box, one of ['cylindrical', 'slit', 'cubic'].
    direction : str
        Direction along which to compute distributions.
    normalized : bool, optional
        Normalize distributions to unit area under the curve (default is
        False).

    Attributes
    ----------
    frequencies : np.ndarray
        Frequency values for bins.
    edges : np.ndarray
        Edges of the bins.
    rho : pd.DataFrame
        Local number densities with the same index and column layout as 
        `frequencies`.
    phi : pd.DataFrame
        Local volume fractions, structured similarly to `frequencies`.
    bounds : np.ndarray
        Indices marking lower and upper bin bounds per bead center.
    volume_shares : dict of dict
        Nested dictionary mapping bin indices to the volume share of each bin
        for particles centered within that bin.

    Notes
    -----
    For calculating volume fractions, this class assumes that beads' centers
    reside at bin centers. This assumption reduces computational cost, with
    error decreasing as the number of bins increases for a fixed system volume.

    Bin schemes:
        - 'concentric': Concentric shellsin spherical or cylindrical geometry.
        - 'consecutive': Consecutive slits (e.g., disks along the z direction
        in the cylindrical coordinate system) in Cartesian or cylindrical
        geometry.
        - 'periodic': Bins with periodicity (e.g., the polar direction the
        cylindrical coordinate system) in cylindrical or spherical
        geometry.
        - 'neighboring': Azimuthal sectors in spherical geometry.

    The periodic boundary condition (PBC) is imposed as needed in each
    direction and binning scheme.

    It is assumed all the bins have the same size.

    When `normalized` is set to `True`, the sum of `rho` is not equal to the
    bulk number density when r approaches infinity, i.e. natom/vol_system. This
    arises from the way we discretize the local number density. The sum of
    `phi` is also not equal to the bulk volume fraction when r approaches
    infinity, i.e, (natoms*vol_atom)/vol_system. This arises from the way we
    discretize the local number density.

    """
    _directions = {
        'cubic': ['r', 'theta', 'phi'],
        'slit': ['r', 'phi', 'z'],
        'cylindrical': ['r', 'phi', 'z']
    }
    _integrands: Dict[str, Dict[str, Callable]] = {
        'cubic': {
            'r': lambda r, const: 4 * const * np.pi * r**2,
            # constant is redundant and merely defined to consistently use
            # of 'args' parameter of "scipi.integral.quad" among integrands.
            'theta': lambda theta, dcyl: np.pi * dcyl**3 * np.sin(theta) / 12,
            'phi': lambda theta, dcyl: dcyl**3 / 12
        },
        'slit': {
            'r': lambda r, lcyl: 2 * np.pi * lcyl * r,
            'phi': lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2,
            'z': lambda z, dcyl: 0.25 * np.pi * dcyl**2
        },
        'cylindrical': {
            'r': lambda r, lcyl: 2 * np.pi * lcyl * r,
            'phi': lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2,
            'z': lambda z, dcyl: 0.25 * np.pi * dcyl**2
        }
    }
    _fullnames = {
        'cubic': {
            'r': 'radial',
            'theta': 'azimuthal',
            'phi': 'polar'
        },
        'slit': {
            'r': 'radial',
            'phi': 'polar',
            'z': 'longitudinal'
        },
        'cylindrical': {
            'r': 'radial',
            'phi': 'polar',
            'z': 'longitudinal'
        }
    }

    def __init__(
        self,
        frequencies: np.ndarray,
        edges: np.ndarray,
        hist_info: ParserT,
        radius_attr: str,
        lbox: str,
        dbox: str,
        geometry: str,
        direction: str,
        normalized: bool = False,
    ) -> None:
        """
        Initialize SpatialDistribution with frequency data and geometry settings.

        Parameters
        ----------
        frequencies : np.ndarray
            Array of frequencies in each bin.
        edges : np.ndarray
            Array of bin edges, where `len(edges) = len(frequencies) + 1`.
        hist_info : ParserT
            Parser containing system information.
        radius_attr : str
            Attribute name in `hist_info` with species radius.
        lbox : str
            Attribute name in `hist_info` for box length.
        dbox : str
            Attribute name in `hist_info` for box diameter in cylindrical/slit 
            geometries.
        geometry : str
            Shape of the simulation box, one of ['cylindrical', 'slit', 'cubic'].
        direction : str
            Direction along which distributions are computed.
        normalized : bool, optional
            Normalize distributions to unit area (default is False).

        Raises
        ------
        ValueError
            If `frequencies` or `edges` is not an ndarray.
            If `geometry` or `direction` is invalid.
        """
        # Validate frequency and edge inputs
        if not isinstance(frequencies, np.ndarray):
            raise ValueError("frequencies must be an np.ndarray,"
                             f" got {type(frequencies)}.")
        if not isinstance(edges, np.ndarray):
            raise ValueError("edges must be an np.ndarray,"
                             f" got {type(edges)}.")

        self.frequencies = frequencies
        self.edges = edges
        self.box_length = self.edges[-1] - self.edges[0]
        self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        self.bin_size = self.edges[1] - self.edges[0]

        self.hist_info = hist_info
        self.r_bead = 0.5 * getattr(self.hist_info, radius_attr)
        self.lbox = getattr(self.hist_info, lbox)
        self.dbox = getattr(self.hist_info, dbox) \
            if dbox != 'N/A' else getattr(self.hist_info, lbox)

        # Validate geometry and direction
        valid_geometries = list(self._directions.keys())
        if geometry not in valid_geometries:
            raise ValueError(f"Invalid geometry '{geometry}'. Choose from"
                             f" {valid_geometries}.")
        self.geometry = geometry

        if direction not in self._directions[self.geometry]:
            valid_directions = self._directions[self.geometry]
            raise ValueError(
                f"Invalid direction '{direction}' for '{geometry}'. "
                f"Choose from {valid_directions}."
            )
        self.direction = direction
        self.normalized = normalized

        # Initialize calculations
        self._set_args()
        self._number_density()

        if self.direction == 'z':
            self.bin_handler = \
                ConsecutiveDiskLikeBins(self.edges, self.r_bead)
        elif self.direction == 'r':
            self.bin_handler = \
                ConcentricSphericalShellBins(self.edges, self.r_bead)
        else:
            raise ValueError(f"Unsupported direction '{self.direction}'.")
        self._volume_fraction()

        if self.normalized:
            # Normalize rho and phi
            self.rho /= self.rho.sum()
            self.phi /= self.phi.sum()

    def _set_args(self) -> None:
        """
        Set integration arguments for each geometry and direction.

        Notes
        -----
        The radius of the simulation box should be used as the radial argument
        `r`. A pre-factor of 0.5 is applied to `lbox` and `dbox` in the
        integrands for radial directions to maintain consistency.
        """
        self._args: Dict[str, Dict[str, Tuple[float, ...]]] = {
            'cubic': {
                'r': (1,),
                'theta': (self.lbox,),
                'phi': (self.lbox,)
            },
            'cylindrical': {
                'r': (self.lbox,),
                'theta': (self.lbox, self.dbox),
                'z': (self.lbox,)
            }
        }

    def _number_density(self) -> None:
        """
        Calculate local number density along the specified direction.

        Notes
        -----
        The local number density is averaged over time steps collected in each
        simulation. For example, in a simulation with 7*10^7 total time steps
        sampled every 5000 steps, there are 14001 measurements.
        """
        integrand: Callable = self._integrands[self.geometry][self.direction]
        arguments: Tuple[float, ...] = \
            self._args[self.geometry][self.direction]
        bin_vols = np.array([
            integrate.quad(
                integrand, self.edges[idx], self.edges[idx] + self.bin_size, 
                args=arguments
            )[0] for idx in range(len(self.edges) - 1)
        ])
        self.rho = self.frequencies / bin_vols

    def _volume_fraction(self) -> None:
        """
        Calculate the local volume fraction along the direction of interest.

        Notes
        -----
        Assumes all particles have identical shapes. The local volume fraction 
        is normalized, giving the integral over the volume of interest.
        """
        n_centers = len(self.rho)
        volume_shares = self.bin_handler.volume_shares
        self.phi = np.zeros(n_centers)

        for c_idx in range(n_centers):
            for idx, vol in volume_shares[c_idx].items():
                self.phi[c_idx] += self.rho[idx] * vol


def distributions_generator(
    freqs: FreqDataT,
    bin_edges: EdgeDataT,
    group: str,
    species: str,
    geometry: str,
    topology: str,
    direction: str,
    parser: Callable,
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
        distributions = SpatialDistribution(
            freq,
            bin_edges[whole],
            whole_info,
            radius_attrs[parser_name][species],
            lbox_attrs[geometry][parser_name],
            dbox_attrs[geometry][parser_name],
            geometry,
            direction,
            normalized=normalized
        )
        densities[whole] = distributions.rho
        vol_fractions[whole] = distributions.phi
        if save_to is not None:
            filename = whole + '-' + group + '-' + direction
            np.save(
                save_to + filename + 'Rho' + species + '.npy',
                distributions.rho
            )
            np.save(
                save_to + filename + 'Phi' + species + '.npy',
                distributions.phi
            )
    return densities, vol_fractions


def entropic_energy(
    histo_collections: ArrayLike,
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


def looping_p(
    histo_collections: np.ndarray,
    bin_edges: np.ndarray,
    r_max: float,
    r_min: float
) -> float:
    """
    Looping entropy P_L is defined as the integral from r_min to r_max over
    P(R)*4pi*R^2*dR; since P(R)*4pi*R^2 is equivalent to hitso_collection[i]
    for the bin between bin_edges[i] and bin_edges[i+1] we have
    P_L = integral from r_min to r_max over P(R)*4pi*R^2*dR ~ sum from r_min
    to r_max over P(bin_center_i)*4pi*bin_centers_i^2*(bin_edges_i+1 -
    bin_edges_i)=sum from r_min to r_max over hitso_collection_i*(bin_edges_i+1
    - bin_edges_i)

    Since the sizes of the bins are equal ((bin_edges_i+1 - bin_edges_i)
    =constant), the effect of (bin_edges_i+1 - bin_edges_i) is cancelled out
    upon normalization.

    Parameters
    ----------
    histo_collections: a numpy array of the histograms of the Flory radius
    or the end-to-end distance during the whole run (the total number of
    data points is equal to the total number of time steps.)
    bin_edges: the edges of the bins used in histograms.
    r_max: the minimum center-to-center distance the two end-of-chain
    monomers can have; Since they have the same size, this is usually equal
    to their diameter (or size).
    r_min: this is the size of the crowders that can reside between two
    monomers. When the center-to-center distance between two monomers is
    less than r_max+r_min, no crowders can be between two monomers and the
    depletion force is non-zero.

    Return
    ------
    The probability of looping.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    histo_collections = np.multiply(
        histo_collections, (bin_edges[1:] - bin_edges[:-1])
        )
    histo_collections = histo_collections / np.sum(histo_collections)
    looped_probability = 0
    for i in range(len(bin_centers)):
        if bin_centers[i] <= (r_max + r_min):
            looped_probability += histo_collections[i]
    return looped_probability


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
    parser: ParserT,
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
        property_info: ParserT = parser(
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
                round_up_nearest, args=[divisor, round_to]
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
