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


def spherical_segment(
    r: float,
    a: float,
    b: float
) -> float:
    """
    computes the volume of a spherical segment.

    Parameters
    ----------
    r : float
        Radius of the sphere.
    a : float
        Distance of the first base from the center of the sphere.
    b : float
        Distance of the second base from the center of the sphere.

    Return
    ------
    vol: float
        Volume of the spherical segment.

    Notes
    -----
    1. a and b can be negative or positive, depending on the position
    with respect to the center of sphere in the range [-r,r].
    2. When a=r or b=r, we have a spherical cap.

    References
    ----------
    https://mathworld.wolfram.com/SphericalSegment.html
    """
    if r <= 0:
        raise ValueError(
            "The radius "
            f"'{r}'"
            " is not a positive number."
        )
    # The axis, on which left, center, and right reside, points to right.
    lower = min(a, b)  # The lower bound of the volume integral
    upper = max(a, b)  # The upper bound of the volume integral
    if np.abs(upper) >= r:  # upper bound cannot larger than r.
        upper = np.sign(upper) * r
    if np.abs(lower) >= r:  # lower limit cannot less than -1*r.
        lower = np.sign(lower) * r
    vol = np.pi * (r**2 * (upper-lower) - (upper**3 - lower**3) / 3)
    # if lower and upper bounds are both outside of [-r,r] and have the
    # same signs, then vol=0:
    if lower*upper >= r**2:
        vol = 0.0
    return vol


def sphere_sphere_intersection(
        r1: float,
        r2: float,
        d: float
) -> float:
    """
    computes the volume of intersection of two spheres. The sphere with
    radius `r1` has distance `d`from the other one with radius `r2`
    along axis x, so the vector form of distance is indeed (d,0,0)
    in cartesian coordinate system.

    Parameters
    ----------
    r1: float
        The radius of the one sphere.
    r2: float
        The radius of the other sphere.
    d: float
        The distance of the two spheres from each other along axis x.

    Returns
    -------
    vol: float
        The intersection volume.

    References
    ----------
    https://mathworld.wolfram.com/Sphere-SphereIntersection.html
    """
    # By define r_max and r_min, we handle the situations in which d = 0:
    r_max = max(r1, r2)
    r_min = min(r1, r2)
    # when one of the sphere has 0 radius, vol=0:
    if r1 == 0.0 or r2 == 0.0:
        vol = 0.0
    # Spheres are either tangential to each other or do not intersect.
    elif d >= r_min + r_max:
        vol = 0.0
    # Spheres intersect:
    # Small sphere resides completely in the large one.
    elif d <= r_max - r_min:
        vol = 4 * np.pi * r_min**3 / 3
    # Other scenarios:
    else:
        vol = (np.pi / (12*d)) * (r_max + r_min - d)**2 * (
            d**2
            + 2 * d * (r_max + r_min)
            - 3 * (r_min**2 + r_max**2)
            + 6 * r_min * r_max
        )
    return vol


class SpatialDistribution(object):
    freqs: np.ndarray
    edges: np.ndarray
    hist_info: ParserT
    radius_attr: str
    lbox: str
    dbox: str
    geometry: str
    direction: str
    normalized: bool = False
    """
    Computes the local number density of any species and the local volume
    fraction of spherical species in 1D (cylinder with the periodic boundary
    condition (PBC) along the longitudinal axis (z axis)), 2D (slit with the
    PBC along the x and y axes (radial direction)), and 3D (cube with the PBC
    along x, y, and z axes (spherical)).

    The radius of the species is needed for compute local volume fraction.

    To-do List
    ----------
    Vectorize the for-loops.

    Parameters
    ----------
    freqs: np.ndarray
        A 1D  array of frequencies in bins.
    edges: np.ndarray
        A 1D array of the bin edge where `len(edge)=len(freq)`.
    hist_info: ParserT
        A ParserT object that contains information about the system to which
        `freq` belongs.
    radius_attr: str
        The name of an attribute of `hist_info` that contains the radius of the
        species to which `freq` belongs.
    radius_attr: str
        The name of an attribute of `hist_info` that contains the radius of the
        species to which `freq` belongs.
    radius_attr: str
        The name of an attribute of `hist_info` that contains the radius of the
        species to which `freq` belongs.
    geometry: str t
        The shape of the simulation box
    direction: str
        The direction of interest in the `geometry` of interest.
    normalized: bool
        Whether normalizing the distributions to have a unit value total under
        the distribution curve or not.

    Attributes
    ----------
    self.frequencies: numpy.ndarray
        Frequencies in different bins
    self.edges: numpy.ndarray
        The edges of the bins.
    self.radius_attr: float
        The radius of particle species for which distributions are computed.
    self.geometry: str
        The shape of the simulation box.
    self.direction: str
        The direction along which distribution are computed.
    self.lbox: float
        The side of a cubic box, the radius or side length of a rectangular or
        cylindrical slit, or the longitudinal length of a cylinder.
    self.dbox
        The diameter of a cylindrical slit or a cylinder.
    self.rho (pandas dataframe):
        the dataframe of the local number
        densities which has the same index, number of columns and names of
        columns as self.frequencies. However, the column values are the local
        number densities.
    self.phi (pandas dataframe):
        the dataframe of the local volume
        fractions which has the same index, number of columns and names of
        columns as self.frequencies. However, the column values are the local
        volume fractions.
    self.bounds: numpy.ndarray
        Indices of the lower and upper boundaries of a bead that resides at
        each of centers in `self.center`.
    self.volume_shares: A two-fold nested dict
        A nested dictionary in which keys are indices of `self.centers`
        (the positions of a bead's center) and the values are themselves
        dictionaries. In each of this dictionary, each key is one of the
        indices of `self.centers` within the bounds (inclusively, i.e.
        [A,B] includes A and B as well given by `self.bounds` and the
        value of that key is the share of the bin of that center from the
        volume of the bead. The sum of all the values in one inner
        dictionary is equal to the volume of the bead.
    self._args: dict
        The arguments of each of the nine _integrands.

    Notes
    -----
    To calculate the distribution of local volume fraction, we have to
    consider how the volume of a bead (particle) intersect with the bins
    defined in a given direction in a given coordinate system. This
    problem has two steps:
        1. Finding the bins (more precisely, the edges) by which the boundary
        (or volume) of a particle intersects. In practice, the lower and
        upper limits of these edges are found.
        2. Calculating the share of each bin from the volume of a bead.

    While COGs or COMs of beads can be anywhere in the system, it is assumed
    that the COGs or COMs of all the beads reside in a bin are the center of
    that bin. While this assumption produces some error in the local volume
    distribution, it significantly reduces the computational cost. Moreover,
    the error decrease by increase the number of bins at a fixed volume of the
    system. By this assumption, all the beads in a bin are equivalent and it
    is sufficient to do the two above steps only for one bead in a given
    direction in a given geometry.

    To find the upper and lower bounds on edges, the center of a bead with
    radius `self.r_bead` is placed at one of `self.centers`. Using the
    above assumption, for a given `self.r_bead`, a pair of lower and upper
    bounds (`self.bounds`) can be found for each center. It is assumed
    that a lower or upper bound is the leftmost, rightmost, innermost, or
    outermost edges smaller than the boundary of bead in a given direction.

    There are four bin schemes:

        'concentric':
            bins, edges, or centers are non-negative strictly increasing
            arrays, creating concentric circles, spherical shells, or
            cylindrical shells along the radial direction in the polar,
            spherical, or cylindrical coordinate system respectively.

        'consecutive':
            bins, edges, or centers are strictly increasing arrays,
            creating consecutive slits with rectangular or circular cross
            sections along the different directions in the cartesian
            coordinate system or the longitudinal direction in the
            cylindrical coordinate system respectively.

        'periodic':
            bins, edges, or centers are strictly increasing but bounded
            by the period of the direction of interest, thus creating
            periodic circular (or disk) sectors or spherical wedges along
            the polar direction with period [0,2*Pi] in the polar and
            cylindrical coordinate systems or the spherical coordinate
            system respectively, or generating concentric.

        'neighboring:
            bins, edges, or centers are strictly increasing but bounded
            by the period of the direction of interest, thus generating
            sectors created by two neighboring azimuth (conic) planes
            along the azimuth direction in the spherical coordinate system.

    For each of this bin type, an algorithm is defined below to find the
    lower and upper bounds.

    Finding bond edges along a given direction is a 1D problem. Moreover,
    the following assumptions about positions of beads (particles), bin
    centers and edges are used in this class:
    1. If a bead is in a bin, its position is assumed to be the center
        of that bin, not its actual position.
    2. len(self.edges) = len(self.centers) + 1
    3. For a bin, we have: left_edge <= center < right_edge.
    4. There is only one center between two consecutive edges.
    5. The edge and center arrays are strictly increasing functions of
        their indices.
    6. The bins are equally spaced, son bin_size is the same between any
        two consecutive edges.

    To conserve the total volume of the beads, the periodic boundary condition
    (PBS) is imposed if required (For instance, the z direction in the
    cylindrical geometry). The PBC is impose in both step: finding the bounds
    and calculating the volume shares.
    """
    _directions = {
        'cubic': ['r', 'theta', 'phi'],
        'slit': ['r', 'phi', 'z'],
        'cylindrical': ['r', 'phi', 'z']
    }
    _integrands = {
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
        if isinstance(frequencies, np.ndarray):
            self.frequencies = frequencies
        else:
            raise ValueError(
                f"'{frequencies}'"
                " is not a numpy.ndarray.")
        if isinstance(edges, np.ndarray):
            self.edges = edges
            self.box_length = self.edges[-1] - self.edges[0]
            self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
            # assuming bin.edges are equally-spaced
            self.bin_size = self.edges[1] - self.edges[0]
        else:
            raise ValueError(
                f"'{edges}'"
                " is not a numpy.ndarray.")
        self.hist_info = hist_info
        self.r_bead = 0.5 * getattr(self.hist_info, radius_attr)
        self.lbox = getattr(self.hist_info, lbox)
        if dbox != 'N/A':
            self.dbox = getattr(self.hist_info, dbox)
        else:  # This is meaningless, fix it!
            self.dbox = getattr(self.hist_info, lbox)
        if geometry in ['cylindrical', 'slit']:
            self.dbox = getattr(self.hist_info, dbox)
        if geometry in self._directions.keys():
            self.geometry = geometry
        else:
            geomteries_string = (
                "'" + "', '".join(self._directions.keys()) + "'"
                )
            raise ValueError(
                f"'{geometry}' is not a valid geometry for the simulation box."
                f" Please select one of {geomteries_string} geometries.")
        if direction in self._directions[self.geometry]:
            self.direction = direction
        else:
            directions_string = "'" + "', '".join(
                self._directions[self.geometry]) + "'"
            raise ValueError(
                f"'{direction}' "
                "is not a valid direction for "
                f"'{self.geometry}' geometry. Please select one of "
                f"{directions_string} directions.")
        self.normalized = normalized
        self._set_args()
        self._number_density()
        self._vol_shares_type()
        self._volume_fraction()
        if self.normalized is True:
            # the sum of rho is not equal to the bulk number density when r
            # approaches infinity, i.e. natom/vol_system. This arises from
            # the way we discretize the local number density.
            self.rho = self.rho / self.rho.sum()
            # time averaging: the sum of histograms = natoms * nframes.
            # normalization: the sum of the number density is now 1.
            # the sum of phi is not equal to the bulk volume fraction when r
            # approaches infinity, i.e, (natoms*vol_atom)/vol_system. This
            # arises from the way we discretize the local number density.
            self.phi = self.phi / self.phi.sum()

    def _set_args(self) -> None:
        """
        sets the arguments for the integrands along different directions
        in different geometries.

        Note
        ----
        In a given geometry, The radius of the simulation box should be used as
        the radial argument 'r'. Here, instead of halving 'lcube' or 'dcyl', we
        introduced a 0.5 pre-factor in the definition of the integrand.
        """
        self._args = {
            'cubic': {
                'r': (1, ),
                # For concentric spherical shells, the constant, i.e. 1,
                # is redundant and merely defined to make the use of args
                # parameter of scipi.integral.quad function consistent among
                # integrands.
                'theta': (self.lbox, ),
                # In a box cubic or free space, the radius of the space
                # is half of the length of simulation box.
                'phi': (self.lbox, ),
                # In a box cubic or free space, the radius of the space
                # is half of the length of simulation box.
            },
            'cylindrical': {
                'r': (self.lbox, ),
                'theta': (self.lbox, self.dbox, ),
                'z': (self.lbox, )
            }
        }

    def _number_density(self) -> None:
        """
        Calculates the local number density along the given direction in
        the given geometry.

        Notes
        -----
        The number density in each simulation is an average over the number
        densities collected every X time steps, so there are N=L/X
        measurements of the local number density in each simulation where L
        is total number of time steps in the simulation; for instance, in
        the cylindrical sum rule project X=5000 and the total number of
        time steps is 7*10^7, so N=14001.
        """
        integrand = self._integrands[self.geometry][self.direction]
        arguments = self._args[self.geometry][self.direction]
        bin_vols = np.array([integrate.quad(
            integrand,
            self.edges[idx],
            self.edges[idx] + self.bin_size,
            args=arguments)[0] for idx in range(len(self.edges[:-1]))
        ])
        # histogram[col_name] and bin_vols have the same size, so:
        self.rho = self.frequencies / bin_vols

    def _consecutive_bounds(self):
        """
        Finds the indices of smallest and largest bin edges by which a
        sphere, that resides at the center of that bin, intersects in a
        'consecutive' bin scheme.

        Attributes
        ----------
        self.bounds: np.ndarray
            Indices of the leftmost and rightmost boundaries of a bead
            that resides at each of centers in `self.center`. For each
            center, the following four variables are used to
            generate `self.bounds`:

            'leftmost' is he value of the center of the bin in which the
            leftmost boundary of a bead is located.

            'left_bound' is he index of the 'smaller' edge of the bin in
            which the leftmost boundary of a bead is located.

            'rightmost' is the value of the center of the bin in which the
            rightmost boundary of a bead is located.

            'right_bound' is the index of the 'smaller' edge of the bin in
            which the rightmost boundary of a bead is located.
        """
        # The left- and right-most bounds by which a bead is limited:
        leftmost = self.centers - self.r_bead
        rightmost = self.centers + self.r_bead
        # Ensuring PBC to conserve total volume of all the beads:
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
        for idx, leftmost_value in enumerate(leftmost):
            for edge_idx in range(len(self.edges) - 1):
                if (leftmost_value >= self.edges[edge_idx]) and \
                        (leftmost_value < self.edges[edge_idx + 1]):
                    # left_bound is the index of the smaller edges of
                    # the leftmost bin
                    left_bound[idx] = edge_idx
                if (rightmost[idx] >= self.edges[edge_idx]) and \
                        (rightmost[idx] < self.edges[edge_idx + 1]):
                    # right_bound is the index of the smaller edges of
                    # the rightmost bin
                    right_bound[idx] = edge_idx
        self.bounds = np.column_stack((left_bound, right_bound))

    def _consecutive_vol_shares(self):
        """
        Computes the portion of the volume of a sphere in each of the
        consecutive disk-like bins by which the bead intersects along the
        longitudinal direction with the PBC in a cylindrical geometry.

        Attributes
        ----------
        self.volume_shares: A two-fold nested dict
            A nested dictionary in which keys are indices of `self.centers`
            (the positions of a bead's center) and the values are themselves
            dictionaries. In each of this dictionary, each key is one of the
            indices of `self.centers` within the bounds (inclusively, i.e.
            [A,B] includes A and B as well given by `self.bounds` and the
            value of that key is the share of the bin of that center from the
            volume of the bead. The sum of all the values in one inner
            dictionary is equal to the volume of the bead.
        """
        self.volume_shares = {}
        for center_idx, bounds_minmax in enumerate(self.bounds):
            # A bead's center is the center of the bin in which it resides:
            center = self.centers[center_idx]
            self.volume_shares[center_idx] = {}
            # A bead near the end of a system along a periodic direction
            # contributes to bins at both ends of that direction. If
            # bounds_minmax[0] > bounds_minmax[1], the BPC should be imposed:
            if bounds_minmax[0] > bounds_minmax[1]:
                # If the center is smaller than the half the box length, then
                # the bead's center is moved to the right, so the volume of
                # intersection can be calculated correctly:
                if center_idx <= len(self.centers)//2:
                    center = self.centers[center_idx] + self.box_length
                # The first element of the bounds_minmax is larger:
                for edge_idx in range(bounds_minmax[0], len(self.edges) - 1):
                    # Distances of a bin's left and right edges of from the
                    # center of the bead (in reality, the center of another
                    # bin):
                    left_distance = self.edges[edge_idx] - center
                    right_distance = self.edges[edge_idx + 1] - center
                    # the most right bound can be a spherical cap or
                    # segment; the `spherical_segment` is used to find
                    # the volume of intersection:
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(
                            self.r_bead,
                            left_distance,
                            right_distance
                        )
                # Reset center to its initial value to calculate the rest
                # of shares:
                center = self.centers[center_idx]
                # If the center is larger the half the box length, then
                # the bead's center is moved to the left, so the volume of
                # intersection can be calculated correctly.
                if center_idx > len(self.centers)//2:
                    center = self.centers[center_idx] - self.box_length
                for edge_idx in range(bounds_minmax[1] + 1):
                    left_distance = self.edges[edge_idx] - center
                    right_distance = self.edges[edge_idx+1] - center
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(
                            self.r_bead,
                            left_distance,
                            right_distance
                        )
            # When bounds_minmax[0] <= bounds_minmax[1], everything is normalL
            else:
                for edge_idx in range(bounds_minmax[0], bounds_minmax[1] + 1):
                    left_distance = self.edges[edge_idx] - center
                    right_distance = self.edges[edge_idx+1] - center
                    self.volume_shares[center_idx][edge_idx] = \
                        spherical_segment(
                            self.r_bead,
                            left_distance,
                            right_distance
                        )

    def _concentric_bounds(self):
        """
        Finds the indices of smallest and largest bin edges by which a
        bead, that resides at the center of that bin, intersects in a
        'concentric' bin.

        Instead of the leftmost/rightmost pairs (which are more appropriate
        for the longitudinal direction), the innermost/outermost is used.
        """
        # The left- and right-most bounds by which a bead is limited:
        innermost = self.centers - self.r_bead
        outermost = self.centers + self.r_bead
        # Ensuring the bounds in the interval[edges[0], edges[-1]]:
        innermost = np.where(
            innermost < self.edges[0],
            self.edges[0],
            innermost
        )
        outermost = np.where(
            outermost > self.edges[-1],
            self.edges[-1],
            outermost
        )
        # Initiate the innermost bound with the lowest possible bound
        inner_bound = np.zeros(len(innermost), dtype=int)
        # Initiate the outermost bound with the highest possible bound
        outer_bound = (len(outermost) - 1) * np.ones(len(outermost), dtype=int)
        for idx, innermost_value in enumerate(innermost):
            for edge_idx in range(len(self.edges) - 1):
                if (innermost_value >= self.edges[edge_idx]) and \
                        (innermost_value < self.edges[edge_idx + 1]):
                    # inner_bound is the index of the larger edge of the
                    # innermost bin, since the intersection of the bead
                    # and that edge is important:
                    inner_bound[idx] = edge_idx + 1
                if (outermost[idx] >= self.edges[edge_idx]) and \
                        (outermost[idx] < self.edges[edge_idx + 1]):
                    # outer_bound is the index of the smaller edge of the
                    # outermost bin, since the intersection of the bead
                    # and that edge is important:
                    outer_bound[idx] = edge_idx
        self.bounds = np.column_stack((inner_bound, outer_bound))

    def _concentric_vol_shares(self):
        """
        Computes the portion of the volume of a bead in each of
        the concentric cylindrical-shell-like or spherical-shell-like
        bins by which the bead intersects along the radial direction in
        a cylindrical geometry.

        Notes
        -----
        For the sake of simplicity, the sphere-sphere intersection is
        also used for the sphere-cylinder intersection in the radial
        direction in the cylindrical and slit geometries. Hence, the
        concentric_vol_shares can be used in all these three different
        radial directions.

        This function can be used in calculation of the local volume
        fraction of beads in two different situation:
        1. The radial direction in the cubic geometry with the
        sphere_sphere_intersection function as the intersection_calculator.
        2. The radial direction in the cylindrical geometry with the
        sphere_cylinder_intersection function as the intersection_calculator.
        """
        self.volume_shares = {}
        for center_idx, bound_minmax in enumerate(self.bounds):
            self.volume_shares[center_idx] = {}
            # The intersection volume of the previous bin; for the innermost
            # bin, there is no previous bin, so this quantity is initiated
            # by 0:
            intersect_vol_previous = 0
            for edge_idx in range(bound_minmax[0], bound_minmax[1] + 2):
                # The difference between the intersection volume of a bead
                # with a bin edge with index idx and the previous bin is the
                # volume share of the bin (or center) with index idx-1:
                intersect_vol = \
                    sphere_sphere_intersection(
                        self.r_bead,
                        self.edges[edge_idx],
                        self.centers[center_idx]
                    )
                self.volume_shares[center_idx][edge_idx-1] = (
                    intersect_vol - intersect_vol_previous)
                intersect_vol_previous = intersect_vol

    def _vol_shares_type(self):
        """
        chooses how the volume_shares should be measured based on the
        given direction.

        Notes
        -----
        Currently, the `_vol_shares` method is implemented in the 'radial'
        direction in all the geometries (see the notes for
        `_concentric_vol_shares`) and the 'longitudinal' direction in the
        cylindrical geometry.
        """
        if self.direction == 'r':
            self._concentric_bounds()
            self._concentric_vol_shares()
        elif self.direction == 'z':
            self._consecutive_bounds()
            self._consecutive_vol_shares()
        else:
            raise ValueError(
                "'_vol_shares_type'"
                f" is not defined in the {self.direction} direction."
                )

    def _volume_fraction(self):
        """
        Computes the local volume fraction along the direction of interest
        in the given geometry.

        Notes
        -----
        It is assumed that all the particles have the same shape. The local
        volume fraction is normalized to give the integral of the local
        volume fraction  along the direction of interest in the region of
        interest. (See the explanation for the `_number_density` method and
        `SpatialDistribution` class).
        """
        n_centers = len(self.rho)
        self.phi = np.zeros(n_centers)
        for c_idx in range(n_centers):
            for idx, vol in self.volume_shares[c_idx].items():
                self.phi[c_idx] = self.phi[c_idx] + (self.rho[idx] * vol)


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
