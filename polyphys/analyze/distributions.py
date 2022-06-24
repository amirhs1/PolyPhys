"""
A sub-module to calculate various 1D spatial distributions of finite (not
point-like) particles (i.e. beads) in spherical, cylindrical, cartesian, or
polar coordinate system in a box, slit, or biaxial geomtery.

Finite particle can be hard or soft. The volume of a hard particle is set by
its geometrical shape. However, the volume of a soft particle is given by the
probability distribution function (PDF) of its mass around its center of mass
(COM) or its center of geoemetry (COG). Such a PDF can have some connection
with the soft potential defined for the interaction of this particle with
other particles in a system.

Currently, tehe following distribution are implemeted:
    1. The local number desnity of hard beads.
    2. The local volume fraction of hard beads.
"""
from typing import Callable, Dict, Tuple, Optional, Union
import numpy as np
import scipy.integrate as integrate
from polyphys.manage.parser import SumRule, TransFoci
import warnings

ParserT = Union[SumRule, TransFoci]


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
    1. a and b can be negative or postive, depending on the position
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
            " is not a postive number."
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
    redius `r1` has distance `d`from the other one with radius `r2`
    along axis x, so the vector form of distance is indeed (d,0,0)
    in cartesian coordinate system.

    Parameters
    ----------
    r1: float
        The radius of the one sphere.
    r2: flaot
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
    # By define rmax and rmin, we handle the situtations in which d = 0:
    rmax = max(r1, r2)
    rmin = min(r1, r2)
    # when one of the sphere has 0 radius, vol=0:
    if r1 == 0 or r2 == 0:
        vol = 0
    # Spheres are either tengential to each other or do not intersect.
    elif d >= rmin + rmax:
        vol = 0.0
    # Spheres intersect:
    # Small sphere resides completely in the large one.
    elif d <= rmax - rmin:
        vol = 4 * np.pi * rmin**3 / 3
    # Other scnarios:
    else:
        vol = (np.pi / (12*d)) * (rmax + rmin - d)**2 * (
            d**2
            + 2 * d * (rmax + rmin)
            - 3 * (rmin**2 + rmax**2)
            + 6 * rmin * rmax
        )
    return vol


class Distribution(object):
    hist_data: np.ndarray
    edges: Union[np.ndarray, None]
    centers: Union[np.ndarray, None]
    hist_info: ParserT
    integrand: Optional[Callable] = None
    normalized: bool = False
    """NOT FINISHED YET. Generates a Probability Distribution Function from
    a given histgrams by taking into account the sizes of multi-dimentional
    bins.

    Parameters
    ----------

    Attributes
    ----------
    self.histogram:
        histogram
    self.centers (numppy array):
        the centers of the bins: the index of histogram
    self.edges (numpy array): the edges of the bins of the same size. The
        range is [A-0.5*bin_size, B+0.5*bin_size]
    """
    def __init__(
        self,
        hists: np.ndarray,
        edges: Union[np.ndarray, None],
        centers: Union[np.ndarray, None],
        hist_info: ParserT,
        integrand: Optional[Callable] = None,
        normalized: bool = False,
    ):
        if isinstance(hists, np.ndarray):
            self.histogram = hists
        else:
            raise ValueError(
                f"'{hists}'"
                " is not a numpy.ndarray.")
        if isinstance(edges, np.ndarray) and (centers is None):
            self.edges = edges
            self.centers = 0.5 * (self.edges[:-1] + self.edges[1:])
        elif isinstance(centers, np.ndarray) and (edges is None):
            self.centers = centers
            # It is assume the bins are of the same size and the centers are
            # the acutal centers of the bins:
            _bin_size = centers[1] - centers[0]
            _edges = centers - 0.5 * _bin_size
            _edges = np.append(_edges, _edges[-1] + _bin_size)
            self.edges = _edges
            warnings.warn(
                "It is assumed that the bins are of the same size, so 'edges'"
                "are inferred from 'centers'.",
                UserWarning
                )
        elif isinstance(edges, np.ndarray) and isinstance(centers, np.ndarray):
            self.edges = edges
            self.centers = centers
        else:
            raise ValueError(
                f"Either '{edges}' and"
                f" '{centers}'"
                " are notnumpy.ndarray, or"
                " none of them are passed.")
        self.hist_info = hist_info
        self.integrand = integrand
        self.normalized = normalized
        if self.normalized is True:
            # the sum of rho is not equal to the bulk number density when r
            # approaches infinity, i.e. natom/vol_system. This arises from
            # the way we descritize the local number desnity.
            self.dist = self.dist / self.dist.sum()


class SpatialDistribution(object):
    hist_data: np.ndarray
    edges: np.ndarray
    hist_info: ParserT
    radius_attr: str
    geometry: str = 'biaxial'
    direction: str = 'r'
    normalized: bool = False
    """computes the local number density of any type of particle and the
    local volume fraction of bead (spheres) in 1D (cylinder with the
    periodic boundary condition (PBC) along the longitudinal axis (z axis)),
    2D (slit with the PBC along the x and y axes (radial direction)), and 3D
    (cube with the PBC along x, y, and z axes (spheric)).

    Parameters:
    histogram (pandas dataframe): a dataframe of an ensemble or
    ensemble-averaged or group simulation in which the index is the bin
    centers, the names of the columns are the name of simulations in an
    ensemble or the name of ensemble-averaged group or the name of
    simulation group, the columns are the frequenies of partilces of the
    type of interest (see radius_type) in each of the bins, and the number
    of columns is the number of simulations in a ensemble, one in an
    ensemble-averaged group, or the number of ensmebles (=the number of
    ensemble-averaged) groups an a simulations group.
    properties (pandas dataframe): the properties of each simulation/
    ensemble-averaged simulation/ensemble-averaged simulations of the
    simualtions/ensemble-averaged simulation/ensemble-averaged simulations
    in a ensemble/ensemble-averged group/simulation group.
    raduis_type (str): name of the column in the properties in which the size
    (or diameter) of the particles are tabled. The particle type is the same
    for all the particles that their frequencies are counted in histogram.
    geometry (str): the shape of the simulation box
    direction (str): the direction of interest in the geometry of interest.

    Attributes
    ----------
    self.histogram:
        histogram
    self.centers (numppy array):
        the centers of the bins: the index of histogram
    self.geometry:
        geometry
    self.direction:
        direction
    self.radius_type:
        raduis_type
    self.r_bead:
        radius of particle for whihch distribution are calculated.
    self.short_name_rho:
        the short name of number density distribution
    self.short_name_phi:
        the short name of volume fraction distribution
    self.r_bead:
        the radius of the particles that their type is set by
        raduis_type
    self.bin_size (float):
        the size of the equaly-sized bins.
    self.edges (numpy array): the edges of the bins of the same size. The
        range is [A-0.5*bin_size, B+0.5*bin_size]
    self.rho (pandas dataframe):
        the dataframe of the local number
        densities which has the same index, number of columns and names of
        columns as self.histogram. However, the column values are the local
        number densities.
    self.phi (pandas dataframe):
        the dataframe of the local volume
        fractions which has the same index, number of columns and names of
        columns as self.histogram. However, the column values are the local
        volume fractions.
    self.bounds: numpy.ndarray
        Indices of the lower and upper boundaries of a bead that resides at
        each of centers in `self.center`.
    self.volume_shares: A two-fold nested dict
        A nested dictionary in which keys are indices of `self.centers`
        (the positions of a bead's center) and the values are themselves
        dictionaries. In each of this dictionary, each key is one of the
        indices of `self.centers` within the bounds (inclusivily, i.e.
        [A,B] includes A and B as well given by `self.bounds` and the
        value of that key is the share of the bin of that center from the
        volume of the bead. The sum of all the values in one inner
        dictionary is equal to the volume of the bead.
    self._args: dict
        The arugments of eahc of the nine _integrands.

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

    While COGs or COMs of beads can be anywhere in thesystem, it is assumed
    that the COGs or COMs of all the beads reside in a bin are the center of
    that bin. While this assumption produces some error in the local volume
    distirbution, it significantly reduces the computaional cost. Moreover, the
    error decrease by increase the number of bins at a fixed volume of the
    systen. By this assumption, all the beads in a bin are equivalent and it
    is sufficent to do the two above steps only for one bead in a given
    direction in a given geomtery.

    To find the upper and lower bounds on edges, the center of a bead with
    radius `self.r_bead` is placed at one of `self.centers`. Using the
    above assumption, for a given `self.r_bead`, a pair of lower and upper
    bounds (`self.bounds`) can be found for each center. It is assumed
    that a lower or upper bound is the leftmost, rightmost, innermost, or
    outermost edges smaller than the boundary of bead in a given direction.

    There are four bin schemes:

        'concenteric':
            bins, edges, or centers are non-negative strictly increasing
            arrays, creating concentirc circles, spherical shells, or
            cylindircal shells along the radial direction in the polar,
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
            cylindircal coordinate systems or the spherical coordiante
            system respectively, or generating concentric.

        'neighboring:
            bins, edges, or centers are strictly increasing but bounded
            by the period of the direction of interest, thus generating
            sectors created by two neighboring azimuth (conic) planes
            along the azimuth direction in the spherical coordiante system.

    For each of this bin type, an algorithm is defined below to find the
    lower and upper bounds.

    Finding bond edges along a given direction is a 1D problem. Moreover,
    the following assumtions about positions of beads (particles), bin
    centers and edges are used in this class:
    1. If a bead is in a bin, its position is assumed to be the center
        of that bin, not its actual position.
    2. len(self.edges) = len(self.centers) + 1
    3. For a bin, we have: left_edge <= center < right_edge.
    4. There is only one center between two consecutive edges.
    5. The edge and center arrays are strictly increasing funcions of
        their indices.
    6. The bins are equally spaced, son bin_size is the same between any
        two consecutive edges.

    To coserve the total volume of the beads, the periodic boundary condition
    (PBS) is imposed if required (For instance, the z direction in the
    biaxial geometry). The PBC is impose in both step: finding the bounds
    and calculating the volume shares.
    """
    _directions = {
        'box': ['r', 'theta', 'phi'],
        'slit': ['r', 'phi', 'z'],
        'biaxial': ['r', 'phi', 'z']
    }
    _integrands = {
        'box': {
            'r': lambda r, const: 4 * const * np.pi * r**2,
            # constant is redundant and merely defined to consistantly use
            # of 'args' parameter of "scipi.integral.quad" among integrands.
            'theta': lambda theta, dcyl: np.pi * dcyl**3 * np.sin(theta) / 12,
            'phi': lambda theta, dcyl: dcyl**3 / 12
        },
        'slit': {
            'r': lambda r, lcyl: 2 * np.pi * lcyl * r,
            'phi': lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2,
            'z': lambda z, dcyl: 0.25 * np.pi * dcyl**2
        },
        'biaxial': {
            'r': lambda r, lcyl: 2 * np.pi * lcyl * r,
            'phi': lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2,
            'z': lambda z, dcyl: 0.25 * np.pi * dcyl**2
        }
    }
    _fullnames = {
        'box': {
            'r': 'radial',
            'theta': 'azimuthal',
            'phi': 'polar'
        },
        'slit': {
            'r': 'radial',
            'phi': 'polar',
            'z': 'longitudinal'
        },
        'biaxial': {
            'r': 'radial',
            'phi': 'polar',
            'z': 'longitudinal'
        }
    }

    def __init__(
        self,
        hist_data: np.ndarray,
        edges: np.ndarray,
        hist_info: ParserT,
        radius_attr: str,
        geometry: str = 'biaxial',
        direction: str = 'r',
        normalized: bool = False,
    ):
        if isinstance(hist_data, np.ndarray):
            self.histogram = hist_data
        else:
            raise ValueError(
                f"'{hist_data}'"
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
        self.hist_info = hist_info
        self.r_bead = 0.5 * getattr(self.hist_info, radius_attr)
        self.normalized = normalized
        self.short_name_rho = \
            self._fullnames[self.geometry][self.direction]+'Rho'
        self.short_name_phi = \
            self._fullnames[self.geometry][self.direction]+'Phis'
        self._vol_shares_type()
        self._set_args()
        self._number_density()
        self._volume_fraction()
        if self.normalized is True:
            # the sum of rho is not equal to the bulk number density when r
            # approaches infinity, i.e. natom/vol_system. This arises from
            # the way we descritize the local number desnity.
            self.rho = self.rho / self.rho.sum()
            # time averaging: the sum of histograms = natoms * nframes.
            # normalization: the sum of the number density is now 1.
            # the sum of phi is not equal to the bulk volume fraction when r
            # approaches infinity, i.e, (natoms*vol_atom)/vol_system. This
            # arises from the way we descritize the local number desnity.
            self.phi = self.phi / self.phi.sum()

    def _consecutive_bounds(self):
        """
        finds the indices of smallest and largest bin edges by which a
        bead, that resides at the center of that bin, intersects in a
        'consecutive' bin scheme.

        Attributes
        ----------
        self.bounds: np.ndarray
            Indices of the leftmost and rightmost boundaries of a bead
            that resides at each of centers in `self.center`. For each
            each center, the following four varaibles are used to
            generate `self.bounds`:

            leftmost:
                The value of the center of the bin in which the leftmost
                boundary of a bead is located.

            left_bound:
                The index of the 'smaller' edge of the bin in which the
                leftmost boundary of a bead is located.

            rightmost:
                The value of the center of the bin in which the rightmost
                boundary of a bead is located.

            right_bound:
                The index of the 'smaller' edge of the bin in which the
                rightmost boundary of a bead is located.
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
                    # of the rightmost bin
                    right_bound[idx] = edge_idx
        self.bounds = np.column_stack((left_bound, right_bound))

    def _consecutive_vol_shares(self):
        """
        computes the portion of the volume of a bead (a sphere) in each of
        the consecutive disk-like bins by which the bead intersects along
        the longitudinal direction with the PBC in a cylindrical geometry.

        To-do List
        ----------
        Vectorize the for-loops.

        Attributes
        ----------
        self.volume_shares: A two-fold nested dict
            A nested dictionary in which keys are indices of `self.centers`
            (the positions of a bead's center) and the values are themselves
            dictionaries. In each of this dictionary, each key is one of the
            indices of `self.centers` within the bounds (inclusivily, i.e.
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
            # A bead near the end of a system along a peroidic direction
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
                    # Distacnes of a bin's left and right edges of from the
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
                # Reset center to its inintal value to calculate the rest
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
        finds the indices of smallest and largest bin edges by which a
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
        computes the portion of the volume of a bead in each of
        the conncentric cylindrical-shell-like or spherical-shell-like
        bins by which the bead intersects along the radial direction in
        a cylindrical geometry.

        To-do List
        ----------
        Vectorize the for-loops.

        Notes
        -----
        For the sake of simpilicity, the sphere-sphere interesction is
        also used for the sphere-cylinder intersection in the radial
        direction in the cylindrical and slit goemteries. Hence, the
        concentric_vol_shares can be used in all these three different
        radial directions.

        This function can be used in calculation of the local volume
        fraction of beads in two different situation:
        1. The radial direction in the cubic geometry with the
        sphere_sphere_intersection function as the intersection_calculator.
        2. The radial dirction in the cylindrical geometry with the
        sphere_cylinder_intersection function as the intersection_calculator.
        """
        self.volume_shares = {}
        for center_idx, bound_minxax in enumerate(self.bounds):
            self.volume_shares[center_idx] = {}
            # The intersection volume of the previous bin; for the innermost
            # bin, there is no previous bin, so this quanyoyu is initiated
            # by 0:
            intersect_vol_previous = 0
            for edge_idx in range(bound_minxax[0], bound_minxax[1] + 2):
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

    def _set_args(self):
        """
        sets the arguments for the integrads along different directions
        in different geometries.
        """
        self._args = {
            'box': {
                'r': (1, ),
                # For concentric spherical shells, the constant, i.e. 1,
                # is redundant and merely defined to make the use of args
                # parameter of scipi.integral.quad function consistent among
                # integrands.
                'theta': (self.hist_info.lcyl, ),
                # In a box cubic or free space, the radius of the space
                # is half of the length of simulation box.
                'phi': (self.hist_info.lcyl, ),
                # In a box cubic or free space, the radius of the space
                # is half of the length of simulation box.
            },
            'slit': {
                'r': (self.hist_info.lcyl, ),
                'theta': (self.hist_info.lcyl, self.hist_info.dcyl, ),
                'z': (self.hist_info.dcyl, )
            },
            'biaxial': {
                'r': (self.hist_info.lcyl, ),
                'theta': (self.hist_info.lcyl, self.hist_info.dcyl, ),
                'z': (self.hist_info.dcyl, )
            }
        }

    def _number_density(self):
        """
        calculates the local number density along the given direction in
        the given geometry.

        If `self.normalized=True`, the the local number density is
        normalized to give the area under the curve equal to one.

        Parameters
        ----------
        col_name: the name of column for which the number density is
        calculated.

        Notes
        -----
        The number density in each simulation is an average over the number
        densities collected every X time steps, so there are N=L/X
        measurements of the local number desnity in each simulation where L
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
        self.rho = self.histogram / bin_vols

    def _volume_fraction(self):
        """
        computes the local volume fraction along the direction of interest
        in the given goemetry.

        Notes
        -----
        It is assumed that all the particles have the same shape. The local
        volume fraction is normalized to give the integral of the local
        volume fraction  along the direction of interest in the region of
        interest. (See the explnation for the `_number_density` method and
        `Distribution` class).
        """
        n_centers = len(self.rho)
        self.phi = np.zeros(n_centers)
        for c_idx in range(n_centers):
            for idx, vol in self.volume_shares[c_idx].items():
                self.phi[c_idx] = self.phi[c_idx] + (self.rho[idx] * vol)


def distributions_generator(
    wholes: Dict[str, np.ndarray],
    bin_edges: Dict[str, np.ndarray],
    group: str,
    species: str,
    geometry: str,
    direction: str,
    parser: Callable,
    save_to: Optional[str] = None,
    normalized: bool = False
) -> Tuple[Dict, Dict]:
    """
    generates the local number density (rho) and volume fraction (phi)
    of the bead_name woth bead_type column name in properties
    dataframe along the direction of interest in the geometry of interest.

    Caution:
    For the sumrule problem in the cylindrical goemetry, a simulation group
    is a collection of simulations that all have the same values for the
    number of monomers, the diameter of cylinder, and the size of crowders
    (assuming size of monomers is 1). An ensemble (usually with M number of
    simulations) is a collection of themodynamically-equivalent simulations
    that all have the same values for the number of monomers, the diameter
    of cylinder, the size of crowders, and the same number of crowders
    (assuming size of monomers is 1).  In standard statitatical mechanical
    approach, an ensmeble is equivalent a simulation, but here we used it
    to reffer to all the thermodynamically-equivalent simulations.
    In an ensemble-average simulation group, each ensemble is replaced with
    the average of all its simulation, so if we have M simulation in an
    ensemble and P ensembles in a group, then we have M*P simulations in a
    simulation group and P ensemble-averaged simulations in a
    ensemble-averaged simulation group.

    Parameters
    ----------
    histogram: dict
        A dictionary of an ensemble or ensemble-averaged or group
        simulation in which the keys are the name of ensembles/
        ensemble-averaged groups/groups and the keys of the histogram
        dataframes. The names of the columns in each dataframe are the name of
        simulations in an ensemble or the name of ensemble-averaged group or
        the name of simulation group, the columns are the frequenies of
        partilces of the type of interest (see radius_type) in each of the
        bins, and the number of columns is the number of simulations in a
        ensemble, one in an ensemble-averaged group, or the number of
        ensmebles (=the number of ensemble-averaged) groups an a simulations
        group.
    properties: pd.dataframe
        The properties of each simulation/ensemble-averaged
        simulation/ensemble-averaged simulations of the
        simualtions/ensemble-averaged simulation/ensemble-averaged simulations
        in a ensemble/ensemble-averged group/simulation group.
    raduis_type: str
        the name of the column in the properties in which the size (or
        diameter) of the particles are tabled. The particle type is the same
        for all the particles that their frequencies are counted in histogram.
    geometry: str
        the shape of the simulation box
    direction: str
        the direction of interest in the geometry of interest.
    bead_name: the name of paticle type.
    parser: Callable
        A class from 'PolyPhys.manage.parser' moduel that parses filenames
        or filepathes to infer information about a file.
    save_to: str
        path to which the output saved.
    normalized: bool
        whether normalized the distributions or not.

    Return
    ------
    densities: dict
       A dictionary in which the keys are 'whole' names and
       the values are the local number densities.
    vol_fractions: dict
    A dictionary in which the keys are 'whole' names and
       the values are the local volume fractions.
    """
    densities = {}
    vol_fractions = {}
    radius_attrs = {
        'Mon': 'dmon',
        'Crd':  'dcrowd'
    }
    for whole, histogram in wholes.items():
        whole_info = parser(
            whole,
            geometry,
            group,
            'whole',
            ispath=False
        )
        distributions = SpatialDistribution(
            histogram,
            bin_edges[whole],
            whole_info,
            radius_attr=radius_attrs[species],
            geometry=geometry,
            direction=direction,
            normalized=normalized
        )
        densities[whole] = distributions.rho
        vol_fractions[whole] = distributions.phi
        if save_to is not None:
            filename = whole + '-' + group + '-' + direction
            np.save(
                filename + 'Rho' + species + '.npy', distributions.rho
            )
            np.save(
                filename + 'Phi' + species + '.npy', distributions.phi
            )
    return densities, vol_fractions


def entropic_energy(histo_collections, beta=1.0):
    """
    this is the effective free energy potential as the function of the
    end-to-end distance. Here, we use the end-to-end distance as the
    action parameters and write a Langevin equation for that.

    parameters:
    histo_collections: radial distribution function or the probability
    distribution function of *the magnitude of the end-to-end vector *
    (or *the end-to-end distance*) between R and R+dR.
    beta: k_B*T is the inverse of the thermal energy.

    return:
    the dimensionlessfree energy potential of the end-to-end distance
    -beta * Ln(P(r)*4pi*r^2) = -beta * Ln(histo_collections) where
    beta=k_BT is set 1.0.
    """
    histo_collections = histo_collections / np.sum(histo_collections)
    free_energy = -1.0 * (
        np.log(histo_collections) - np.log(np.sum(histo_collections))) / beta
    return free_energy


def looping_p(histo_collections, bin_edges, Rmax, Rmin):
    """looping entropy P_L is defined as the integral from Rmin to Rmax over
    P(R)*4pi*R^2*dR; since P(R)*4pi*R^2 is equivalent to hitso_collection[i]
    for the bin between bin_edges[i] and bin_edges[i+1] we have
    P_L = integral from Rmin to Rmax over P(R)*4pi*R^2*dR ~ sum from Rmin
    to Rmax over P(bin_center_i)*4pi*bin_centers_i^2*(bin_edges_i+1 -
    bin_edges_i)=sum from Rmin to Rmax over hitso_collection_i*(bin_edges_i+1
    - bin_edges_i)

    Since the sizes of the bins are equal ((bin_edges_i+1 - bin_edges_i)
    =constant), the effect of (bin_edges_i+1 - bin_edges_i) is cancell out
    upon normalization.

    Inputs:
    histo_collections: a numpy array of the histograms of the Flory radius
    or the end-to-end distance during the whole run (the total number of
    data points is equal to the total number of time steps.)
    bin_edges: the edges of the bins used in histograming.
    Rmax: the minimum center-to-center distance the two end-of-chain
    monomers can have; Since they have the same size, this is usually equal
    to their diameter (or size).
    Rmin: this is the size of the crowders that can reside between two
    monomers. When the center-to-center distance between two monomers is
    less than Rmax+Rmin, no crowders can be between two monomers and the
    depletion force is non-zero.

    Returns:
    The probability of looping.
    """
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    histo_collections = np.multiply(
        histo_collections, (bin_edges[1:] - bin_edges[:-1]))
    histo_collections = histo_collections / np.sum(histo_collections)
    looped_probablity = 0
    for i in range(len(bin_centers)):
        if bin_centers[i] <= (Rmax + Rmin):
            looped_probablity += histo_collections[i]
    return looped_probablity
