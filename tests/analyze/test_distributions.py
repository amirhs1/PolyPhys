import warnings
import pathlib
from glob import glob
import numpy as np
import pandas as pd
import scipy.integrate as integrate

from ..manage import organizer

# Analysis based on the dtat extracted from trajectories; in other words, indirect analysis of trajectories:
class Distributions():
    """
    Distributions computes the local number density of any type of particle and the local volume fraction of bead (spheres) in 1D (cylinder with the eperiodic boundary condition (PBC) along the longitudinal axis (z axis)), 2D (slit with the PBC along the x and y axes (radial direction)), and 3D (cube with the PBC along x, y, and z axes (spheric)).
        Here it is assumed that particles are hard so their volume of intersections with bins are numerically calculated in full details. However, particles can also be soft, so their volume can follw a particle distribution such as the normal distribution. Indeed, the soft particle approcimation is a more reasonbale choice for Lennard-Jones particles. If the sof particle distribution is desirable, this class should be extended to have a functionality to take an arbitrary function as the volume distribution function for particles. The soft-particle approximation is indeed similar to methods in other branches of STEM to smooth a histogram.
    There are three types of bins: concenteric bins (and thus concentric edges) along the radial direction, consecutive bins along the longitudinal direction, and cyclic consecutive bins along the polar and azimuthal directions. This method correctly finds the bounds in the first two cases.
    The histogram (frequecny of occurance) and number denisty are functions of one bin/bin center/a pair of consecutive edges:
    If H_i the freiqency of particle in in bin or at bin center i, then rho_i=H_i/V_i is the number density and V_i is the volume of the bin i. 
    However, the volume fraction is a function of several bins. More precisely, the number of bins is set by the ratio of the bead size (diameter) R divided by the bin size h: M_j<=R/h is the number of bins with which the bead intersect/overlaps where j is the index of a bin or bin center. the smaller relation (<) holds when the bead which is close to the boundaries (the left and right boundary along longitudinal direction or the outermost boundary along the radial direction). In other cases, the equality (=) holds. In other words, the volume of a bead is shared amoing M_j bins, so it is need to find the volume of intersection of the bead and a bin. For the sake of similipcity, it is assumed all the particles that are in bin i with frequency H_i are located at the bin center, so we have a bin/test bead that is located at the center of bin. There are N bins or bin centers, test or bin bead, frequnecies H_i,  number densities rho, and volume fraction phi while there is  N+1 bin edges.
    The volume fraction phi_i of bin i is the sum over all the shares that bin i have from the volumes of N bin/test beads; in other words, phi_i=The sum over rho_j*v_j where P_i is the number of shares (the number of terms in the sum) that bin i have from the volumes of N bin/test beas, v_j is the volume share of the bin_i from bin/test bead j and rho_j is the number density of the bin j in which the bin/test bead j resides.
    Due to the above-described assumption and the one-dimensional nature of the problem (frequencies, densities, and volume fractions are measured along one direction), it is needed to find the volume share v_j once for each bin/test bead.
    The distribution functions (histogras, number densities, and volume fractions) are averaged twice: First, over all the measurements in one simulation (time-averaging), this steps result in the distributions of an individual simulation. Then, the time averages distribution of interest is avaraged over all the thermoynamiclly-equivalent simulation in an ensemble and an ensemble-averaged distribution is generated.
    len(histogram.index)=len(bin_size)=len(bin_vols)=len(rho)=len(bin_edges)-1
    
    Issue: Resolved on 20211016
    The problem with some choice of bin size such a 0.2a for monomer size a=1.0 in the local volume fraction of monomers, or 0.4a for crowder size a_c=6.0a with monomer size a=1.0 in the local volume fraction of crowders has not been resolved yet.
    
    Caution:
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    
    Parameters:
    histogram (pandas dataframe): a dataframe of an ensemble or ensemble-averaged or group simulation in which the index is the bin centers, the names of the columns are the name of simulations in an ensemble or the name of ensemble-averaged group or the name of simulation group, the columns are the frequenies of partilces of the type of interest (see radius_type) in each of the bins, and the number of columns is the number of simulations in a ensemble, one in an ensemble-averaged group, or the number of ensmebles (=the number of ensemble-averaged) groups an a simulations group. 
    properties (pandas dataframe): the properties of each simulation/ensemble-averaged simulation/ensemble-averaged simulations of the simualtions/ensemble-averaged simulation/ensemble-averaged simulations in a ensemble/ensemble-averged group/simulation group.
    raduis_type (str): the name of the column in the properties in which the size (or diameter) of the particles are tabled. The particle type is the same for all the particles that their frequencies are counted in histogram. 
    geometry (str): the shape of the simulation box
    direction (str): the direction of interest in the geometry of interest.
    
    Properties:
    self.histogram: histogram
    self.centers (numppy array): the centers of the bins: the index of histogram
    self.properties: properties
    self.geometry: geometry
    self.direction: direction
    #self.radius_type: raduis_type
    self.r_particle: radius of particle for whihch distribution are calculated.
    self.short_name_rho: the short name of number density distribution
    self.short_name_phi: the short name of volume fraction distribution
    self.r_particle: the radius of the particles that their type is set by raduis_type
    self.bin_size (float): the size of the equaly-sized bins.
    self.edges (numpy array): the edges of the bins of the same size. The range is [A-0.5*bin_size, B+0.5*bin_size]
    self.rho (pandas dataframe): the dataframe of the local number densities which has the same index, number of columns and names of columns as self.histogram. However, the column values are the local number densities.
    self.phi (pandas dataframe): the dataframe of the local volume fractions which has the same index, number of columns and names of columns as self.histogram. However, the column values are the local volume fractions.
    self.particle_bounds: an array of dimension (len(centers),2) containing the pairs of lowest and highest bounds (the edges) at which the bead of radius r_particle located at centers interesect with.
    self.volume_shares (nested dict): A two-fold nested dict in which the (outmost) keys of the outmost dict are the centers and the (outmost) values of that are dicts with an unequal numbers of keys. The keys of the inner dict are the edge indexes between the atom_bounds (inclusivily, i.e. [A,B] includes A and B as well) and the values of that are the volume of intersection between the bead and the bin for which the edge is the left or inner index. This is a representation of the volume_share
        self.volume_shares={center1:{edge_index1: intersection_volume1, edge_index2: intersection_volume2 ...} , ...}
    self._args (dict): the arugments of eahc of the nine _integrands.
    
    Requirements:
    panads, numpy, integrate function from scipy.integrate.
    """
    _geometries = ["cubic","slit","cylindrical"] # different types of simulation boxs: cubic (free space or pbc: periodic boundary condition in all three directions), slit (pbc in x and y directions), cylinder (pbc in z direction)
    _directions = {
        "cubic":["radial","polar","azimuthal"],
        "slit":["radial","polar","longitudinal"],
        "cylindrical":["radial","polar","longitudinal"]
    }
    _integrands = {
        "cubic":{"radial" : lambda r, const: 4 * const * np.pi * r**2, # cocentric spherical shells; constant is redundant and merely defined to make the use of args parameter of scipi.integral.quad function consistent among integrands.
                 "polar": lambda phi, dcyl: dcyl**3 / 12, # spherical sections of diameter dcyl
                 "azimuthal": lambda theta, dcyl: np.pi * dcyl**3 * np.sin(theta) / 12 # spherical wedges of diameter dcyl
                },
        "slit":{
            "radial" : lambda r, lcyl: 2 * np.pi * lcyl * r, # cocentric cyliderical shells with length lcyl
            "polar" : lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2, # cylindrical sectors of length lcyl and diameter dcyl
            "longitudinal": lambda z, dcyl: 0.25 * np.pi * dcyl**2 # disks of diamter dcyl
        },
        "cylindrical":{
            "radial": lambda r, lcyl: 2 * np.pi * lcyl * r, # cocentric cyliderical shells with length lcyl
            "polar": lambda theta, lcyl, dcyl: 0.25 * lcyl * dcyl**2, # cylindrical sectors of length lcyl and diameter dcyl
            "longitudinal": lambda z, dcyl: 0.25 * np.pi * dcyl**2 # disks of diamter dcyl
        }
    }
    _short_names = {
        "cubic":{"radial" : "r",
                 "polar": "theta",
                 "azimuthal": "phi"
                },
        "slit":{
            "radial" : "r",
            "polar" : "theta",
            "longitudinal": "z"
        },
        "cylindrical":{
            "radial": "r",
            "polar": "theta",
            "longitudinal": "z"
        }
    }
    
    def __init__(self, histogram, properties, r_particle, geometry, direction, normalized=True):
        if isinstance(histogram, pd.DataFrame):
            self.histogram = histogram
            self.centers = np.around(histogram.index.to_numpy(),decimals=2) # since centers are index it is important to have them round. decimals=2 work well for a_c/a<1 and a_c>a>=1 where a is monomer diamter and a_c is  crowder diamter.
        else:
            raise TypeError(f"'{histogram}' is not a Pandas Dataframe. Currently, Pandas Dataframes are only supported.")    
        if isinstance(properties, pd.DataFrame):
            self.properties = properties
        else:
            raise TypeError(f"'{properties}' is not a Pandas Dataframe. Currently, Pandas Dataframes are only supported.")
        if geometry in self._geometries:
            self.geometry = geometry
        else:
            geomteries_string = "'" + "', '".join(self._geometries) + "'"
            raise ValueError(f"'{geometry}'"
                             " is not a valid geometry for the simulation box. Please select one of "
                             f"{geomteries_string} geometries.")
        if direction in self._directions[self.geometry]:
            self.direction = direction
        else:
            directions_string = "'" + "', '".join(self._directions[self.geometry]) + "'"
            raise ValueError(f"'{direction}'"
                             " is not a valid direction for "
                             f"'{self.geometry}' geometry. Please select one of "
                             f"{directions_string} directions.")
        #self.radius_type = radius_type
        #self.r_particle = 0.5 * self.properties[self.properties.filename==filename][self.radius_type].values[0] # Assumed the sizes of the particle of interest are the same for all the columns.
        #self.r_particle = 0.5 * self.properties[self.radius_type].values[0] # Assumed the sizes of the particle of interest are the same for all the columns.
        self.r_particle = r_particle
        self.normalized = normalized
        self.short_name_rho = self._short_names[self.geometry][self.direction]+'Rhos'
        self.short_name_phi = self._short_names[self.geometry][self.direction]+'Phis'
        self._edges_from_centers()
        self._initiate_distribution()
        self._vol_shares_type()
        self._run()
    
    def _edges_from_centers(self):
        """
        _edges_from_centers creates bin edges from bin centers, assuming:
        1. all the bins have the same size, so the distance between two consecutive edges and two consecutive centers are equal.
        2. The centers linearly increase over range is [A,B]. ] means including B.
        
        Caution:
        len(edges) = len(centers) + 1
        """
        
        self.bin_size = np.around(self.centers[1] - self.centers[0],decimals=2) # bin_size: it shound be rounded otherwise, the bin centers (the histograms indexes) and bin edges does not have the same noumber of meaningful digits.
        self.edges = self.centers - 0.5 * self.bin_size # all the edges except the last
        self.edges = np.append(self.edges,self.centers[-1]+ 0.5*self.bin_size) # the last edge      
    
    def _initiate_distribution(self):
        """
        _initiate_distribution creates empty dataframe with the same index and columns as the given one.
    
        """
        self.rho = pd.DataFrame(index=self.histogram.index, columns=self.histogram.columns) # number density
        self.phi = pd.DataFrame(index=self.histogram.index, columns=self.histogram.columns) # local volume fraction
    
    def spherical_segment(self, r, a, b):
        """
        spherical_segment computes the volume of a spherical segment. This function is written 
        based on equation 3 in the following webpage:
        https://mathworld.wolfram.com/SphericalSegment.html

        Caution: 
        1. a and b can be negative or postive, depending on the position with respect to the center of sphere in the range [-r,r].
        2. When a=r or b=r, we have a spherical cap.

        Parameters:
        r (float): Radius of the sphere
        a (float): Distance of the first base from the center of the sphere
        b (float): Distance of the second base from the center of the sphere

        Returns:
        V (float): Volume of the spherical segment

        Requierment: 
        Numpy.
        """
        if r <= 0:
            raise ValueError(" r should be within a non-zero postive number.")
        
        lower = min(a,b) # The lower bound of the volume integral
        upper = max(a,b) # The upper bound of the volume integral
        # always upper>lower since we assume the axis over which left, center, right are located points to right.       
        if np.abs(upper) >= r :
            #print("The absolut value of the upper bound is larger r, so it is set to r.")
            upper = upper * r / np.abs(upper)
        elif np.abs(lower) >= r :
            #print("The absolute value of the lower bound is larger r, so it is set to r.")
            lower = lower * r / np.abs(lower)
        V = np.pi * (r**2*(upper-lower) - (upper**3 - lower**3)/3)
        if (np.abs(lower) >= r and np.abs(upper) >= r) and lower * upper > 0:
        #if lower * upper > r**2:
            #print(" Both a and b are either smaller than -r or larger than r.")
            V = 0
        return V

    def _consecutive_bounds(self):
        """
        _find_bounds_consecutive finds the lowest and highest bin edges by which a spherical bead residing at the center of a bin intersects. There are three types of bins: concenteric bins (and thus concentric edges) along the radial direction, consecutive bins along the longitudinal direction, and cyclic consecutive bins along the polar and azimuthal directions. This method correctly finds the bounds in the first two cases. This method find the bounds in a consecutive-bins shceme along the longitudinal direction.
        
        Cuation: 
        This is a 1D problems since the bins, bin centers and bin edges are along the same axis. In reality, the center of a spherical bead can be anywhere in a bin, but it is placed at the center of the bin along the axis of interest.
        The center of a particle can be in any of the bins, so 
        While this assumption makes the value of the local volume fraction slightly inaccurate, it singnificantly reduce the computational cost of measure the local volume fraction in a system with a large number of particles.
        The positions and r_concentrics have some relations:
        1. len(edges) = len(centers) + 1
        2. There is only one center between two consecutive edges.
        3. All the arrays are increasing funcions of their indices.
        
        leftmost/leftmost_idx is the value/index of the bin center in which the leftmost boundary of a particle is located.
        rightmost/rightmost_idx is the value/index of the bin center in which the rightmost boundary of a particle is located.
        
        Instead of the leftmost/rightmost pairs (which are more appropriate for the longitudinal direction), the innermost/outermost can also be used.
        """
        box_length = self.edges[-1] - self.edges[0] # length of simulation domain along the direction of interest.
        leftmost = self.centers - self.r_particle # The minimum distance of the r_atoms perimeter from the origin
        leftmost = np.where(leftmost < self.edges[0],leftmost + box_length,leftmost)
        leftmost_idx = np.zeros(len(leftmost),dtype=int) # Initiate the leftmost bound with the lowest possible bound
        rightmost = self.centers + self.r_particle # The maximum distance of the r_atoms perimeter from the origin
        rightmost = np.where(rightmost > self.edges[-1],rightmost - box_length,rightmost)
        rightmost_idx = (len(rightmost)-1) * np.ones(len(rightmost),dtype=int) # Initiate the rigtmost bound with the highest possible bound
        for idx, leftmost_value in enumerate(leftmost):
            for edge_idx in range(len(self.edges[:-1])):
                if (leftmost_value >= self.edges[edge_idx]) and (leftmost_value < self.edges[edge_idx+1]): # the index of the leftmost bin (or the index of the **left edges** of leftmost bin) is set as the index of the bin in which the leftmost side of the bead is located.
                    leftmost_idx[idx] = edge_idx
                if (rightmost[idx] > self.edges[edge_idx]) and (rightmost[idx] <= self.edges[edge_idx+1]): # the index of the rightmost bin (or the index of the **left edge** of rigthmost bin) is set as the index of the bin in which the rightmost side of the bead is located. Keep the difference in <= ith the leftmost in mind.
                    rightmost_idx[idx] = edge_idx
        self.particle_bounds = np.column_stack((leftmost_idx,rightmost_idx)) 

    def _consecutive_vol_shares(self): # make this method better -- see the commented effort:
        """
        _concentric_vol_shares computes the portion of the volume of a bead (a sphere) in consecutive disk-like bins along the longitudinal direction in a cylindrical geometry. The center of the test particle with radius self.r_particle is placed at the center of each bin. Depending on its radius, the particle can intersect with more than one bins. so its volume can contributes to the local volume fraction in more than one bin. In the self.volume_shares nested dictionary below,  self.centers are keys, each showing the bin at which the center of particle resides and dictionaries of the volume shares are the values. In inner dictionary, the bin center of a bin by which the particle intersects is the key and the volume of intersection of the particle with the bin is the value.
        
        In this algorithm, the total volume of a particle is conserved.

        Caution:
        The centers and edges have these relations:
        1. len(edges) = len(centers) + 1
        2. There is only one center between two consecutive edges.
        3. All the arrays are increasing funcions of their indices.
        4. To solve the problem with PBC in z direction, there is if-elif-else statement below in which the if part ensures the volume shares in the bins, that are close to the left side of the simulation box, are counted correctly and sum up to the volume of the particle. The elif part does the same for the bins close to the right side of the simulation box. the else part handles the rest of bins.

        Parameters:
        intersection_calculator: the function computes the volume of intersection between the bead with r_particle as its radius centerred at centers  and the consecute disk with diameter dcyl centered at edges. It is the self.spherical_segment in this method.

        Return: 
        volume_shares (a thwo-fold nested dict): A three-fold nested dict in which the (outmost) keys of the outmost dict are the centers and the (outmost)values of that are the dicts with an unequal number of keys. The keys of the inner dict are the edge indexes between the atom_bounds (inclusivily, i.e. [A,B] includes A and B as well) and the values of that are the volume of intersection between the bead and the bin for which the edge is the left or inner index.
        volume_shares={center1:{edge_index1: intersection_volume1, edge_index2: intersection_volume2 ...} , ...}  
        
        Effort to make this method better:
               for center_idx, bounds_minmax in enumerate(self.particle_bounds):
            self.volume_shares[center_idx] = {}
            # share of the middle bins from the volume of the sphere:
            # edge value is an increasing function of edge index
            if (bounds_minmax[0] > bounds_minmax[1]) and (center_idx <= len(self.centers)//2):
                lower_bound = bounds_minmax[0]
                upper_bound = len(self.edges)-1
                lower_bound_2 = 0
                upper_bound_2 = bounds_minmax[1]+1
                center_1 = self.centers[center_idx] + box_length # Due to the PBC, the bin center is moved to right side, so the volume of intgersection can be calculated correctly
                center_2 = self.centers[center_idx]
            elif (bounds_minmax[0] > bounds_minmax[1]) and (center_idx > len(self.centers)//2):
                lower_bound = bounds_minmax[0]
                upper_bound = len(self.edges)-1
                lower_bound_2 = 0
                upper_bound_2 = bounds_minmax[1]+1
                center_1 = self.centers[center_idx]           
                center_2 = self.centers[center_idx] - box_length # Due to the PBC, the bin center is moved to right side, so the volume of intgersection can be calculated correctly.
            else:
                lower_bound = bounds_minmax[0]
                upper_bound = bounds_minmax[1]+1
                center= self.centers[center_idx] 
            for edge_idx in range(bounds_minmax[0],bounds_minmax[1]+1,1):
                    left_distance = self.edges[edge_idx]-self.centers[center_idx]
                    right_distance =  self.edges[edge_idx+1]-self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    #if np.abs(right_distance) >= self.r_particle: 
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
        """
        box_length = self.edges[-1] - self.edges[0]
        self.volume_shares = {}
        for center_idx, bounds_minmax in enumerate(self.particle_bounds):
            self.volume_shares[center_idx] = {}
            # share of the middle bins from the volume of the sphere:
            # edge value is an increasing function of edge index
            
            if (bounds_minmax[0] > bounds_minmax[1]) and (center_idx <= len(self.centers)//2): 
                for edge_idx in range(bounds_minmax[0],len(self.edges)-1,1):
                    center = self.centers[center_idx] + box_length # Due to the PBC, the bin center is moved to right side, so the volume of intgersection can be calculated correctly.
                    left_distance = self.edges[edge_idx] - center
                    right_distance =  self.edges[edge_idx+1] - center # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
                
                for edge_idx in range(0, bounds_minmax[1]+1,1):
                    left_distance = self.edges[edge_idx] - self.centers[center_idx]
                    right_distance =  self.edges[edge_idx+1] - self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)

            elif (bounds_minmax[0] > bounds_minmax[1]) and (center_idx > len(self.centers)//2):
                for edge_idx in range(bounds_minmax[0],len(self.edges)-1,1):
                    left_distance = self.edges[edge_idx]-self.centers[center_idx]
                    #if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    right_distance =  self.edges[edge_idx+1]-self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
                
                for edge_idx in range(0, bounds_minmax[1]+1,1):
                    center = self.centers[center_idx] - box_length # Due to the PBC, the bin center is moved to left side, so the volume of intgersection can be calculated correctly.
                    left_distance = self.edges[edge_idx] - center
                    right_distance =  self.edges[edge_idx+1] - center # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    #if np.abs(right_distance) >= self.r_particle: 

                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
            
            else:
                for edge_idx in range(bounds_minmax[0],bounds_minmax[1]+1,1):
                    left_distance = self.edges[edge_idx]-self.centers[center_idx]
                    right_distance =  self.edges[edge_idx+1]-self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    #if np.abs(right_distance) >= self.r_particle: 
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)

    def sphere_sphere_intersection(self, r, R, d):
        """
        sphere_sphere_intersction computes the volume of intersection of two spheres. The sphere with redius R
        is at the origin (0,0,0) while the other one with radius r is located along axis x at x=d (d,0,0). This function can be used to find the local volume fraction of a spherical beads in the radial direction of in a space with spherical
        symmetry.

        Reference: https://mathworld.wolfram.com/Sphere-SphereIntersection.html

        Inputs:
        r: the radius of the sphere locared along axis x.
        R: the radius of the sphere located at the origin.
        d: the distance of the the off-origin sphere from the origin along axis x.

        Returns:
        V: volume of intersection.

        Requirements:
        numpy package for constant Pi.
        """
        # By define Rmax and Rmin, we handlthe situtations in which d = 0:
        Rmax = max(R,r)
        Rmin = min(R,r)
        if r ==0 or R == 0:
            V = 0 
        else : 
            if d == 0: # the small sphere resides completely in the large one.
                V = 4*np.pi*Rmin**3 / 3 

            elif d >= Rmin+Rmax: # the spheres are either tengential to eachother or not interesting.
                V = 0
            else:
                if d <= Rmax-Rmin: # the small sphere resides completely in the large one.
                    V = 4*np.pi*Rmin**3 / 3

                else :
                    V = np.pi * (Rmax + Rmin - d)**2 * (d**2 + 2*d*Rmin - 3*Rmin**2 + 2*d*Rmax + 6*Rmin*Rmax - 3*Rmax**2) / (12*d)
        return V
    
    def _concentric_bounds(self):
        """
        _find_concentric_bounds finds the lowest and highest bin edges by which a spherical bead residing at the center of a bin intersects. There are three types of bins: concenteric bins (and thus concentric edges) along the radial direction, consecutive bins along the longitudinal direction, and cyclic consecutive bins along the polar and azimuthal directions. This method correctly finds the bounds in the first two cases. This method find the bounds in a concentric-bins shceme along the  ridial direction.
        
        Cuation: 
        This is a 1D problems since the bins, bin centers and bin edges are along the same axis. In reality, the center of a spherical bead can be anywhere in a bin, but it is placed at the center of the bin along the axis of interest.
        The center of a particle can be in any of the bins, so 
        While this assumption makes the value of the local volume fraction slightly inaccurate, it singnificantly reduce the computational cost of measure the local volume fraction in a system with a large number of particles.
        The positions and r_concentrics have some relations:
        1. len(edges) = len(centers) + 1
        2. There is only one center between two consecutive edges.
        3. All the arrays are increasing funcions of their indices.
        
        Instead of the leftmost/rightmost pairs (which are more appropriate for the longitudinal direction), the innermost/outermost is used.
        """
        innermost = self.centers - self.r_particle # The minimum distance of the r_atoms perimeter from the origin
        innermost_idx = np.zeros(len(innermost),dtype=int) # Initiate the leftmost bound with the lowest possible bound
        outermost = self.centers + self.r_particle # The maximum distance of the r_atoms perimeter from the origin
        outermost_idx = (len(outermost)-1) * np.ones(len(outermost),dtype=int) # Initiate the rigtmost bound with the highest possible bound
        for idx, innermost_value in enumerate(innermost):
            for edge_idx in range(len(self.edges[:-1])):
                if (innermost_value >= self.edges[edge_idx]) and (innermost_value < self.edges[edge_idx+1]): # the inner edge index of the bin is set as the index of the bin by which the innermost side of the bead intersects.
                    innermost_idx[idx] = edge_idx + 1 # For the innermost bond, the intersection of the the bead and the outer edge of the innermost bin is important!
                if (outermost[idx] >= self.edges[edge_idx]) and (outermost[idx] < self.edges[edge_idx+1]): # the outer edge index of the bin is set as the index of the bin by which the outermost side of the bead intersects.
                    outermost_idx[idx] = edge_idx # For the outermost bond, the intersection of the the bead and the inner edge of the outermost bin is important!
        self.particle_bounds = np.column_stack((innermost_idx,outermost_idx)) 

    def _concentric_vol_shares(self):
        """
        _concentric_vol_shares computes the portion of the volume of a bead (or a sphere) in different spherical or cylindrical concenteric bins. The intersection_calculator computes the volume of intersection between the bead and the bin (sphere or cylinder) of radius r_concentric located at origin.
        For the sake of simpilicity, the sphere-sphere interesction is also used for the sphere-cylinder intersection in the radial direction in the cylindrical and slit goemteries. Hence, the concentric_vol_shares can be used in all the three different radial directions.

        Caution:
        The volume share of the innermost bin is equal to the volume of the intersection given by the intersection_calculator.
        The volume share of the other bins is equal to the volume of the intersection of the bin with the inner edge of radius r_concentric substracted by the volume of the intersection of the previous bin with the inner edge given by its radius r_concentric.

        This function can be used in calculation of the local volume fraction of beads in two different situation:
        1. The radial direction in the cubic geometry with the sphere_sphere_intersection function as the intersection_calculator.
        2. The radial dirction in the cylindrical geometry with the sphere_cylinder_intersection function as the intersection_calculator.
        
        The centers and edges have these relations:
        1. len(edges) = len(centers) + 1
        2. There is only one center between two consecutive edges.
        3. All the arrays are increasing funcions of their indices.

        Parameters:
        intersection_calculator: the function computes the volume of intersection between the bead with r_particle as its radius centerred at centers  and the concentric spherical or cylindrical shells with edges as their radii all centered at the origin. It is self.sphere_sphere_intersection for this method.

        """
        self.volume_shares = {}
        for center_idx, bound_minxax in enumerate(self.particle_bounds):
            self.volume_shares[center_idx] = {}
            intersect_vol_previous = 0 # The volume share of the lowest bin all comes from itself.
            #for edge_idx in range(bound_minxax[0]+1,bound_minxax[1]+1,1):
            for edge_idx in range(bound_minxax[0],bound_minxax[1]+2,1): # The index of upper limit is increased by 2 units since 1 units becuase of the range function and the other one is because of the share of the last outmost bin. The share of outmost bin is the volume of sphere minus the volume of the interestion of the bead with the laregest bin edge samller than the outmost edge of the bead.
                intersect_vol = self.sphere_sphere_intersection(self.r_particle, self.edges[edge_idx], self.centers[center_idx])
                self.volume_shares[center_idx][edge_idx-1]= intersect_vol - intersect_vol_previous # The intersection volume between bead and edge i belongs to bin i-1. The volume share of the previous bin i-1 should be subsracted from bin i; draw a figure to realize it!.
                intersect_vol_previous = intersect_vol # set this intersection_vol as the old one.

    def _vol_shares_type(self):
        """
        _vol_shares_type chooses how the volume_shares should be measured based on the given direction. Currently, the vol_shares method is implemented in the radial direction in all the geometries (see the notes for _concentric_vol_shares) and the longitudinal direction in the cylindrical geometry.
        """
        if self.direction == "radial":
            self._concentric_bounds()
            self._concentric_vol_shares()
        elif self.direction == "longitudinal":
            self._consecutive_bounds()
            self._consecutive_vol_shares()
        else:
            raise ValueError(f"'volume_shares' is not defined in the {self.direction} direction.")

    def _set_args(self, col_name):
        """
        _set_args set the arguments for the integrads along different directions in different geometries.
        
        Parameters:
        col_name: the name of column for which the arguments are set.
        """
        chosen_properties = self.properties[self.properties.filename==col_name]
        self._args = {
            "cubic": {
                "radial" : (1,), # cocentric spherical shells; constant is redundant and merely defined to make the use of args parameter of scipi.integral.quad function consistent among integrands.
                "polar": (0.5 * chosen_properties['lcyl'].values[0], ), # in a cubic or free space, the radius of the space is half of the length of simulation box
                 "azimuthal": (0.5 * chosen_properties['lcyl'].values[0], ),
            },
            "slit": {
                "radial" : (chosen_properties['lcyl'].values[0], ),
                "polar" : (chosen_properties['lcyl'].values[0], chosen_properties['dcyl'].values[0], ),
                "longitudinal": (chosen_properties['dcyl'].values[0], )
            },
            "cylindrical": {
                "radial": (chosen_properties['lcyl'].values[0], ),
                "polar": (chosen_properties['lcyl'].values[0], chosen_properties['dcyl'].values[0], ),
                "longitudinal": (chosen_properties['dcyl'].values[0], )
            }
        }
    
    def _number_density(self, col_name):
        """
        _number_density calculate the local number density along the given direction in the given geometry. The local number density is normalized to give the area under the curve equal to one. The number density in each simulation is an average over the number densities collected every X time steps, so there are N=L/X measurements of the local number desnity in each simulation where L is total number of time steps in the simulation. For the cylindrical sum rule project X=5000 and the total number of time steps is 7*10^7, so N=14001. For the ensemble-averages, each local number desnity is averages over the M simulations in an ensemble. FOr the cylindrical sum rule project, M=8.
        
        Parameters:
        col_name: the name of column for which the number density is calculated.
        """
        integrand = self._integrands[self.geometry][self.direction]
        arguments = self._args[self.geometry][self.direction]
        bin_vols = np.array([integrate.quad(integrand, self.edges[idx] ,self.edges[idx]+self.bin_size, args=arguments)[0] for idx in range(len(self.edges[:-1]))])
        #self.rho[col_name] = np.divide(self.histogram[col_name], bin_vols)
        self.rho[col_name] = self.histogram[col_name].divide(bin_vols) # elf.histogram[col_name] and bin_vols have the same size.
        # the sum of rho is not equal to the bulk number density (r=infiity) natom/cell_vo. This arises from the way we descritize the local number desnity.
        if self.normalized:
            self.rho[col_name] = self.rho[col_name] / self.rho[col_name].sum() # time averaging: the sum of histograms = natoms * nframes. normalization: the sum of the number density is now 1.

    def _volume_fraction(self, col_name):
        """
        _volume_fraction computes the local volume fraction along the direction of interest in the given goemtetry. All the particles have the same shape. The local volume fraction is normalized to give the integral of p(r or theta or z)=1 along the direction of interest in the region of interest. See the explnation for the _number_density method and Distributions class.
        
        Parameters:
        col_name: the name of column for which the volume fraction is calculated.
        """
        rho = self.rho[col_name].to_numpy()
        n_centers = len(rho)
        phi = np.zeros(n_centers)
        for center_idx in range(n_centers):
            for vol_share_idx, vol_share in self.volume_shares[center_idx].items():
                #phi[vol_share_idx] = phi[vol_share_idx] + (rho[center_idx] * (vol_share)) # this sum and the brelow one are equivalent.
                #phi[center_idx] = phi[center_idx] + (rho[vol_share_idx] * (vol_share))
                phi[center_idx] = phi[center_idx] + (rho[vol_share_idx] * (vol_share))
        # the sum of phi is not equal to the bulk volume fraction (r=infiity) natom*vol_per_atom/cell_vol. this arises from the way we descritize the local volume fraction and the way we relate it to the local number density.
        if self.normalized:
            self.phi[col_name] = phi / np.sum(phi)
        else:
            self.phi[col_name] = phi

    def _run(self):
        """
        _run perform a list of operation over the columns in the input histogram.
        """
        for col_name in self.histogram.columns:
            self._set_args(col_name)
            self._number_density(col_name)
            self._volume_fraction(col_name)      