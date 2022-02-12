 def sphere_sphere_intersection(self, r, R, d):
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
        self.volume_shares = {}
        for center_idx, bound_minxax in enumerate(self.particle_bounds):
            self.volume_shares[center_idx] = {}
            intersect_vol_previous = 0 # The volume share of the lowest bin all comes from itself.
            #for edge_idx in range(bound_minxax[0]+1,bound_minxax[1]+1,1):
            for edge_idx in range(bound_minxax[0],bound_minxax[1]+2,1): # The index of upper limit is increased by 2 units since 1 units becuase of the range function and the other one is because of the share of the last outmost bin. The share of outmost bin is the volume of sphere minus the volume of the interestion of the bead with the laregest bin edge samller than the outmost edge of the bead.
                intersect_vol = self.sphere_sphere_intersection(self.r_particle, self.edges[edge_idx], self.centers[center_idx])
                self.volume_shares[center_idx][edge_idx-1]= intersect_vol - intersect_vol_previous # The intersection volume between bead and edge i belongs to bin i-1. The volume share of the previous bin i-1 should be subsracted from bin i; draw a figure to realize it!.
                intersect_vol_previous = intersect_vol # set this intersection_vol as the old one.