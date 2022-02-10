import datetime
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 
from MDAnalysis.coordinates.LAMMPS import DumpReader
from MDAnalysis import transformations as mdatransform
from MDAnalysis.analysis import (diffusionmap, align, rms)

from ..polyphys import CellAttributes

def log_outputs(log_file, geometry):
    cell_attrs = CellAttributes(log_file,geometry)
    output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}'
    # if dangerous build needed 
    details_out = output + "-log_details.csv"
    with open(details_out, "w") as detailsfile:
        # neig_modify delay NUM every NUM check YES/NO:
        detailsfile.write('groupname,filename,ens,run_seg,rskin,delay,every,check,epsilon,dcrowd,')
        detailsfile.write('ncrowd,lcyl,dcyl,nmon,total_time_s,cores,timestep,atoms,ts_per_sec,')                                           
        # Section columns: min time, avg time, max time, %varavg, %total"
        # Section rows: Pair, Bond, Neigh, Comm, Output, Modify, Other
        detailsfile.write('pair_avg_s,pair_pct,bond_avg_s,bond_pct,neigh_avg_s,neigh_pct,comm_avg_s,')
        detailsfile.write('comm_pct,output_avg_s,output_pct,modify_avg_s,modify_pct,other_avg_s,other_pct,dangerous\n')
    
    runtime_out = output + "-log_runtime.csv"
    with open(runtime_out, "w") as runfile:
        runfile.write('groupname,filename,ncores,natoms,wall_time\n')
      
    return details_out , runtime_out

def lammps_log_details(log_files, details_out , runtime_out):
    for file in log_files:
        filename = file[0].split('.log')
        filename = filename[0]
        filename = filename.split('/')[-1]
        ens = filename.split('ens')[-1]
        groupname = filename.split('ens')[0]
        with open(file[0],'r') as log,\
        open(details_out, "a") as detailsfile,\
        open(runtime_out, "a") as runfile:   
            line = log.readline()
            # The other of while loop are important
            #neigh_modify delay NUM every NUM check YES/NO page NUM one NUM: 
            j = 1
            while line:
                while line.startswith('variable'):
                    words = line.split()
                    line = log.readline()
                    if words[1] == 'epsilon1':
                        epsilon = words[3]
                    if words[1] == 'sig2':
                        dcrowd = words[3]  
                    if words[1] == 'n_crowd':
                        ncrowd = words[3] 
                    if words[1] == 'lz':
                        lcyl = str(2*float(words[3]))  
                    if words[1] == 'r':
                        dcyl = str(2*float(words[3]))
                    if words[1] == 'n_bug':
                        nmon = words[3]
                if line.startswith('neighbor'):
                    words = line.split()
                    rskin = words[1].strip() # rskin
                #neigh_modify delay NUM every NUM check YES/NO page NUM one NUM:  
                if line.startswith('neigh_modify'):
                    words = line.split()
                    # picking the NUMs and Yes/No from neigh_modify command
                    delay = words[2].strip()
                    every = words[4].strip()
                    check = words[6].strip()
                if line.startswith('Loop time'):
                    detailsfile.write(groupname)
                    detailsfile.write(",")
                    detailsfile.write(filename)
                    detailsfile.write(",")
                    detailsfile.write(ens)
                    detailsfile.write(",")
                    detailsfile.write(str(j))#total time
                    detailsfile.write(",")
                    j += 1
                    # neighbor and neigh_modify ocurres  one time but other occures 15 times.
                    detailsfile.write(rskin) # rskin
                    detailsfile.write(",")
                    detailsfile.write(delay) # delay
                    detailsfile.write(",")
                    detailsfile.write(every) # every
                    detailsfile.write(",")
                    detailsfile.write(check) # check
                    detailsfile.write(",")
                    detailsfile.write(epsilon) # epsilon
                    detailsfile.write(",")
                    detailsfile.write(dcrowd) # dcrowd
                    detailsfile.write(",")
                    detailsfile.write(ncrowd) # ncrowd
                    detailsfile.write(",")
                    detailsfile.write(lcyl) # lcyl
                    detailsfile.write(",")
                    detailsfile.write(dcyl) # dcyl
                    detailsfile.write(",")
                    detailsfile.write(nmon) # nmon
                    detailsfile.write(",")
                    words = line.split()
                    detailsfile.write(words[3].strip())#total time
                    detailsfile.write(",")
                    ncores = words[5].strip()
                    detailsfile.write(ncores)# # of cores
                    detailsfile.write(",")
                    detailsfile.write(words[8].strip())# total timesteps
                    detailsfile.write(",")
                    natoms = words[11].strip()
                    detailsfile.write(natoms)# total atoms
                    detailsfile.write(",")                          
                if line.startswith('Performance:'):
                    words = line.split()
                    detailsfile.write(words[3].strip())# timesteps per second
                    detailsfile.write(",")
                if line.startswith('Section'):
                    _ = log.readline()
                    for i in range(6): # Section rows: Pair, Bond, Neigh, Comm, Output, Modify, Other
                        # Section columns: min time, avg time, max time, %varavg, %total"
                        line = log.readline()
                        sect_min = line.split('|')[2].strip()
                        detailsfile.write(sect_min)
                        detailsfile.write(",")

                        sect_pct = line.split()[-1] # Pair pct of total time
                        detailsfile.write(sect_pct)
                        detailsfile.write(",")
                    line = log.readline()
                    sect_min = line.split('|')[2].strip()
                    detailsfile.write(sect_min)
                    detailsfile.write(",")
                    sect_pct = line.split()[-1] # Pair pct of total time
                    detailsfile.write(sect_pct)
                    detailsfile.write(",")
                if line.startswith('Dangerous'):
                    words = line.split()
                    detailsfile.write(str(int(words[-1]))) # # number of dangerous builds
                    #detailsfile.write(",")
                    detailsfile.write("\n")
                # runtime files
                if line.startswith('Total wall time'):
                    runfile.write(groupname)
                    runfile.write(",")
                    runfile.write(filename)
                    runfile.write(",")
                    runfile.write(ncores)
                    runfile.write(",")
                    runfile.write(natoms)
                    runfile.write(",")
                    words = line.split()
                    runfile.write(words[-1]) # total wall time
                    runfile.write("\n")
                line = log.readline()

# Handling I/O:
def cylinder_write(output,cellAttrs):
    """
    This function writes the information stored in a cellAttrs object into a file.
    Parameters: 
    cellAttrs (an object from cellAtttribute class):this object contains information about the parameters used as inout in a simulation.
    ensemble (python file object in "a" mode): the output file that contains the statistics about a simulation.
    
    Returns: ---
    
    Requirements:
    CellAttributes class
    
    """
    output.write("{},{},{},{},{},".format(cellAttrs.filename,cellAttrs.dmon,cellAttrs.nmon,cellAttrs.dcyl,cellAttrs.lcyl))
    output.write("{},{},{},{},{},".format(cellAttrs.phi_m_bulk,cellAttrs.rho_m_bulk,cellAttrs.eps_mon_wall,cellAttrs.ncrowd,cellAttrs.dcrowd))
    output.write("{},{},{},".format(cellAttrs.phi_c_bulk,cellAttrs.rho_c_bulk,cellAttrs.ens))



# Direct analysis on trajectory files:
def chain_stats(array, output):
    """
    chain_stats computes the mean, standard deviation, variance, and standard error of the mean for a vector array and write the results to an output file.
    
    Caution: 
    anom stands for analysis of the mean
    this function has not been finished yet; other statistical analysis should be added.
    this function can be combine with anova (analysis of varaince) and create a class.

    Parameters:
    array (a numpy ndarray): a vector of the same lenght as the total number of frames (snapshots) in the trajectory.
    output (python file object in "a" mode): the output file that contains the statistics about this vector.
    
    Return:
    No retrun, but write some information to file.
    --
    
    Requirements:
    numpy
    """
    mean = np.mean(array) # mean
    std = np.std(array,ddof=1)# The bias-corrected sample deviation: N-ddof instead of N with ddof=1
    var = np.var(array,ddof=1)# The bias-corrected sample variance: N-ddof instead of N with ddof=1
    sem = std / np.sqrt(len(array)) # Standard error of mean
    output.write("{},{},{},{},".format(mean, std, var, sem))

def error_calc_block(data, filename):    
    """
    This function computes the statistical inefficiency (si) and errors associated with a physical quantity via simulations data.
    
    Using Flybbjerg-Peterson block average method, the plateau should be evident after 6-8 transformations.

    This function is written based on the related python scripts in the following github reppsitories:
    https://github.com/Allen-Tildesley/examples/blob/master/python_examples/error_calc.py
    https://github.com/MDAnalysis/MDAnalysisCookbook/blob/master/examples/blocks.py

    Both of these code snippets are based on: "Error estimates on averages of correlated data",HG Flyvbjerg and
    H Petersen. J Chem Phys, 91, 461 (1989).
    https://doi.org/10.1063/1.457480

    data: the array of calculated physical quantity in numpy narry format.
    save_to: wherther write the graphs to file or not.
    filename: the name of output file
    
    error_calc uses a copy of the input data.
    if error_calc finds matplotlib, it also draws a diagram to estimate se.

    Parameters
    r (natoms*ndim numpy ndarray):the positions of the N atoms within an atom group in a frame. This positions are arreanged by atom number
    form 1 to N.
    filename (string): the name of the simulations for which the block analyzing is done.
    
    Return:
    a plot of the statistical inefficiency vs the number of block transformations. 
    
    Requirements:
    matplotlib, math, and numpy.
    """
        
    nframe = len(data) # size of the data

    data_mean = np.mean(data)       # sample average (mean)

    data_var = np.var(data, ddof=1)  # Bias-corrected sample variance
    data_var_err = np.sqrt(2/(nframe-1)) * data_var  # error in Bias-corrected sample variance

    data_sem = np.sqrt(data_var / nframe) # standard error of the mean (SEM) neglecting any correlations
    data_sem_err = np.sqrt(1/(2*(nframe-1)))* data_sem # error in SEM

    bdata = data.copy()

    cblock = np.zeros(0,dtype=np.int8) # Count the number of block tranformartion used.
    cblock = np.append(cblock,0) # initializing with 0

    nblock = np.zeros(0,dtype=np.int8) # Number of block
    nblock = np.append(nblock,nframe) # # initializing with nframe

    tblock = np.zeros(0,dtype=np.int8) # Size of each block
    tblock = np.append(tblock,1) # initializing with 1

    bvar = np.zeros(0) # variance
    bvar = np.append(bvar,data_var) # initializing with the data_var

    bvar_err = np.zeros(0) # error in sem
    bvar_err = np.append(bvar_err,data_var_err) # initializing with data_var_error

    bsem = np.zeros(0) # variance of the mean or the standrad error of mean (sem)
    bsem = np.append(bsem,data_sem) # initializing with data_sem

    bsem_err = np.zeros(0) # error in sem
    bsem_err = np.append(bsem_err,data_sem_err) # initializing with data_sem_error

    si_initial = tblock[-1] * bvar[-1] / data_var # inintial statistical inefficiency (SI)
    si_initial_err = np.sqrt(2/(nframe-1)) * si_initial # error in inintial SI

    si = np.zeros(0) # statistical inefficiency (SI)
    si = np.append(si,si_initial) # initializing with initial si which is 1!

    si_err = np.zeros(0) # error in SI
    si_err = np.append(si_err,si_initial_err) # initializing with si_err which is 0 since si has the trivial value 1.

    outerror = filename+'_error_analyze.csv'
    
    with open(outerror, mode="w") as errorfile:
        #errorfile.write("# Blocking method for estimating statistical inefficiency:\n")
        errorfile.write("{:>25},{:>25},{:>25},{:>15},{:>25},{:>25},{:>25},{:>25},{:>25}\n".format('cblock','tblock', 'nblock', 'Var','Var_err', 'SEM','SEM_err', 'SI','SI_err') )
        errorfile.write("{:25d},{:25d},{:25d},{:25.15f},{:25.15f},{:25.15f},{:25.15f},{:25.20f},{:25.20f}\n".format(cblock[-1],tblock[-1], nblock[-1], bvar[-1], bvar_err[-1], bsem[-1], bsem_err[-1], si[-1], si_err[-1]))
        while True: # Loop over number, and hence length, of blocks
            nblock = np.append(nblock, nblock[-1] // 2) # Halve the number of blocks, rounding down if nblock is odd
            tblock = np.append(tblock, tblock[-1] * 2)    # Double the block length

            bdata[0:nblock[-1]] = ( bdata[0:2*nblock[-1]-1:2] + bdata[1:2*nblock[-1]:2] ) / 2.0 # Blocking transformation, halving the data set
            bvar = np.append(bvar, np.var(bdata[0:nblock[-1]], ddof=1))  # Bias-corrected variance of block averages
            bvar_err = np.append(bvar_err, np.sqrt(2/(nblock[-1]-1))*bvar[-1] )
            bsem = np.append(bsem, np.sqrt( bvar[-1] / nblock[-1] ) )   # Estimate of SEM
            bsem_err = np.append(bsem_err, np.sqrt((1/(2*(nblock[-1]-1))))*bsem[-1] )  # error in SEM
            si = np.append(si, tblock[-1] * bvar[-1] / data_var) # Statistical inefficiency (SI)
            si_err = np.append(si_err, np.sqrt((1/(2*(nblock[-1]-1))))*si[-1] ) # error in SI
            errorfile.write( "{:25d},{:25d},{:25d},{:2515f},{:25.15f},{:25.15f},{:25.15f},{:25.20f},{:25.20f}\n".format(cblock[-1], tblock[-1], nblock[-1], bvar[-1], bvar_err[-1], bsem[-1], bsem_err[-1], si[-1], si_err[-1]))
            cblock = np.append(cblock,cblock[-1]+1)
            if nblock[-1] <= 3: # loop counter
                break
        #errorfile.write('# First row is the original data set statistics. This row should not be used to estimate si by fitting\n')
        #errorfile.write('# Plateau at large tblock (small nblock).\n')
        #errorfile.write('# Plateau should be evident after 6-8 transformation.\n')
        #errorfile.write('# This method should agree quite well with exact error estimate.\n')
        #errorfile.write('# Plot SI or error**2 (or SEM**2) against 1/tblock or log2(tblock) to find the plateau .\n')

    # si vs nblock
    fig, axes = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
    axes.grid(True, which="both")
    axes.errorbar(cblock, si,yerr=si_err, fmt='--o')
    axes.set_xlabel(r"Number of transformation: $n_{block}}$")
    axes.set_ylabel(r"Statistical inefficiency: $s(n_{block})$")
    picname = filename+'_plot_err.pdf'
    plt.savefig(picname,dpi=400)
    plt.close() # not showing the figure in jupyter.

def end_to_end(r):
    """
    calculates the end-to-end distance of a polymer.
    
    Caution:
    Since all the atoms composing the polymer have the same mass, the masses of them are not considered in
    the following lines. If the masses are diffrent, you need to use average method from the numPy.

    Parameters:
    r (natoms*ndim numpy ndarray):the positions of the N atoms within an atom group in a frame. This positions are arreanged by atom number
    form 1 to N.

    Return:
    the end-to-end distance: numpy.float
    
    Requirements:
    MDAnalysis and numpy. This functions does not use MDAnalysis directly, but it use the fact that the vector r is ordered by the
    monomoer id from one end of the chain to the other end of that.
    """
    
    com = np.mean(r, axis=0) # center of mass of the chain.
    r = r - com # positions of monomers in the com reference frame.
    r = r[-1] - r[0] # end-to-end vector; the position vector r is arranged based on the id of the monomers (see Caution above).
    return np.linalg.norm(r)


def max_distance(r):
    import numpy as np
    """
    max_distance returns the maximum ditance in each of the three Cartesian direction. The minimum adn maximum valaues by which the maximum distance is computed are not for the same particle but are in the same frame or snapshot or time step.

    Caution:
    Since all the atoms composing the polymer have the same mass, the masses of them are not considered in
    the following lines. If the masses are diffrent, you need to use average method from the numPy.

    Parameters
    r (natoms*ndim numpy ndarray):the positions of the N atoms within an atom group in a frame. This positions are arreanged by atom number
    form 1 to N.
    
    Return:
    xmax, ymax, zmax in each timestep (snapshot): numpy.float
    
    Requirements:
    MDAnalysis and numpy. This functions does not use MDAnalysis directly, but it use the fact that the vector r is ordered by the
    monomoer id from one end of the chain to the other end of that.
    """
    print("The max/min distance measured in the COM frame of reference. Is this correct?")
    com = np.mean(r, axis=0) # center of mass of the chain.
    r = r - com # positions of monomers in the com reference frame.
    xmax = np.abs(np.amax(r[:,0])-np.amin(r[:,0]))
    ymax = np.abs(np.amax(r[:,1])-np.amin(r[:,1]))
    zmax = np.abs(np.amax(r[:,2])-np.amin(r[:,2]))
    return [xmax, ymax, zmax]

def fsd(coordinates, axis=2):
    """
    The function calculate the average size/diameter of a polymer confined in a cylindrical geometry based on
    the farthermost distance concept in one frame.
    fsd stands for Feret's statistical diameter: other names are the mean span dimension,
    the farthermost distance, or the mean caliper diameter.
    See: "A Theoretical Study of the Separation Principle in Size Exclusion Chromatography",Wang Y Teraoka
    I Hansen FY Peters GH Ole H. Macromolecules 2010, 43, 3, 1651-1659
    https://doi.org/10.1021/ma902377g

    coordinates is the positions of the atoms within an atom group in a frame.
    coordinates is a natoms*ndim numpy array.
    axis is the index of the axis of the cylinder; by default it is the z axis.

    this function need numpy package.

    Since all the atoms composing the polymer have the same mass, the masses of them are not considered in
    the following lines. If the masses are diffrent, you need to use average method from the numPy.
    """
    com = np.mean(coordinates, axis=0) # center of mass of the chain.
    coordinates = coordinates - com # positions of monomers in the com reference frame.
    return np.ptp(coordinates[:,axis])

def bin_ceate(sim_name, bin_size, lmin, lmax, edge_name, save_to):
    """
    bin_ceate produces arrays of bins and histograms

    inputs:
    sim_name: the name of the simulation.
    bin_size (float): size of each bin.
    lmin (float): lower bound of the system in the direction of interest.
    lmax (float): upper bound of the system in the direction of interest.
    edge_name: name of the variable for which the histogram is computed.
    save_to (str): address to which the output is saved
    
    returns:
    save a txt files with name sim_name+'_'+edge_name+'.txt'
    bin_edges: an array of bins' edges.
    histogram: an int-tyoe array of zero values collecting the number of occurance in each bin in
    all the time-steps of interest.

    Requirements:
    numpy
    """
    bin_edges = np.arange(lmin,lmax+bin_size,bin_size)
    hist_collectors = np.zeros(len(bin_edges)-1,dtype=np.int16)
    np.savetxt(save_to+sim_name+'-'+edge_name+'.csv', bin_edges, delimiter=',')
    return bin_edges, hist_collectors

class AmirDumpReader(DumpReader):
    
    format = 'AMIRLAMMPSDUMP'
    
    def __init__(self, filename, fractional=False, **kwargs):
        self.fractional = fractional
        super().__init__(filename, **kwargs)
        
    def _read_next_timestep(self):
        f = self._file
        ts = self.ts
        ts.frame += 1
        if ts.frame >= len(self):
            raise EOFError

        f.readline() # ITEM TIMESTEP
        step_num = int(f.readline())
        ts.data['step'] = step_num

        f.readline() # ITEM NUMBER OF ATOMS
        n_atoms = int(f.readline())
        if n_atoms != self.n_atoms:
            raise ValueError("Number of atoms in trajectory changed "
                             "this is not suported in MDAnalysis")

        triclinic = len(f.readline().split()) == 9  # ITEM BOX BOUNDS
        if triclinic:
            xlo, xhi, xy = map(float, f.readline().split())
            ylo, yhi, xz = map(float, f.readline().split())
            zlo, zhi, yz = map(float, f.readline().split())

            box = np.zeros((3, 3), dtype=np.float64)
            box[0] = xhi - xlo, 0.0, 0.0
            box[1] = xy, yhi - ylo, 0.0
            box[2] = xz, yz, zhi - zlo

            xlen, ylen, zlen, alpha, beta, gamma = mdamath.triclinic_box(*box)
        else:
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())
            xlen = xhi - xlo
            ylen = yhi - ylo
            zlen = zhi - zlo
            alpha = beta = gamma = 90.
        ts.dimensions = xlen, ylen, zlen, alpha, beta, gamma

        indices = np.zeros(self.n_atoms, dtype=int)

        f.readline()  # ITEM ATOMS etc
        for i in range(self.n_atoms):
            idx, _, xs, ys, zs = f.readline().split()

            indices[i] = idx
            ts.positions[i] = xs, ys, zs

        order = np.argsort(indices)
        ts.positions = ts.positions[order]
        # by default coordinates are given in scaled format, undo that
        if self.fractional:
            ts.positions = distances.transform_StoR(ts.positions, ts.dimensions)

        return ts



def extract_trj_bug(simulation_pair, geometry, save_to="./"):
    """
    extract_trj_bug does all the analysis on the a whole trajectory (simulation),

    Parameters:
    simulation_pair (a tuple of size 2): the pair of topology(first argument) and trajectory (second argument) of a simulation.
    geometry (string): the name of the geometry of the simulation box.
    save_to (str): address to which the output is saved

    Requirements:
    All the above-defined fucntions and classes, MDAnalysis.
    
    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.
    For other goemteries you need to redfined them.
    """  
    print("This script for cylindrical geomery")
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = CellAttributes(simulation_pair[1],geometry)
    sim_name = cellAttrs.filename
    #ensGroup = "N"+str(cellAttrs.nmon)+'D'+str(cellAttrs.dcyl)+"_"
    outfile = save_to+sim_name+"-properties.csv"

    print("")
    print(sim_name+" is analyzing...")
    print("")
    with open(outfile, mode="w") as ensemble:
        ensemble.write(cellAttrs.cols)
        cylinder_write(ensemble,cellAttrs)

    print("Analyzing...")
    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(simulation_pair[0], simulation_pair[1], format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # radial direction of the cylindrical coordinate system
    edge_name = 'rEdges'
    bin_size = 0.1 * min(cellAttrs.dmon,cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.1*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    redges, rhist_collectors = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name, save_to)

    # z direction of the cylindrical coordinate system
    edge_name = 'zEdges'
    bin_size = 0.5 * min(cellAttrs.dmon,cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.5*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = cellAttrs.lcyl / 2.
    zedges, zhist_collectors = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name, save_to)

    # theta of the cylindrical coordinate system
    edge_name = 'thetaEdges'
    bin_size = np.degrees(np.pi/36) #bin size 5 degrees
    print(f"Bin size for theta in the cylindrical geometry is set to {bin_size} in degrees.")
    thetaedges, thetahist_collectors = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degres(np.pi), edge_name, save_to)

    # Chain end-to-end size distribution
    edge_name = 'rFloryEdges'
    lmax = cellAttrs.nmon * cellAttrs.dmon # The contour size of the chain
    bin_size = 0.5 * cellAttrs.dmon # in a_c=1.0 unit; for a_c<a_mon_small; check this
    print(f"Bin size for the PDF of end-to-end vector is set to 0.5*cellAttrs.dmon={bin_size} in a_m=1 units.")
    rFloryedges, rFloryhist_collectors = bin_ceate(sim_name, bin_size, 0, lmax, edge_name, save_to)

    #instantaneous quantity
    fsd_t = np.empty([0]) # the furthermost size using my own function fsd_cylinder.
    rFlory_t = np.empty([0]) # the end-to-end size, using my own function end_to_end.
    gyr_t = np.empty([0]) # radius of gyration

    if  any([rhist_collectors.any() != 0, zhist_collectors.any() != 0, thetahist_collectors.any() != 0, rFloryhist_collectors.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #chain size : the furthermost distance:
        fsd_t = np.append(fsd_t, np.array([fsd(chrm.positions)]), axis = 0)

        # radius of gyration
        gyr = chrm.radius_of_gyration()
        gyr_t = np.append(gyr_t, np.array([gyr]), axis = 0) # radius of gyration in each time step

        # end-to-end size distribution in eahc frame.
        rms = end_to_end(chrm.positions)
        rFlory_t = np.append(rFlory_t, np.array([rms]), axis = 0) # end-to-end vectors in each timestep

        # histogram in r direction
        dummy_hist, _ = np.histogram(rms, rFloryedges) # For one timestep, rms is equal to norm of end-to-end, so we use the later for histogram
        rFloryhist_collectors = np.add(rFloryhist_collectors, dummy_hist) 

        #number density in the cell's frame of reference
        # histogram in r direction
        rmon = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rmon, redges)
        rhist_collectors = np.add(rhist_collectors, dummy_hist)

        # histogram in z direction
        zmon = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zmon, zedges)
        zhist_collectors = np.add(zhist_collectors, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahist_collectors = np.add(thetahist_collectors, dummy_hist)

    np.savetxt(save_to+sim_name+'-rHists.csv', rhist_collectors, delimiter = ',')
    np.savetxt(save_to+sim_name+'-zHists.csv', zhist_collectors, delimiter = ',')
    np.savetxt(save_to+sim_name+'-thetaHists.csv', thetahist_collectors, delimiter = ',')
    np.savetxt(save_to+sim_name+'-rFloryHists.csv', rFloryhist_collectors, delimiter = ',')
    np.savetxt(save_to+sim_name+'-fsd_t.csv', fsd_t, delimiter = ',')
    np.savetxt(save_to+sim_name+'-rFlory_t.csv', rFlory_t, delimiter = ',')
    np.savetxt(save_to+sim_name+'-gyr_t.csv', gyr_t, delimiter = ',')

    with open(outfile, mode="a") as ensemble:
        chain_stats(fsd_t, ensemble) # fsd,fsd_std,fsd_var,fsd_sem,
        chain_stats(rFlory_t, ensemble) # rflory,rflory_std,rflory_var,rflory_sem,
        chain_stats(gyr_t, ensemble) # gyr,gyr_std,gyr_var,gyr_sem,
        ensemble.write("{}\n".format(cell.trajectory.n_frames))
    print('done.') 
    
def rmsd_trj_bug(fname, geometry, save_to):
    """
    This function comptues the rmsd of a bug's trajectory.
    
    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.
    For other goemteries you need to redfined them.
    
    Parameters:
    fname (string): the name of trajectory or simulation.
    geometry (string): the name of the geometry of the simulation box.
    save_to (str): address to which the output is saved

    Requirements:
    All the above-defined fucntions and classes, MDAnalysis.
    """
    print("Setting the name of analyze file...\n")
    print("Doing RMSD analysis...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = CellAttributes(fname[0], geometry = geometry)
    
    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit
    
    sim_name = cellAttrs.filename
    cell = mda.Universe(fname[1], fname[0], format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    cell.transfer_to_memory(step=50,verbose=False)
    
    aligncell = align.AlignTraj(cell, cell, select = 'resid 1', filename = sim_name + '.dcd').run()
    matrix = diffusionmap.DistanceMatrix(cell, select = 'resid 1').run()
    np.save(save_to+sim_name+'-rmsdMatrix.npy',matrix.dist_matrix) 
                
    
def extract_trj_all(all_data, all_trj, geometry, save_to):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    save_to (str): address to which the output is saved
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = CellAttributes(all_trj,geometry, splitter='.lammpstrj')
    sim_name = cellAttrs.filename

    print("")
    print(sim_name+" is analyzing...")
    print("")

    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    crds = cell.select_atoms('resid 0') # resid 1 is the atoms creating the chain, resid 0 is crowders
    chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain, resid 0 is crowders


    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    
    # radial direction of the cylindrical coordinate system
    bin_size = 0.1 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.1*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name, save_to)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name, save_to)
    
    # z direction of the cylindrical coordinate system
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.5*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.lcyl
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name, save_to)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name, save_to)
    
    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) #bin size 5 degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name, save_to)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi), edge_name, save_to)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")
               
    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)
        
        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_mon = np.add(zhists_mon, dummy_hist)

        # histogram in theta 
        theta = np.degrees(np.arctan2(crds.positions[:,1], crds.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_crd = np.add(thetahists_crd, dummy_hist)
        
        theta = np.degrees(np.arctan2(chrm.positions[:,1], chrm.positions[:,0])) # in degrees
        dummy_hist, _ = np.histogram(theta, thetaedges)
        thetahists_mon = np.add(thetahists_mon, dummy_hist)
        
    lastname = 'Crd'
    np.savetxt(save_to+sim_name+'-rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(save_to+sim_name+'-zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(save_to+sim_name+'-thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(save_to+sim_name+'-rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(save_to+sim_name+'-zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(save_to+sim_name+'-thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')