import re
import os
from glob import glob

import math
import numpy as np
import datetime
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams, cycler
import matplotlib.lines as lines
from collections import OrderedDict
from matplotlib.legend import Legend
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


import seaborn as sns

import scipy.integrate as integrate
#import scipy.special as special

#from sympy import elliptic_pi

import MDAnalysis as mda 
from MDAnalysis.coordinates.LAMMPS import DumpReader
from MDAnalysis import transformations as mdatransform
from MDAnalysis.analysis import diffusionmap, align, rms

def file_reader(f_names,extensions = [("lammpstrj","trj"),"data"]):
    """
    file_reader returns a sorted list of sorted file names based on the agg function. A trajectroy (trj and lammpstrj extention) constists of N snapshots and a topology (data extention) contains information about simulation box, particle and bond types.

    Parameters:
    f_names (A list of strings): A list of pathes to Lammps trajectroies (trj and lammpstrj extention) and topology (data extention) files.
    
    Return:
    a sorted list of tuples where each tuple has len(extensions) members.

    Requirements:
    re
    """
    f_names_sorted = [None] * len(extensions) # an empty list of len(extensions)
    # a nested list where each sublist has all the files with the same extension:
    for idx, extension in enumerate(extensions): 
        f_names_sorted[idx] = [f_name for f_name in f_names if f_name.endswith(extension)]
    
    def sort_by_int(alphanumeric):
        """
        sort_by_int split an alphanumeric into words and integers.
        
        Reference:
        http://code.activestate.com/recipes/135435-sort-a-string-using-numeric-order/
        
        Parameters:
        alphanumeric (char): an alphanumeric string.

        Return:
        a mixed list of words and integers.
        
        Requirement:
        re
        """

        pattern = re.compile('(\d+)')
        split = pattern.split(alphanumeric)
        mixed = [int(word) if word.isdigit() else word for word in split]
        return mixed
    
    for idx, names_same_ext in enumerate(f_names_sorted):
        f_names_sorted[idx] = sorted(names_same_ext, key = sort_by_int) # Sorting pathnames alphanumerically
    
    f_names = list(zip(*f_names_sorted))
    print("Total number of files is ", len(f_names))
    print("Path to the first tuple of the sorted file: ", f_names[0])
    return f_names

def log_outputs(log_file, geometry):
    cell_attrs = cellAttributes(log_file,geometry)
    output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{strcell_attrs.dcrowd}'
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

class cellAttributes:
    """
    an instance from cellAttributes contains the information about a simulation based on the name of the
    simulation output.
    
    geometry: 'cubic', 'cylindrical', 'slit'
    """

    def __init__(self, filename, geometry, splitter='.bug', printname = False, warning= True):
        self.filepath = filename
        self.warning = warning
        self.filename = 'N'+filename.split('N')[-1].split(splitter)[0]
        if printname:
            print("Filename: "+self.filename)
        self.geometry = geometry
        self.mmon = 1.0 # mass of a monomer
        self.dmon = 1.0 # monomers size/diameter
        self.eps_others = 1.0 # other LJ interaction strengths
        self.mcrowd = 1.0 # mass of a crowder        
        self.nmon_large = 0
        self.dmon_large = 0
        
        self.nmon = 0
        self._attributes()
        
        if self.filepath.endswith(('trj','lammpstrj')):# or (self.filename.endswtich('lammpstrj'))
            self._box_size()

    def _attributes(self):
        """
        I found the followling script for the internet:
        http://code.activestate.com/recipes/135435-sort-a-string-using-numeric-order/
        This function returns the number of crowders written in a file name.
        The problem here is that all the number in the filename except the ncrowd should be the same,
        otherwise the soritng will fail.

        fname is the path to the terjactory file.
        this functions return a list of words devided intro strings and integers.
        """
        pattern = re.compile('([a-zA-Z\_]+)')
        words = pattern.split(self.filename)
        fword = [words.index(word) for word in words if 'N' in str(word)] # Index of the words element containing N
        self.nmon = int(words[words.index('N')+1]) # number of monomers

        if self.geometry  == 'cylindrical':
            # example file name :N#epsilon#r#lz#sig#nc#dtbdump#adump#ens#.bug.lammpstrj[trj,data]
            self.dcrowd = float(words[words.index('sig')+1]) # crowders size/diameter
            self.ncrowd = int(words[words.index('nc')+1]) # number of crowders
            
            try :
                self.ens = int(float(words[words.index('ens')+1])) # ensemble number
            except :
                self.ens = 0
                if self.warning:
                    print("This is an ensemble-averaged file, so ensemble number: {}".format(self.ens))
            
            try :
                self.segnum = int(words[words.index('j')+1]) # ensemble number
            except :
                self.segnum = 0
                if self.warning:
                    print("This is NOT an ALL file, so segnum: {}".format(self.segnum))
        
            self.dwall = 1.0 # wall-forming particles diameter
            self.eps_mon_wall = float(words[words.index('epsilon')+1]) # wall-monomer LJ interaction strength         
            # Be careful about the float number smaller than 1 when use sort_by_int:
            self.dcyl = 2. * float(words[words.index('r')+1]) - 1 # diameter of the cylinder 
            self.lcyl = 2. * float(words[words.index('lz')+1]) # length of the cylinder
            try :
                self.bdump = int(words[words.index('bdump')+1]) # bug dump timestep
            except: 
                self.bdump = 1000 # length of the cylinder
            
            try :
                self.adump = int(words[words.index('adump')+1]) # all dump timestep
            except: 
                self.adump = 1000 # length of the cylinder
                      
            try :
                self.dt = float(words[words.index('dt')+1]) # data sampling timestep
            except:
                self.dt = 0.005 # data sampling timestep
                
            self.cols = "filename,dmon,nmon,dcyl,lcyl,phi_m,rho_m,epsilon,ncrowd,dcrowd,phi_crowd,rho_crowd,ens,fsd,fsd_std,fsd_var,fsd_sem,rflory,rflory_std,rflory_var,rflory_sem,gyr,gyr_std,gyr_var,gyr_sem,nframes\n"    

            self.vol_cell = np.pi * self.dcyl**2 * self.lcyl / 4.0
            
            self.vol_mon = np.pi * self.dmon**3 / 6
            self.phi_mon = (self.nmon * self.vol_mon) / self.vol_cell
            self.rho_mon = self.phi_mon / self.vol_mon
            
            self.vol_crowd = np.pi * self.dcrowd**3 / 6
            self.phi_crowd = (self.ncrowd * self.vol_crowd) / self.vol_cell
            self.rho_crowd = self.phi_crowd / self.vol_crowd
            #self.dcyl_mon = self.dcyl - (self.dwall + self.dmon)/2 # effective diameter of the cylinder for monomers
            #self.dcyl_crd = self.dcyl - (self.dwall + self.dcrowd)/2 # effective diameter of the cylinder for crowder
            
        elif self.geometry  == 'cube':
            self.nmon_large = int(words[words.index('nl')+1]) # nmon large

            self.dmon_large = float(words[words.index('dl')+1]) # dmon large

            self.dcrowd = float(words[words.index('dc')+1]) # crowders size/diameter

            self.ncrowd = int(words[words.index('nc')+1]) # number of crowders

            self.ens = int(words[words.index('ens')+1]) # ensemble number
               
            # Be careful about the float number smaller than 1 when use sort_by_int:
            self.lcube = 2. * float(words[words.index('l')+1]) # diameter of the cylinder    
            #self.phi_c = (4 * np.pi * (self.dcrowd/2.0)**3 * self.ncrowd) / (3 * self.lcube**3) # volume fraction of crowders in bulk (or infinity)
            #self.phi_m = (4 * np.pi * (self.dmon/2.0)**3 * self.nmon) / (3 * self.lcube**3) # volume fraction of monomers
            #self.rho_c = self.ncrowd / self.lcube**3 # density of crowders
            #self.rho_m = self.nmon / self.lcube**3 # density of monomers


        
        #cube_looping_cols='filename','lcube','nmon','nmon_small','dmon_small','nmon_large','dmon_large','phi_m','ncrowd','dcrowd',
        #'phi_crowd','ens','max_dist_x','max_dist_y','max_dist_z','rflory','rflory_std','rflory_std_err','gyr','gyr_std','gyr_std_err'
        else: 
            print("the slit geometry attributes is not defined yet")

    def _box_size(self):
        """
        This piece of script returns the bounds of the simulations box together with the first timestep and
        the total numbero of atoms.
        """
        with open(self.filepath,mode='r') as file: # The same as Chanil's code
            file.readline()
            self.firstframe = int(file.readline()) # timestep
            #print("The first timestep in the file:", self.firstframe)
            file.readline()
            self.natoms = int(file.readline()) # total number of atoms=(monomers+crowders)
            #print("The total number of particles in the first frame:", self.natoms)
            file.readline()
            self.xmin, self.xmax = map(float, file.readline().split(' ')) # x bounds
            #print("x ranges from ",self.xmin," to ", self.xmax)
            self.ymin, self.ymax = map(float, file.readline().split(' ')) # y bounds
            #print("y ranges from ",self.ymin," to ", self.ymax)
            self.zmin, self.zmax = map(float, file.readline().split(' ')) # z bounds
            #print("z ranges from ",self.zmin," to ", self.zmax)
        self.xlen = self.xmax - self.xmin
        self.ylen = self.ymax - self.ymin
        self.zlen = self.zmax - self.zmin
        #else:
         #   print("the file is not a lammps data file.")
          #  print("But, it is assumed to be an  MD-related filename.")
            #self.config_name = self.filename[:-4] # the config name does not have the ens#.
            #del(self.filename)

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

def chain_stats(array, output):
    """
    chain_stats computes the mean, standard deviation, variance, and standard error of the mean of the vector array and write them to the output file.
    
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

def data_file_generator(template, fnames):
    
    if not template[0].endswith('data'):
        raise Exception("The template file is not Lammps data file.")
    
    for fname in fnames:
        if not fname.endswith('lammpstrj'):
            raise Exception("The input files are not Lammps trajectroy files.")
    
    firstline = "LAMMPS data file via PipeLine.data_file_generator, based on LAMMPS (template) data file via write_data, version 29 Oct 2020, timestep = 71000000"
    
    for fname in fnames:
        output = fname.split('/')[-1].split('.lammpstrj')[0]+'.data'
        cell_attrs = cellAttributes(fname,'cylindrical')
        xmax = (cell_attrs.dcyl+1.0) / 2.0
        ymax = (cell_attrs.dcyl+1.0) / 2.0
        zmax = cell_attrs.zlen / 2.0
        with open(template[0],'r') as data_in,\
        open(output,'w') as data_out:
            _ = data_in.readline()
            _ = data_in.readline()
            data_out.write(firstline)
            data_out.write('\n')
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            _ = data_in.readline()
            _ = data_in.readline()
            _ = data_in.readline()
            _ = data_in.readline()
            data_out.write('\n')
            data_out.write(str(-1.0*xmax)+' '+str(xmax)+' xlo xhi \n')
            data_out.write(str(-1.0*ymax)+' '+str(ymax)+' ylo yhi \n')
            data_out.write(str(-1.0*zmax)+' '+str(zmax)+' zlo zhi \n')
            line = data_in.readline()
            while(line):
                data_out.write(line)
                line = data_in.readline()
                
def ensemble(simulations, property_name, geometry, single=True, to_file=True, sep='_', **kwargs):
    """
    ensemble generates an ensemble (dataframe) of simulations (columns) for the given property_name in the geometry of interest. 
    
    Previous names:
    dict_of_ens_dfs
    
    Caution:
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    By default, the indexes of an ensemble (dataframe) are timesteps. For some properties such as histograms, however, the indexes are bin centers. As a result, simulations are tuples where each tuple can have one or more members. For histogram-like properties, each tuple has two members where the second member is the bin edges.
    
    Parameters:
    simulations (list of tuples): a sorted liste of tuples where each tuple at least has one member (the path to the simulation data for the property_name). For histogram-like proeprties, the second member of each tuple is the path to the bin_edges.
    property_name (str): name of the property measured in each of the simulations.
    geometry (str): shape of the simulation box
    single (bool): whether each tuple is a single file or more.
    to_file (bool): whether save to file or not.
    sep (str): a delimiter that used to extract a simulation name from a simulation path.
    
    Return:
    ensembles (dict): a dictionary of ensembles where keys are ensemble names and values are ensembles (dataframes). In each ensemble (value of dict/dataframe), the columns are the simulations of that ensemble.
    
    Requirements:
    Pandas, Numpy, my own cellAttributes class.
    """
    # Unique ensembles with simailar initial paramters    
    ens_names = [cellAttributes(simulation[0],geometry,warning=False).filename.split('ens')[0] for simulation in  simulations] # name of simulaion groups
    ensembles = dict.fromkeys(ens_names)
    if not single :
        bin_edges = dict.fromkeys(ens_names) # bin edges create bin centers which in turn used as indexes of histograms dataframes.
        
    for ens_name in ensembles.keys():
        ens_simulations = []
        for sim_tuple in simulations: #files are organized based on their name
            sim_name = sim_tuple[0].split("/")[-1].split(sep+property_name)[0]
            if sim_name.split('ens')[0] == ens_name:
                simulation = pd.read_csv(sim_tuple[0],names=[sim_name],**kwargs)
                ens_simulations.append(simulation)    
        ensembles[ens_name] = ens_simulations # each key (ensemble) has a list of dataframes (simulations)
        ensembles[ens_name] = pd.concat(ensembles[ens_name],axis=1) # convert the list of dataframes to columns 
        if not single :
            bin_edges[ens_name] = np.around(np.loadtxt(sim_tuple[1]),decimals=2) # warning bin edge rounded based on the binsize variable in extract trj functions
            print("Warning: bin edges, bin centers, and the like are rounded to 2 decimals.")
            bin_centers = np.around((bin_edges[ens_name][:-1] + bin_edges[ens_name][1:]) / 2.0,decimals=2)
            ens_index = pd.Index(bin_centers, dtype = float) # use bin_centers as index
            ensembles[ens_name].set_index(ens_index,inplace=True)
        if to_file:
            cell_attrs = cellAttributes(ens_name,geometry,warning=False)
            output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}nc{cell_attrs.ncrowd}-{property_name}.csv'
            ensembles[ens_name].to_csv(output)
    return ensembles           
            
def group(ensembles, property_name, geometry, to_file=True):
    """
    ens_avg_group averages over the property_name in all the simulations of each ensemble in the ensembles.
    
    Previous names:
    ensemble_avg
    
    Caution: 
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group. 
    
    Parameters:
    ensembles (dict): a dict of ensembles (dataframes) where the number/names of columns in each ensemble is euqal ot the number/names of simulations in that ensemble. The values in each column are the measurements of a given property in each simulation.
    property_name (str): name of the property measured in each of the simulations.
    geometry (str): shape of the simulation box.
    to_file (bool): whether save to file as csv or not
    
    Return:
    ens_avgs (dict): a dict of ensemble-averaged groups where keys are names of ensemble-averaged groups and values are dataframes of ensemble-averaged
    
    Requirements:
    Pandas, Numpy, PipeLine
    
    """
    # Averging over ensembles with simailar initial paramters
    ens_avgs = {} # ensemble-average groups: each group is a bunch of ensemble-averaged simulations.
    for ens_name, ensemble in ensembles.items():
        ens_avg = np.zeros(len(ensemble)) # the length is the number of measurements on the property which is the same in all the simulations in an ensemble.
        ens_count = 0
        for sim_name in ensemble.columns:
            if sim_name.split('ens')[0] == ens_name:
                ens_count += 1
                ens_avg += ensemble[sim_name].values
        ens_avg = ens_avg / ens_count # ensemble-averaged property
        ens_avg = pd.DataFrame(np.array(ens_avg),columns=[f"{ens_name}-ens_avg"]) # converting to a dataframe with ens_avg as values of a column and the ens_name as the name of that column.
        ens_avg.set_index(ensemble.index.values,inplace=True)
        if to_file:
            cell_attrs = cellAttributes(ensemble.columns[0],geometry,warning=False)
            output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}nc{cell_attrs.ncrowd}-{property_name}-ens_avg.csv'
            ens_avg.to_csv(output)
        ens_avgs[ens_name] = ens_avg
    return ens_avgs

def properties(simulations, geometry, to_file=True, **kwargs):
    """
    properties creates a dataframe in whihc there are the input parameters and equilibrium properties of all the simulation in a simulation group.
    
    Previous names:
    df_of_properties

    Caution:
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    
    Parameters:
    simulations (list): pathes to the properties of each simulations 
    geometry (str): shape of the simulation box
    to_file (bool): whether save to file or not
    
    Return:
    A dataframe of the input parameters and equilibrium properties of all the siulations.
    
    Requirements:
    Pandas, PipeLine
    """
    group_properties = []
    for simulation in simulations: #files are organized based on their name
        group_properties.append(pd.read_csv(simulation[0],comment='#',**kwargs)) # a list of individual dataframe for each simulation.
    group_properties = pd.concat(group_properties) # convert the list of dataframe to a dataframe.
    group_properties.reset_index(inplace=True,drop=True)
    print(group_properties.filename[0])
    cell_attrs = cellAttributes(group_properties.filename[0],geometry) # use the first row in the
    print(cell_attrs.dcyl)
    if to_file:
        group_properties.to_csv(f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}-properties.csv')
    return group_properties

def group_properties(properties, geometry, to_file=True, old_naming=True):
    """
    group_properties computes the ensemble-averaged properties of all the ensembles in a group simulation.
    
    Previous names:
    ens_avg_properties
    
    Caution:
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    
    Parameters:
    properties (dataframe): a dataframe of all the individual properties of all the simulations in a simulation group
    geometry (str): shape of the simulation box
    to_file (bool): whether save to file or not
    
    Return:
    a dataframe of ensemble-averged properies of the simulations in a simulation group. 
    
    Requirements:
    Pandas, my own cellAttributes class
    """
    # Names of the ensembles over each of them the ensemble-average is done:
    ens_avg_names = [sim_name.split('ens')[0] for sim_name in  properties.filename] # There many simulation with the same ensemble name.
    ens_avg_names = list(dict.fromkeys(ens_avg_names)) # create a list of unique ensemble names
    
    properties_names = properties.columns.drop(['filename','ens']) # names of properties in the ensemble-averaged properties.
    if old_naming:
        input_parameters = ['dmon','nmon','dcyl','lcyl','vfrc_m','rho_m','epsilon','ncrowd','dcrowd','vfrc_crowd','rho_crowd','nframes'] # input names
    else :
        input_parameters = ['dmon','nmon','dcyl','lcyl','phi_m','rho_m','epsilon','ncrowd','dcrowd','phi_crowd','rho_crowd','nframes'] # input names
    aggs_equilibrium_properties = {sim_property: 'mean' for sim_property in properties.columns.drop(input_parameters+['filename','ens'])} # properties over which ansemble-averaging is done,
    ens_avg_groups = []
    for ens_avg_name in ens_avg_names:
        ens_count = 0
        ens_avg_group = pd.DataFrame(columns=properties_names)
        for sim_name in properties.filename: # append the properties of the simulations that belongs to the ens_avg_group to that group
            if sim_name.split('ens')[0] == ens_avg_name: 
                ens_avg_group = ens_avg_group.append(properties[properties.filename==sim_name][properties_names])
                ens_count += 1
        ens_avg_group = ens_avg_group.groupby(input_parameters).agg(aggs_equilibrium_properties)
        ens_avg_group.reset_index(inplace=True)
        ens_avg_group['filename'] = ens_avg_name+'-ens_avg'
        ens_avg_group['number_of_ens'] = ens_count
        ens_avg_groups.append(ens_avg_group) # a list of ensemble-averaged groups
    ens_avg_groups = pd.concat(ens_avg_groups,axis=0) # convert the list to a dataframe.
    new_columns = [col.strip() for col in ens_avg_groups.columns] # remove extra white spaces
    ens_avg_groups.columns = new_columns # change columns' names to stripped ones.
    new_columns = new_columns[-2:] + new_columns[:-2] # re-order columns to have filename and number_of_ens ase the 1st and 2nd columns
    ens_avg_groups = ens_avg_groups[new_columns]
    ens_avg_groups.reset_index(drop=True, inplace=True)
    if to_file:
        cell_attrs = cellAttributes(properties['filename'][0],geometry)
        output = f'N{cell_attrs.nmon}D{cell_attrs.dcyl}ac{cell_attrs.dcrowd}-properties-ens_avg.csv'
        ens_avg_groups.to_csv(output)
    return ens_avg_groups

def cyl_sumrule_norms_ens_evg(all_in_one):
    """
    cyl_sumrule_norms_ens_evg rescales some of the physical properties in an ensemble-average properties dataframe of the cylindrical
    sum rule project.
    
    Parameters:
    all_in_one: a dataframe of all the ensemble-averaged properties.
    
    Return:
    all_in_one with several new columns.
    """
    
    all_in_one['fsd_normalized'] = 0.0
    all_in_one['gyr_normalized'] = 0.0
    all_in_one['rflory_normalized'] = 0.0
    all_in_one['vfrc_c_normalized'] = (all_in_one['dmon'] * all_in_one['vfrc_crowd'])/all_in_one['dcrowd'] 
    group_properties = ['nmon','dcyl','dcrowd'] # this is the way ensemble group defined in the cylindrical sum rule project
    groups = [list(group_name) for group_name in list(all_in_one.groupby(group_properties).size().index)]
    for nmon, dcyl, dcrowd in groups:
        condition = (all_in_one['nmon'] == nmon) & (all_in_one['dcyl'] == dcyl) & (all_in_one['dcrowd'] == dcrowd) & (all_in_one['vfrc_crowd'] == 0)
        fsd = all_in_one[condition]['fsd'].values[0]
        gyr = all_in_one[condition]['gyr'].values[0]
        rflory = all_in_one[condition]['rflory'].values[0]
        cond_normal = (all_in_one['nmon'] == nmon) & (all_in_one['dcyl'] == dcyl) & (all_in_one['dcrowd'] == dcrowd)
        all_in_one.loc[cond_normal,'fsd_normalized'] = all_in_one.loc[cond_normal,'fsd'] / fsd
        all_in_one.loc[cond_normal,'gyr_normalized'] = all_in_one.loc[cond_normal,'gyr'] / gyr
        all_in_one.loc[cond_normal,'rflory_normalized'] = all_in_one.loc[cond_normal,'rflory'] / rflory
        
    all_in_one['vfrc_crowd'] = all_in_one['vfrc_crowd'].round(decimals=3)
    all_in_one['vfrc_crowd_eff'] = (np.pi * all_in_one['dcrowd'] ** 3 / 6) * all_in_one['ncrowd'] / ((np.pi / 4 * (all_in_one['dcyl']-all_in_one['dcrowd']) ** 2) * all_in_one['lcyl'])
    all_in_one['vfrc_c_eff_normalized'] = (all_in_one['dmon'] * all_in_one['vfrc_crowd_eff']) / all_in_one['dcrowd'] 
    return all_in_one

def all_dfs_in_one(property_files, property_name, ens_avg=False, to_file=True, norm_func=None, **kwargs):
    """
    all_dfs_in_one merges the all the dataframes into one dataframe. df_files can be a list of ensmeble csvs or ensemble_average group csv.
    
    Caution: 
    A simulation group usually results in a graph or curve for the project and refers to a collection of simulations that all have the same values for one or several input parameters of the project.
    An ensemble is a collection of themodynamically-equivalent simulations that differs only in their random number seeds, initial conditions, or boundary conditions but have the same input parameters. In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    An ensemble-averaged group is an average over all the simulations in an ensemble and usually gives a data point.
    If there are N esmebles, each with M simulations, then there are N ensemble-average groups and N*M simulations in the simulation group.
    
    Paremeters:
    property_files: filenames of dataframes.
    property_name (str): name of the property.
    ens_avg (bool): whether ensemble files or ensemble-averaged files.
    to_file (bool): whether save to file or not
    
    Return:
    A dataframe of all the properties of all the ensembles.
    
    Requirements:
    Pandas, my own cellAttributes class
    """
    all_in_one = []
    for property_file in property_files: #files are organized based on their name
        all_in_one.append(pd.read_csv(property_file[0],comment='#',**kwargs))# each list member is a tuple with one member; see file_reader function.
    all_in_one = pd.concat(all_in_one)
    all_in_one.reset_index(inplace=True,drop=True)
    if (ens_avg == True) and (norm_func != None) and (property_name == 'properties'):
        all_in_one = norm_func(all_in_one)
        df_type = 'all_in_one-ens_avg-normalized.csv'
    elif ens_avg:
        df_type = 'all_in_one-ens_avg.csv'
    else:
        df_type = 'all_in_one.csv'
    if to_file:
        all_in_one.to_csv(property_name+'-'+df_type)
    return all_in_one

def hist_segments_reader(segments, geometry, splitter='all'):
    """
    hist_segments_reader returns a dict in which the keys are the histogram names and the values are the dataframes of the histograms with bin cneters as their indexes. In each dataframe, the number of columns is equal to the number of extracted segments. Histogram can be done for any propeorty of interest.
    
    Caution: 
    Extraction is a phase in which a physical property is calculated directly from Lammps simulation dump files, usually with the help of Python MDAnalysis package. 
    There is a difference between the extracted segments of a static quantity such the histgram or the chain length and a dynamic quantity such as rmsd. For dyanmic quanitties, operations are done on a sequence of time steps so the segments should also be in sequence. 
    There is also a difference between a counting problem such as histograms (oe number densities, or end-to-end distribution) and other problem such measuring the chain size or average entropy. For counting problems, the measurments in different time steps can collected/accomumlated in one or severla bins so an individual measurement in a time step is not important but the distirbution of measurmetn plays role. However, in other measurment such a chain size, we need the measurments in all time step to calculate a mean or variance, for instance; here, the order of time step is still unimportnat (in contast to a dynamic property). Putting differently, the index of a counting property is bin centers but the index of other proeprties is a time step. 
    In counting problems, each segment is a list with two member, the first member is the distribution of measurement in bins and the second one is the bin edges. The bin edges are the same for all the segments of a simulation extracted property.
    For counting problems the measurments in the segments (columns) are summed (vertically) to give the measurments in the whole simulation. However, this is not correct for other problems; in this cases, they should merged horizontally. 
    
    Works for future:
    Combine hist_segments_merger and hist_segments_reader and add functionality to combine them vertically or horizontally.

    Parameters:
    segments (list of tuples): a sorted list of tuples where each tuple has two members (bin edges (2nd memeber) and the counts in bin edges(1st memeber)). Each tuple is an extracted segment. 
    geometry (str): shape of the simulation box.
    
    Return:
    a dictionary (segments_of_simulations) that has the names of simulations as its keys and the pandas dataframes (the histogram of that segments) as its value. In each dataframe (dict value), the columns are the segments of a simulation (dict key). 
    
    Requirements:
    Pandas, Numpy, my own cellAttributes class.
    """
    # Segment names in an ensemble
    sim_names = []
    for segment in segments:
        cell_attrs = cellAttributes(segment[0],geometry, splitter='.all',warning=False) # the distribution
        sim_name = cell_attrs.filename.split('ens')[0]+'ens'+str(cell_attrs.ens) # finding name of the simulation the segment belongs to.
        sim_names.append(sim_name) 
    segments_of_simulations = dict.fromkeys(sim_names) # create a dict of unique simulation names.
    bin_edges = dict.fromkeys(sim_names) # the ditribution and the bine edges has the same name. 
    for sim_name in segments_of_simulations.keys():
        segments_per_simulation = []
        for segment in segments: #files are organized based on their name
            cell_attrs = cellAttributes(segment[0],geometry, splitter='.all',warning=False)
            seg_name = cell_attrs.filename # the name of a segment.
            if sim_name in seg_name: # check whether the segment belongs to the simulation with the name sim_name
                segment_df = pd.read_csv(segment[0],names=[seg_name]) # df of a segment with its name as the column name
                segments_per_simulation.append(segment_df)
                bin_edges[sim_name] = np.around(np.loadtxt(segment[1]),decimals=2) # warning bin edges rounded based on the binsize variable in extraction phase. Since the bin edge is the same for all the segments of a simulation. The bin_edge dataframe is overwrite several times (ewual to the number of segments) and has only one column.
        segments_of_simulations[sim_name] = pd.concat(segments_per_simulation,axis=1) # Each key (simulation) has a list of dataframes as its value. These dataframes are the segments of a simulation. Segments are concatnated along vertically, so the number of colums are equal to the number of segemtns. For problems other than counting, axis can be set to 0 but then the index should reset.
        bin_centers = np.around((bin_edges[sim_name][:-1] + bin_edges[sim_name][1:]) / 2.0, decimals=2) # warning bin centeres rounded based on the bin edges variable in extract trj functions
        segments_of_simulations[sim_name].set_index(bin_centers,inplace=True)
    return segments_of_simulations

def hist_segments_merger(segments_of_simulations, property_name, geometry, to_file=True):
    """
    hist_segments_merger merges extracted segments resulted from the extraction done on simulation segments to measure the histogram_name.
    
    Caution:
    Extraction is a phase in which a physical property is calculated directly from Lammps simulation dump files, usually with the help of Python MDAnalysis package. 
    There is a difference between the extracted segments of a static quantity such the histgram or the chain length and a dynamic quantity such as rmsd. For dyanmic quanitties, operations are done on a sequence of time steps so the segments should also be in sequence. 
    There is also a difference between a counting problem such as histograms (oe number densities, or end-to-end distribution) and other problem such measuring the chain size or average entropy. For counting problems, the measurments in different time steps can collected/accomumlated in one or severla bins so an individual measurement in a time step is not important but the distirbution of measurmetn plays role. However, in other measurment such a chain size, we need the measurments in all time step to calculate a mean or variance, for instance; here, the order of time step is still unimportnat (in contast to a dynamic property). Putting differently, the index of a counting property is bin centers but the index of other proeprties is a time step. 
    In counting problems, each segment is a list with two member, the first member is the distribution of measurement in bins and the second one is the bin edges. The bin edges are the same for all the segments of a simulation extracted property.
    For counting problems the measurments in the segments (columns) are summed (vertically) to give the measurments in the whole simulation. However, this is not correct for other problems; in this cases, they should merged horizontally. 
    
    Works for future:
    Combine hist_segments_merger and hist_segments_reader and add functionality to combine them vertically or horizontally.
    
    Parameters:
    segments_of_simulations (dict of dataframes): a dictonary in which the keys are the names of simulations and the values are the dataframes. In each dataframe (dict value), the columns are the segments of a simulation (dict key). 
    property_name (str): the name of property for which merging is done on extraction segments.
    geometry (str): shape of the simulation box
    to_file (bool): whether save to file or not
    
    Return:
    a dict in which simulation names are keys and dataframes of whole porperty_name files are values. 
    
    Requirements:
    Pandas, Numpy, PipeLine.
    """
    # Averging over ensembles with simailar initial paramters
    simulations = {} 
    for sim_name, segments_per_simulation in segments_of_simulations.items():
        simulation = segments_per_simulation.sum(axis=1) # added all the segments along the columns (vertically) because it is a coutning problem.
        simulation = simulation.to_frame(name=sim_name)
        simulation.set_index(segments_per_simulation.index,inplace=True)
        if to_file:
            simulation.to_csv(f'{sim_name}-{property_name}.csv')
        simulations[sim_name] = simulation
    return simulations

def simulation_from_segments(filenames,property_name,geometry,extentions,**kwargs):
    """
    simulation_from_segments wraps the three steps for merging extracted segments resulted from the extraction done on simulation segments to measure the property_segments_name.
    
    Caution: Extraction is a phase in which a physical property is calculated directly from Lammps simulation dump files, usually with the help of Python MDAnalysis package.
    
    Parameters:
    filenames (list): a list of extracted segments (filenames).
    property_name (str): the name of property for which extraction is done on simulation segments.
    geomtetry (str): the geometry of the simulatiob box.
    extensions (list): the list of files extensions for property_name
    **kwargs: the kwargs being passed to segment_reader
    
    Return:
    a dict with simulation names as keys and dataframes of whole porperty_name files as values. 
    
    Requirements:
    PipeLine
    """
    # csvs of segments extracted from csv files and then ordered by name
    segments = file_reader(filenames,extensions=extentions) 
    # df from each segments saved in dict with filename as key:
    # to_file=False means that all-segment-in-one dataframes are not save.
    # a all-segments-in-one dataframe is a simluation dataframe in which the number of columns is equal to the number of the segments of that simulation.
    segments_of_simulations = hist_segments_reader(segments,geometry,**kwargs)
    # df from all the segments of one ensemble. 
    return hist_segments_merger(segments_of_simulations,property_name,geometry)

def cylinder_write(output,cellAttrs):
    """
    This function writes the information stored in a cellAttrs object into a file.
    Parameters: 
    cellAttrs (an object from cellAtttribute class):this object contains information about the parameters used as inout in a simulation.
    ensemble (python file object in "a" mode): the output file that contains the statistics about a simulation.
    
    Returns: ---
    
    Requirements:
    cellAttributes class
    
    """
    output.write("{},{},{},{},{},".format(cellAttrs.filename,cellAttrs.dmon,cellAttrs.nmon,cellAttrs.dcyl,cellAttrs.lcyl))
    output.write("{},{},{},{},{},".format(cellAttrs.phi_mon,cellAttrs.rho_mon,cellAttrs.eps_mon_wall,cellAttrs.ncrowd,cellAttrs.dcrowd))
    output.write("{},{},{},".format(cellAttrs.phi_crowd,cellAttrs.rho_crowd,cellAttrs.ens))

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
    to_file: wherther write the graphs to file or not.
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
    data_var_err = math.sqrt(2/(nframe-1)) * data_var  # error in Bias-corrected sample variance

    data_sem = np.sqrt(data_var / nframe) # standard error of the mean (SEM) neglecting any correlations
    data_sem_err = math.sqrt(1/(2*(nframe-1)))* data_sem # error in SEM

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
    si_initial_err = math.sqrt(2/(nframe-1)) * si_initial # error in inintial SI

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
            bvar_err = np.append(bvar_err, math.sqrt(2/(nblock[-1]-1))*bvar[-1] )
            bsem = np.append(bsem, np.sqrt( bvar[-1] / nblock[-1] ) )   # Estimate of SEM
            bsem_err = np.append(bsem_err, math.sqrt((1/(2*(nblock[-1]-1))))*bsem[-1] )  # error in SEM
            si = np.append(si, tblock[-1] * bvar[-1] / data_var) # Statistical inefficiency (SI)
            si_err = np.append(si_err, math.sqrt((1/(2*(nblock[-1]-1))))*si[-1] ) # error in SI
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
    The function calculate the end-to-end distance of a polymer.
    
    Parameters
    r (natoms*ndim numpy ndarray):the positions of the N atoms within an atom group in a frame. This positions are arreanged by atom number
    form 1 to N.

    
    Return:
    the end-to-end distance: numpy.float
    
    Requirements:
    MDAnalysis and numpy. This functions does not use MDAnalysis directly, but it use the fact that the vector r is ordered by the
    monomoer id from one end of the chain to the other end of that.

    Caution:
    Since all the atoms composing the polymer have the same mass, the masses of them are not considered in
    the following lines. If the masses are diffrent, you need to use average method from the numPy.
    """
    
    com = np.mean(r, axis=0) # center of mass of the chain.
    r = r - com # positions of monomers in the com reference frame.
    r = r[-1] - r[0] # end-to-end vector; the position vector r is arranged based on the id of the monomers (see Caution above).
    return np.linalg.norm(r)

def max_distance(r):
    import numpy as np
    """
    This function returns the maximum ditance in each of the three Cartesian direction.

    Parameters
    r (natoms*ndim numpy ndarray):the positions of the N atoms within an atom group in a frame. This positions are arreanged by atom number
    form 1 to N.
    
    Return:
    xmax, ymax, zmax in each timestep (snapshot): numpy.float
    
    Requirements:
    MDAnalysis and numpy. This functions does not use MDAnalysis directly, but it use the fact that the vector r is ordered by the
    monomoer id from one end of the chain to the other end of that.

    Caution:
    Since all the atoms composing the polymer have the same mass, the masses of them are not considered in
    the following lines. If the masses are diffrent, you need to use average method from the numPy.
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
    return np.ptp(coordinates[:,axis])

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
    
def bin_ceate(sim_name, bin_size, lmin, lmax, edge_name):
    """
    bin_ceate produces arrays of bins and histograms

    inputs:
    sim_name: the name of the simulation.
    bin_size (float): size of each bin.
    lmin (float): lower bound of the system in the direction of interest.
    lmax (float): upper bound of the system in the direction of interest.
    edge_name: name of the variable for which the histogram is computed.
    
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
    np.savetxt(sim_name+'_'+edge_name+'.csv', bin_edges, delimiter=',')
    return bin_edges, hist_collectors

def extract_trj_bug(fname, geometry):
    """
    The function does all the analysis on the a whole trajectory (simulation),

    Parameters:
    fname (string): the name of trajectory or simulation.
    geometry (string): the name of the geometry of the simulation box.

    Requirements:
    All the above-defined fucntions and classes, MDAnalysis.
    
    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.
    For other goemteries you need to redfined them.
    """  
    print("This script for cylindrical geomery")
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(fname[0],geometry)
    sim_name = cellAttrs.filename
    #ensGroup = "N"+str(cellAttrs.nmon)+'D'+str(cellAttrs.dcyl)+"_"
    outfile = sim_name+"_properties.csv"
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

    cell = mda.Universe(fname[1], fname[0], format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    edge_name = 'rEdges'
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.2 * min(cellAttrs.dmon,cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    redges, rhist_collectors = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)

    # z direction of the cylindrical coordinate system
    edge_name = 'zEdges'
    #bin_size = cellAttrs.dcrowd / 2.
    #bin_size = 0.2 # in a_c=1.0 unit for a_c
    bin_size = 0.5 * min(cellAttrs.dmon,cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.5*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = cellAttrs.lcyl / 2.
    zedges, zhist_collectors = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)

    # theta of the cylindrical coordinate system
    edge_name = 'thetaEdges'
    bin_size = np.degrees(np.pi/24) # in degrees
    print(f"Bin size for theta in the cylindrical geometry is set to {bin_size} in degrees.")
    thetaedges, thetahist_collectors = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)

    # Chain end-to-end size distribution
    edge_name = 'rFloryEdges'
    lmax = cellAttrs.nmon * cellAttrs.dmon # The contour size of the chain
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.5 * cellAttrs.dmon # in a_c=1.0 unit; for a_c<a_mon_small; check this
    #print("Bin size for chain PDF of end-to-end vector is set to {} in a_c=1 units.".format(bin_size))
    print(f"Bin size for the PDF of end-to-end vector is set to 0.5*cellAttrs.dmon={bin_size} in a_m=1 units.")
    rFloryedges, rFloryhist_collectors = bin_ceate(sim_name, bin_size, 0, lmax, edge_name)

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

    np.savetxt(sim_name+'_rHists.csv', rhist_collectors, delimiter = ',')
    np.savetxt(sim_name+'_zHists.csv', zhist_collectors, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists.csv', thetahist_collectors, delimiter = ',')
    np.savetxt(sim_name+'_rFloryHists.csv', rFloryhist_collectors, delimiter = ',')
    np.savetxt(sim_name+'_fsd_t.csv', fsd_t, delimiter = ',')
    np.savetxt(sim_name+'_rFlory_t.csv', rFlory_t, delimiter = ',')
    np.savetxt(sim_name+'_gyr_t.csv', gyr_t, delimiter = ',')

    with open(outfile, mode="a") as ensemble:
        chain_stats(fsd_t, ensemble) # fsd,fsd_std,fsd_var,fsd_sem,
        chain_stats(rFlory_t, ensemble) # rflory,rflory_std,rflory_var,rflory_sem,
        chain_stats(gyr_t, ensemble) # gyr,gyr_std,gyr_var,gyr_sem,
        ensemble.write("{}\n".format(cell.trajectory.n_frames))

    error_calc_block(fsd_t, sim_name+'_fsd')
    error_calc_block(rFlory_t, sim_name+'_rFlory')
    error_calc_block(gyr_t, sim_name+'_gyr')
    #This function comptues the rmsd of the trajectory.
    print('done.') 
    
def bug_trj_rmsd(fname, geometry):
    """
    This function comptues the rmsd of the trajectory.
    """
    print("Setting the name of analyze file...\n")
    print("Doing RMSD analysis...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(fname[0], geometry = geometry)
    
    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit
    
    sim_name = cellAttrs.filename
    cell = mda.Universe(fname[1], fname[0], format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    cell.transfer_to_memory(step=50,verbose=False)
    
    aligncell = align.AlignTraj(cell, cell, select = 'resid 1', filename = sim_name + '.dcd').run()
    matrix = diffusionmap.DistanceMatrix(cell, select = 'resid 1').run()
    np.save(sim_name+'_rmsdMatrix.npy',matrix.dist_matrix) 
                
    
def extract_all_trj(all_data, all_trj, geometry):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    The histograms are for the cylindrical coordinate system in the geometry=cylindrical.For other goemteries you need to redfined them.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
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
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.2 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for r in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(sim_name, rbin_size, 0.0, lmax, edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(sim_name, bin_size, 0.0, lmax, edge_name)

    # z direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = 0.5 * min(cellAttrs.dmon, cellAttrs.dcrowd)
    print(f"Bin size for z in the cylindrical geometry is set to 0.2*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.lcyl
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(sim_name, zbin_size, -1.0 * lmax, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(sim_name, bin_size, -1.0 * lmax, lmax, edge_name)

    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/24) # in degrees
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)

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
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')

def extract_all_trj_binsize(all_data, all_trj, geometry, rbin_factor=0.2, zbin_factor=0.2):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.

    Caution:
    1. The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.
    2. The diamters/sizes of monomers, crowders, cylinders, and slits as well as the lengths of cylinders, cubes, and slits, all, are either a multiples of 0.5 or 1.0; as a result, the bin size in z or r direction is better to be a mutiple of either 0.5 or 1.0, or possible 0.1. Otherwise, there are some ossciallry pattern (possibly due to round-off error) inthe volume fraction method of Distribtuion class.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)
    crds = cell.select_atoms('resid 0') # resid 1 is the atoms creating the chain, resid 0 is crowders
    chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain, resid 0 is crowders
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = rbin_factor * min(cellAttrs.dmon, cellAttrs.dcrowd)
    r_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor'+str(rbin_factor)+'ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for r in the cylindrical geometry is set to {rbin_factor}*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    #lmax = 0.5 * (cellAttrs.dcyl - max(cellAttrs.dmon, cellAttrs.dcrowd))
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(r_sim_name, bin_size, 0.0, lmax,edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(r_sim_name, bin_size, 0.0, lmax, edge_name)

    # z direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = zbin_factor * min(cellAttrs.dmon, cellAttrs.dcrowd)
    z_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor'+str(zbin_factor)+'ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for z in the cylindrical geometry is set to {zbin_factor}*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.lcyl
    #lmax = 0.5 * (cellAttrs.lcyl - max(cellAttrs.dmon, cellAttrs.dcrowd))
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(z_sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    #zedges, zhists_crd = bin_ceate(z_sim_name, bin_size, 0, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(z_sim_name, bin_size, -1.0 * lmax, lmax,edge_name)
    #_ , zhists_mon = bin_ceate(z_sim_name, bin_size, 0, lmax,edge_name)

    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # is 6 degrees
    theta_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor1.0ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for theta in the cylindrical geometry is set to 1.0*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(theta_sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(theta_sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)

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
        #zpos = np.abs(crds.positions[:,2]) # absolute value of the z component of position of each crowder - why absolute values? Bacause of PBC along z.
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)
        
        zpos = chrm.positions[:,2] # z component of position of each monomer
        #zpos = np.abs(chrm.positions[:,2]) # absolute value of the z component of position of each monomer - why absolute values? Bacause of PBC along z.
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
    np.savetxt(r_sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(z_sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(theta_sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')
    
    lastname = 'Mon'
    np.savetxt(r_sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(z_sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(theta_sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')


def extract_all_trj_binsize_dynamic(all_data, all_trj, geometry, rbin_factor=0.2, zbin_factor=0.2):
    """
    extract_all_trj does all the measurements on the a bug's whole trajectory (simulation) and returns a group of files as the result of the measurements.
    
    

    Caution:
    1. The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.
    2. The diamters/sizes of monomers, crowders, cylinders, and slits as well as the lengths of cylinders, cubes, and slits, all, are either a multiples of 0.5 or 1.0; as a result, the bin size in z or r direction is better to be a mutiple of either 0.5 or 1.0, or possible 0.1. Otherwise, there are some ossciallry pattern (possibly due to round-off error) inthe volume fraction method of Distribtuion class.
    
    Parameters:
    all_data (str): Lammps data (topology) file of the simulation.
    all_trj (str): Lammps dump (trajectory) file of the simulation.
    geometry (str): Shape of the simulation box.
    
    Returns:
    a group of files with different formats, depending on the type of measurements.

    Requirements:
    MDAnalsis, PipeLine 
    """
    
    print("Setting the name of analyze file...\n")
    today = datetime.date.today().strftime('%Y%m%d')
    cellAttrs = cellAttributes(all_trj,geometry, splitter='.lammpstrj')
    time_unit = 1.0 # time unit = dmon * sqrt(mmon * epmon)
    lj_nstep = cellAttrs.bdump # Sampling steps via dump command in Lammps
    lj_dt = cellAttrs.dt # Timestep during the sampling process in Lammps
    simulation_dt = lj_nstep * lj_dt * time_unit

    cell = mda.Universe(all_data, all_trj, format = 'AMIRLAMMPSDUMP', atom_style = "id resid type x y z", dt = simulation_dt)

    length_step = 10 * max(cellAttrs.dmon, cellAttrs.dcrowd)

    externalRadius = 0.5 * cellAttrs.dcyl
    #zMax = length_step # the radii of the three types of particle.
    #zMin  = -1 * length_step # the radii of the three types of particle.
    #near_chrm_cog = cell.select_atoms(f'cyzone {externalRadius} {zMax} {zMin} resid 1', updating=True)
    chrm_initial = cell.select_atoms('resid 1')
    #near_center = cell.select_atoms(f'prop abs z <= {length_step}', updating=True)

    #crds = cell.select_atoms('resid 0') # resid 1 is the atoms creating the chain, resid 0 is crowders
    #crds = cell.select_atoms('prop abs z <= 50 resid 1', updating=True) # resid 1 is the atoms creating the chain, resid 0 is crowders
    #crds = near_chrm_cog.select_atoms('resid 0', updating=True)
    #chrm = cell.select_atoms('resid 1') # resid 1 is the atoms creating the chain, resid 0 is crowders
    #chrm = cell.select_atoms('prop abs z <= 50 resid 1', updating=True) # resid 1 is the atoms creating the chain, resid 0 is crowders
    #chrm = near_chrm_cog.select_atoms('resid 1', updating=True)
    print ("Caution:The histograms are for the cylindrical coordinate system in the geometry=cylindrical. For other goemteries you need to redfined them.")
    # Monomer local (number) density: 
    # radial direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = rbin_factor * min(cellAttrs.dmon, cellAttrs.dcrowd)
    r_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor'+str(rbin_factor)+'ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for r in the cylindrical geometry is set to {rbin_factor}*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    lmax = 0.5 * cellAttrs.dcyl
    #lmax = 0.5 * (cellAttrs.dcyl - max(cellAttrs.dmon, cellAttrs.dcrowd))
    edge_name = 'rEdgesCrd'
    redges, rhists_crd = bin_ceate(r_sim_name, bin_size, 0.0, lmax,edge_name)
    edge_name = 'rEdgesMon'
    _ , rhists_mon = bin_ceate(r_sim_name, bin_size, 0.0, lmax, edge_name)

    # z direction of the cylindrical coordinate system
    #bin_size = cellAttrs.dcrowd / 2.
    bin_size = zbin_factor * min(cellAttrs.dmon, cellAttrs.dcrowd)
    z_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor'+str(zbin_factor)+'ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for z in the cylindrical geometry is set to {zbin_factor}*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    #lmax = length_step
    lmax = 0.5 * cellAttrs.lcyl
    edge_name = 'zEdgesCrd'
    zedges, zhists_crd = bin_ceate(z_sim_name, bin_size, -1.0 * lmax, lmax, edge_name)
    #zedges, zhists_crd = bin_ceate(z_sim_name, bin_size, 0, lmax, edge_name)
    edge_name = 'zEdgesMon'
    _ , zhists_mon = bin_ceate(z_sim_name, bin_size, -1.0 * lmax, lmax,edge_name)
    #_ , zhists_mon = bin_ceate(z_sim_name, bin_size, 0, lmax,edge_name)

    # theta of the cylindrical coordinate system
    bin_size = np.degrees(np.pi/36) # is 6 degrees
    theta_sim_name = cellAttrs.filename.split("ens")[0]+'binFactor1.0ens'+cellAttrs.filename.split("ens")[1]
    print(f"Bin size for theta in the cylindrical geometry is set to 1.0*min(cellAttrs.dmon,cellAttrs.dcrowd)={bin_size} in a_m=1 units.")
    edge_name = 'thetaEdgesCrd'
    thetaedges, thetahists_crd = bin_ceate(theta_sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)
    edge_name = 'thetaEdgesMon'
    _ , thetahists_mon = bin_ceate(theta_sim_name, bin_size, -1*np.degrees(np.pi), np.degrees(np.pi)-bin_size, edge_name)

    # check if any of the histograms are empty or not.
    if any([rhists_crd.any() != 0, rhists_mon.any() != 0, zhists_crd.any() != 0, zhists_mon.any() != 0, thetahists_crd.any() != 0, thetahists_mon.any() != 0]):
        raise ValueError("One of the histogram collectors are not empty!")

    for ts in cell.trajectory: # the length of the for loop is equal to number of snapshots (configurations or frames)
        #number density in the cell's frame of reference
        # histogram in r direction
        chrm_fsd =  fsd(chrm_initial.positions)
        zMax = 0.5 * (chrm_fsd+cellAttrs.dmon) # the radii of the three types of particle.
        zMin  = -0.5 * (chrm_fsd+cellAttrs.dmon) # the radii of the three types of particle.
        around_chrm = cell.select_atoms(f'cyzone {externalRadius} {zMax} {zMin} resid 1', updating=True)
        crds = around_chrm.select_atoms('resid 0', updating=True)
        chrm = around_chrm.select_atoms('resid 1', updating=True)


        rpos = np.linalg.norm(crds.positions[:,:2], axis = 1) # r component of position of each crowder
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_crd = np.add(rhists_crd, dummy_hist)

        rpos = np.linalg.norm(chrm.positions[:,:2], axis = 1) # r component of position of each monomer
        dummy_hist, _ = np.histogram(rpos, redges)
        rhists_mon = np.add(rhists_mon, dummy_hist)

        # histogram in z direction
        zpos = crds.positions[:,2] # z component of position of each crowder
        #zpos = np.abs(crds.positions[:,2]) # absolute value of the z component of position of each crowder - why absolute values? Bacause of PBC along z.
        dummy_hist, _ = np.histogram(zpos, zedges)
        zhists_crd = np.add(zhists_crd, dummy_hist)

        zpos = chrm.positions[:,2] # z component of position of each monomer
        #zpos = np.abs(chrm.positions[:,2]) # absolute value of the z component of position of each monomer - why absolute values? Bacause of PBC along z.
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
    np.savetxt(r_sim_name+'_rHists'+lastname+'.csv', rhists_crd, delimiter = ',')
    np.savetxt(z_sim_name+'_zHists'+lastname+'.csv', zhists_crd, delimiter = ',')
    np.savetxt(theta_sim_name+'_thetaHists'+lastname+'.csv', thetahists_crd, delimiter = ',')

    lastname = 'Mon'
    np.savetxt(r_sim_name+'_rHists'+lastname+'.csv', rhists_mon, delimiter = ',')
    np.savetxt(z_sim_name+'_zHists'+lastname+'.csv', zhists_mon, delimiter = ',')
    np.savetxt(theta_sim_name+'_thetaHists'+lastname+'.csv', thetahists_mon, delimiter = ',')
    print('done.')
    
def color_patcher(colors):
    from matplotlib import rcParams, cycler
    import matplotlib.patches as mpatches
    """
    return color patches for legend.
    
    Params:
    color: matplotlib color object
    Return"
    color_patches; matplotlib color patch object
    Dependencies:
    matplotlib
    """
    color_patches = []
    for kwargs in cycler(color=colors):
        color_patches.append(mpatches.Patch(**kwargs))
    return color_patches

def rdf_ideal_plotter(ax, bin_edges, nmon, dmon, rdf=True, dmon_large=0, *args, **kwargs):
    """
    plot the probability distribution function (pdf) of the end-to-end vector of an ideal linear chain
    
    inputs:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an ideal linear chainm
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    rdfFactor = 4*np.pi*bin_centers**2
    if rdf == False:
        rdfFactor = 1.0   
    pdf_values = rdfFactor*np.exp((-3.0*bin_centers**2)/(2.0*nmon*dmon**2))
    pdf_values = pdf_values/np.sum(pdf_values) # normalization
    bin_centers = bin_centers / dmon_large # normalized bins
    return ax.plot(bin_centers, pdf_values, **kwargs)

def pdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """
    pdf: probability distribution function of the end-to-end distance vector.
    pdf draws a step plot.

    inputs:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax
    return:
    a plot of (normalized) probability distribution function of the end-to-end distance vector.
    """
    histo_collections = histo_collections / np.sum(histo_collections) # normalization
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    rdfFactor = 4*np.pi*(bin_centers)**2
    histo_collections = np.divide(histo_collections,rdfFactor) 
    histo_collections = histo_collections/np.sum(histo_collections) # normalization
    bin_centers = bin_centers / dmon_large
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)

def rdf_plotter(ax, histo_collections, bin_edges, dmon_large=0, **kwargs):
    """
    rdf: radial distribution function or the probability distribution function of *the magnitude of the end-to-end vector* 
    (or *the end-to-end distance*) in between R and R+dR

    inputs:
    ax: matplotlib axis object.
    histo_collections: the collected histogram in the direction of interest over whole universe.
    bin_edges: an array of bins' edges.
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.
    **kwargs: the parameters that are used in ax 

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    histo_collections = histo_collections / np.sum(histo_collections) # normalization
    bin_centers = bin_centers/dmon_large # normalized bins
    return ax.step(bin_centers, histo_collections, where='mid', **kwargs)

def rdf_real_plotter(ax, bin_edges, rFlory, nmon, dmon,rdf=True, dmon_large=0, *args, **kwargs):
    """
    plot the probability distribution function (pdf) of the end-to-end vector of a real linear chain
    
    inputs:
    ax: matplotlib axis object.
    bin_edges: an array of bins' edges.
    nmon: number of monomers
    dmon: size of monomer (or the kuhn segment)
    rdf: The probability distribution of the end-to-end distance of an ideal linear chainm
    dmon_large = size of the pair of monomers at the two ends of the chain. This is used for scaling bin centers defined below.

    return:
    a plot of (normalized) radial probability distribution function of the magnitude of the end-to-end distance vector.
    """
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    #rFlory = 1.12 * nmon ** 0.588 * dmon  # chain size or the root-mean-square end-to-end distancs or Flory radius
    rdfFactor = 4*np.pi*(bin_centers)**2
    if rdf == False:
        rdfFactor = 1.0   
    pdf_values = rdfFactor*(bin_centers/rFlory)**0.28*np.exp(-1.206*(bin_centers/rFlory)**2.43) # default coeff = 1.206
    pdf_values = pdf_values/np.sum(pdf_values) # normalized
    bin_centers = bin_centers/dmon_large # normalized bins
    return ax.plot(bin_centers, pdf_values, **kwargs)


def entropic_energy(histo_collections, beta = 1.0):
    """
    this is the effective free energy potential as the function of the end-to-end distance. Here, we use the end-to-end distance as the action parameters and write a Langevin equation for that.
    
    inputs:
    histo_collections: radial distribution function or the probability distribution function of *the magnitude of the end-to-end vector* 
    (or *the end-to-end distance*) between R and R+dR
    
    beta = k_B*T is the inverse of the thermal energy.

    return:
    the dimensionlessfree energy potential of the end-to-end distance  -beta * Ln(P(r)*4pi*r^2) = -beta * Ln(histo_collections)
    where beta=k_BT is set 1.0.
    
    requirements: 
    Numpy
    """
    histo_collections = histo_collections / np.sum(histo_collections) # normalization    
    free_energy = -1.0 * (np.log(histo_collections) - np.log(np.sum(histo_collections))) / beta
    return free_energy

def looping_p(histo_collections, bin_edges, Rmax, Rmin):
    """
    looping entropy P_L is defined as the integral from Rmin to Rmax over P(R)*4pi*R^2*dR; since P(R)*4pi*R^2 is equivalent to hitso_collection[i] for 
    the bin between bin_edges[i] and bin_edges[i+1] we have P_L = integral from Rmin to Rmax over P(R)*4pi*R^2*dR ~ sum from Rmin to Rmax over 
    P(bin_center_i)*4pi*bin_centers_i^2*(bin_edges_i+1 - bin_edges_i)=sum from Rmin to Rmax over hitso_collection_i*(bin_edges_i+1 - bin_edges_i)
    
    Since the sizes of the bins are equal ((bin_edges_i+1 - bin_edges_i)=constant), the effect of (bin_edges_i+1 - bin_edges_i) is cancell out upon
    normalization.
    
    Inputs:
    histo_collections: a numpy array of the histograms of the Flory radius or the end-to-end distance during the whole run (the total number of 
    data points is equal to the total number of time steps.)
    bin_edges: the edges of the bins used in histograming.
    Rmax: the minimum center-to-center distance the two end-of-chain monomers can have; Since they have the same size, this is usually equal to their
    diameter (or size). 
    Rmin: this is the size of the crowders that can reside between two monomers. When the center-to-center distance between two monomers is less than 
    Rmax+Rmin, no crowders can be between two monomers and the depletion force is non-zero.
    
    Returns:
    The probability of looping.
    
    requirements: 
    Numpy
    """
       
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 # the centers of the bins being used for plotting
    histo_collections = np.multiply(histo_collections,(bin_edges[1:] - bin_edges[:-1]))
    histo_collections = histo_collections / np.sum(histo_collections) # normalization
    looped_probablity = 0
    for i in range(len(bin_centers)):
        if bin_centers[i] <= (Rmax + Rmin):
            looped_probablity += histo_collections[i]
    return looped_probablity  


def looping_plotter(plotting_df):
    with sns.plotting_context("paper"):
        font_name = "Times New Roman"

        fig, axes = plt.subplots(nrows=2, ncols=1,sharex=True, figsize=(16,9))
        plt.rcParams['font.sans-serif'] = font_name
        plt.rcParams['font.family'] = font_name
        plt.rcParams["font.serif"] = font_name
        plt.rcParams['mathtext.fontset'] = 'dejavuserif'

        colors = sns.color_palette("husl", 3)
        
        nmon = plotting_df.loc[0,'nmon']
        dmon_large = plotting_df.loc[0,'dmon_large']
        dcrowd = plotting_df.loc[0,'dcrowd']

        style = {'marker':'s','markersize':10,'markerfacecolor':'none','linestyle':'-','linewidth':2}
        axes[0].plot('vfrc_crowd','looping_p',data=plotting_df,c=colors[0],**style)
        axes[0].set_title(r'$N={},a_M={},a_c={}$'.format(nmon,dmon_large,dcrowd),fontsize=20,fontname=font_name)
        axes[0].set_xlim(-0.002, 0.302)
        axes[0].set_ylim(-0.002, 0.16)
        axes[0].set_ylabel(r'$P_L(R_{min},R_{max})$',fontsize=18,fontname=font_name)
        axes[0].text(0.002,0.11,r'$\Delta=[0,a_M+a_m]$', fontsize=18, va='bottom',fontname=font_name)
        axes[0].xaxis.set_major_locator(MultipleLocator(0.05))
        axes[0].xaxis.set_minor_locator(MultipleLocator(0.025))
        axes[0].yaxis.set_major_locator(MultipleLocator(0.05))
        axes[0].yaxis.set_minor_locator(MultipleLocator(0.025))
        axes[0].grid(which='major',ls=':',lw=1)
        axes[0].tick_params(axis='both', which='major', labelsize=16)
        axes[0].tick_params(axis='both', which='minor', labelsize=12)


 
        axes[1].set_xlim(-0.002, 0.302)
        axes[1].set_ylim(0.7, 1.05)
        axes[1].plot('vfrc_crowd','rflory_nor',data=plotting_df,c=colors[2],**style)
        axes[1].set_xlabel(r'$\phi_c$',fontsize=20,fontname=font_name)
        axes[1].set_ylabel(r'$R_F/R_0$',fontsize=18,fontname=font_name)
        axes[1].xaxis.set_major_locator(MultipleLocator(0.05))
        axes[1].xaxis.set_minor_locator(MultipleLocator(0.025))
        axes[1].yaxis.set_major_locator(MultipleLocator(0.1))
        axes[1].yaxis.set_minor_locator(MultipleLocator(0.05))
        axes[1].grid(which='major',ls=':',lw=1)
        axes[1].tick_params(axis='both', which='major', labelsize=16)
        axes[1].tick_params(axis='both', which='minor', labelsize=12)
        fig.align_ylabels()
        today = datetime.date.today().strftime('%Y%m%d')
        plt.savefig(today+'_loopP_and_R_vs_vfrc_c.png',dpi=300,format='png')   

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
    
    Issue:
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
    self.radius_type: raduis_type
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
    
    def __init__(self, histogram, properties, radius_type, geometry, direction):
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
        self.radius_type = radius_type
        self.r_particle = 0.5 * self.properties[self.radius_type].values[0] # Assumed the sizes of the particle of interest are the same for all the columns.
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
        V (float) Volume of the spherical segment

        Requierment: 
        Numpy.
        """
        left = min(a,b) # Distance of the left base (lower bond) from the center of the sphere
        right = max(a,b) # Distance of the right base (upper bond) from the center of the sphere
        # always right>left since we assume the axis over which left, center, right are located points to right.
        V = np.pi * (r**2*(right-left) - (right**3 - left**3)/3)
        if np.abs(right) > r or np.abs(left) > r:
            #raise ValueError(" a and b should be within [-r,r].")
            #print(" a and b should be within [-r,r], otherwise the segment volume is 0.")
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
                if (rightmost[idx] >= self.edges[edge_idx]) and (rightmost[idx] < self.edges[edge_idx+1]): # the index of the rightmost bin (or the index of the **left edge** of rigthmost bin) is set as the index of the bin in which the rightmost side of the bead is located.
                    rightmost_idx[idx] = edge_idx
        self.particle_bounds = np.column_stack((leftmost_idx,rightmost_idx)) 

    def _consecutive_vol_shares(self):
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

        
        """
        box_length = self.edges[-1] - self.edges[0]
        self.volume_shares = {}
        for center_idx, bounds_minmax in enumerate(self.particle_bounds):
            self.volume_shares[center_idx] = {}
            # share of the middle bins from the volume of the sphere:
            # edge value is an increasing function of edge index
            
            if (bounds_minmax[0] > bounds_minmax[1]) and (center_idx <= len(self.centers)//2): 
                for edge_idx in range(bounds_minmax[0],len(self.centers),1):
                    center = self.centers[center_idx] + box_length # Due to the PBC, the bin center is moved to right side, so the volume of intgersection can be calculated correctly.
                    left_distance = self.edges[edge_idx] - center
                    if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                        left_distance = -1 * self.r_particle # -1 of r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the left side so we use -1*r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    right_distance =  self.edges[edge_idx+1] - center # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    if np.abs(right_distance) >= self.r_particle: 
                        right_distance = self.r_particle  # r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the right side so we use r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
                
                for edge_idx in range(0, bounds_minmax[1]+1,1):
                    left_distance = self.edges[edge_idx] - self.centers[center_idx]
                    if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                        left_distance = -1 * self.r_particle # -1 of r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the left side so we use -1*r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    right_distance =  self.edges[edge_idx+1] - self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    if np.abs(right_distance) >= self.r_particle: 
                        right_distance = self.r_particle  # r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the right side so we use r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)

            elif (bounds_minmax[0] > bounds_minmax[1]) and (center_idx > len(self.centers)//2):
                for edge_idx in range(bounds_minmax[0],len(self.centers),1):
                    left_distance = self.edges[edge_idx]-self.centers[center_idx]
                    if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                        left_distance = -1 * self.r_particle # -1 of r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the left side so we use -1*r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    right_distance =  self.edges[edge_idx+1]-self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    if np.abs(right_distance) >= self.r_particle: 
                        right_distance = self.r_particle  # r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the right side so we use r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
                
                for edge_idx in range(0, bounds_minmax[1]+1,1):
                    center = self.centers[center_idx] - box_length # Due to the PBC, the bin center is moved to left side, so the volume of intgersection can be calculated correctly.
                    left_distance = self.edges[edge_idx] - center
                    if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                        left_distance = -1 * self.r_particle # -1 of r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the left side so we use -1*r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    right_distance =  self.edges[edge_idx+1] - center # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    if np.abs(right_distance) >= self.r_particle: 
                        right_distance = self.r_particle  # r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the right side so we use r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    self.volume_shares[center_idx][edge_idx] = self.spherical_segment(self.r_particle, left_distance, right_distance)
            
            else:
                for edge_idx in range(bounds_minmax[0],bounds_minmax[1]+1,1):
                    left_distance = self.edges[edge_idx]-self.centers[center_idx]
                    if np.abs(left_distance) >= self.r_particle: # the most left bound can be a spherical cap or a spherical segment; the spherical segments are used for the bins near the bounds of simulation box.
                        left_distance = -1 * self.r_particle # -1 of r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the left side so we use -1*r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
                    right_distance =  self.edges[edge_idx+1]-self.centers[center_idx] # the most right bound can be a spherical cap or shperical segment; the spherical segments are used for the bins near the bounds of simulation box.
                    if np.abs(right_distance) >= self.r_particle: 
                        right_distance = self.r_particle  # r_particle is used since we want to the volume of the smaller sphiercal cap (a sphere can be divide to two spherical caps where the sum of their widths is equal to the radius of the sphere; here we want the amller one and it is always at the right side so we use r_particle to find its volume - This works well with the algorithm we define to calculate the volume of a spherical cap or segment.
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
                "polar": (chosen_properties['dcyl'].values[0], ),
                 "azimuthal": (chosen_properties['dcyl'].values[0], ),
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
        self.phi[col_name] = phi / np.sum(phi)

    def _run(self, to_file=True):
        """
        _run perform a list of operation over the columns in the input histogram.
        """
        for col_name in self.histogram.columns:
            self._set_args(col_name)
            self._number_density(col_name)
            self._volume_fraction(col_name)      

def distributions_generator(histograms, properties, radius_type, geometry, direction, particle_name, to_file=True):
    """
    distributions_generator generate the local number density (rho) and volume fraction (phi) of the particle_name woth particle_type column name in properties dataframe along the direction of interest in the geometry of interest.
    
    Caution: 
    For the sumrule problem in the cylindrical goemetry, a simulation group is a collection of simulations that all have the same values for the number of monomers, the diameter of cylinder, and the size of crowders (assuming size of monomers is 1). An ensemble (usually with M number of simulations) is a collection of themodynamically-equivalent simulations that all have the same values for the number of monomers, the diameter of cylinder, the size of crowders, and the same number of crowders (assuming size of monomers is 1).  In standard statitatical mechanical approach, an ensmeble is equivalent a simulation, but here we used it to reffer to all the thermodynamically-equivalent simulations.
    In an ensemble-average simulation group, each ensemble is replaced with the average of all its simulation, so if we have M simulation in an ensemble and P ensembles in a group, then we have M*P simulations in a simulation group and P ensemble-averaged simulations in a ensemble-averaged simulation group.
    
    Parameters:
    histogram (dict): a dictionary of an ensemble or ensemble-averaged or group simulation in which the keys are the name of ensembles/ensemble-averaged groups/groups and the keys of the histogram dataframes. The names of the columns in each dataframe are the name of simulations in an ensemble or the name of ensemble-averaged group or the name of simulation group, the columns are the frequenies of partilces of the type of interest (see radius_type) in each of the bins, and the number of columns is the number of simulations in a ensemble, one in an ensemble-averaged group, or the number of ensmebles (=the number of ensemble-averaged) groups an a simulations group. 
    properties (pandas dataframe): the properties of each simulation/ensemble-averaged simulation/ensemble-averaged simulations of the simualtions/ensemble-averaged simulation/ensemble-averaged simulations in a ensemble/ensemble-averged group/simulation group.
    raduis_type (str): the name of the column in the properties in which the size (or diameter) of the particles are tabled. The particle type is the same for all the particles that their frequencies are counted in histogram. 
    geometry (str): the shape of the simulation box
    direction (str): the direction of interest in the geometry of interest.
    particle_name: the name of paticle type.
    to_file: whether save to file or not
    
    """
    densities = {} 
    vol_fractions = {}
    for bunch_name, bunch_histogram in histograms.items(): # a bunch of simulations in a group/ensemble/ensemble-averaged group
        distributions = Distributions(bunch_histogram, properties, radius_type, geometry, direction)
        densities[bunch_name] = distributions.rho            
        vol_fractions[bunch_name] = distributions.phi
        if to_file:
            densities[bunch_name].to_csv(bunch_name+'-'+distributions.short_name_rho+particle_name+'.csv')
            vol_fractions[bunch_name].to_csv(bunch_name+'-'+distributions.short_name_phi+particle_name+'.csv')    
    return densities, vol_fractions

def yticks(axis,limits,code=False,decimals=3,**kwargs):
    ylow,yhigh,ystep,ylim = limits
    if code == False:
        axis.set_ylim(ylow-ylim,yhigh+ylim)
        axis.set_yticks([0 if (i < ylim) and (i > -1*ylim) else i for i in np.arange(ylow,yhigh+ystep,ystep)])
        axis.set_yticklabels(labels=[])
        #axis.set_yticklabels(labels=np.around(axis.get_yticks(), decimals=2),fontsize=16)
    else: 
        axis.set_ylim(ylow-ylim,yhigh+ylim)
        axis.set_yticks([0 if (i < ylim) and (i > -1*ylim) else i for i in np.arange(ylow,yhigh+ystep,ystep)])
        axis.set_yticklabels(labels=np.around(axis.get_yticks(), decimals=decimals),**kwargs)

def xticks(axis,limits,code=False,decimals=3,**kwargs):
    xlow,xhigh,xstep,xlim = limits
    if code == False:
        axis.set_xlim(xlow-xlim,xhigh+xlim)
        axis.set_xticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.arange(xlow,xhigh+xstep,xstep)])
        axis.set_xticklabels(labels=[])
        
    else: 
        axis.set_xlim(xlow-xlim,xhigh+xlim)
        axis.set_xticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.arange(xlow,xhigh+xstep,xstep)])
        axis.set_xticklabels(labels=np.around(axis.get_xticks(), decimals=decimals),**kwargs)

def yticks_nstep(axis,ylow,yhigh,nstep,ylim,code=False,decimals=3,**kwargs):
    if code == False:
        axis.set_ylim(ylow-ylim,yhigh+ylim)
        axis.set_yticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.linspace(ylow,yhigh,nstep)])
        axis.set_yticklabels(labels=[])
        #axis.set_yticklabels(labels=np.around(axis.get_yticks(), decimals=2),fontsize=16)
    else: 
        axis.set_ylim(ylow-ylim,yhigh+ylim)
        axis.set_yticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.linspace(ylow,yhigh,nstep)])
        axis.set_yticklabels(labels=np.around(axis.get_yticks(), decimals=decimals),**kwargs)

def xticks_nstep(axis,xlow,xhigh,nstep,xlim,code=False,decimals=3,**kwargs):
    if code == False:
        axis.set_xlim(xlow-xlim,xhigh+xlim)
        axis.set_xticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.linspace(xlow,xhigh,nstep)])
        axis.set_xticklabels(labels=[])
        
    else: 
        axis.set_xlim(xlow-xlim,xhigh+xlim)
        axis.set_xticks([0 if (i < xlim) and (i > -1*xlim) else i for i in np.linspace(xlow,xhigh,nstep)])
        axis.set_xticklabels(labels=np.around(axis.get_xticks(), decimals=decimals),**kwargs)
        
def change_legend_name(line, legend_new_labels):
    legend_old_labels = line.legend(fontsize=16,bbox_to_anchor=(1.005, 1),
                                    loc=2,edgecolor='black',title_fontsize=16,markerscale=1.5).texts
    for new_label, line in zip(legend_old_labels, legend_new_labels): new_label.set_text(line)
        

def chainsize_plot(df, xcol, leg_labels, colors, fontsize=20):
    sns.set_context('paper')
    sns.set_style("ticks")
    fig, axes = plt.subplots(nrows=3,ncols=1,sharex=True,figsize=(16,12))
    line1 = sns.lineplot(x=xcol, y="fsd_normalized", hue='dcyl',style='dcrowd', size='nmon', palette=colors, markers=True, markersize=8, data=df,ax=axes[0])
    line2 = sns.lineplot(x=xcol, y="gyr_normalized", hue='dcyl',style='dcrowd', size='nmon', palette=colors, markers=True, markersize=8, data=df, legend=False,ax=axes[1])
    line3 = sns.lineplot(x=xcol, y="rflory_normalized", hue='dcyl',style='dcrowd', size='nmon', palette=colors, markers=True, markersize=8, data=df, legend=False,ax=axes[2])

    xlabels = {"phi_c":r"$\phi_c$","phi_c_normalized":r"${a\phi_c}/{a_c}$","phi_c_eff":r"$\phi_c^{eff}$",
               "phi_c_eff_normalized":r"${a\phi_c^{eff}}/{a_c}$"}
    
    ylabels = {"phi_c":[r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
               "phi_c_eff":[r'$\frac{L_{FSD}(\phi_c)}{L_{FSD}(0)}$',r'$\frac{R_{ROG}(\phi_c)}{R_{ROG}(0)}$',r'$\frac{R_{Flory}(\phi_c)}{R_{Flory}(0)}$'],
              "phi_c_normalized":[r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$'],
              "phi_c_eff_normalized":[r'$\frac{L_{FSD}({a\phi_c}/{a_c})}{L_{FSD}(0)}$',r'$\frac{R_{ROG}({a\phi_c}/{a_c})}{R_{ROG}(0)}$',r'$\frac{R_{Flory}({a\phi_c}/{a_c})}{R_{Flory}(0)}$']}

    for num, axis in enumerate(axes):
        axis.grid(True,ls=':',lw=1)
        axis.tick_params(axis ='both',direction='inout',width=1)
        axis.set_ylabel(ylabels[xcol][num],fontsize=fontsize)
        yticks(axis,(0,1.0,0.2,0.04),code=True,fontsize=14)

    xticks_dict = {"phi_c":(0.0,0.4,0.05,0.005),
                   "phi_c_normalized":(0.0,0.4,0.05,0.005),
                   "phi_c_eff":(0.0,0.65,0.05,0.005),
                   "phi_c_eff_normalized":(0.0,0.5,0.05,0.005)}
    
    xticks(axes[2],xticks_dict[xcol],code=True,fontsize=14)
    change_legend_name(line1,leg_labels)
    line3.set_xlabel(xlabels[xcol],fontsize=fontsize)

    picname = "chainsize-"+xcol
    plt.savefig(picname+'.pdf',dpi=300,bbox_inches='tight')
    

def organize_hists(csv_files, xname, yname, divider):
    """
    This funcions receives a list of paths to histogram data (each contains two columns; one bins, other counts) in different directions and return a dict of paris of filenames and Pandas dataframes as pairs of key and values.
    """
    dfs = {}
    for idx, csvfile in enumerate(csv_files):
        dummdy_df = pd.read_csv(csvfile[0],index_col=0)
        rphic_sum = dummdy_df.iloc[:,0].sum()
        keyname = list(dummdy_df.columns)[0].split('-')[0]
        cell_attrs = cellAttributes(keyname,warning=False)
        divider_dicts = {'dcyl':cell_attrs.dcyl, 'lcyl':cell_attrs.lcyl, 'degree': 360}
        dummdy_df = pd.read_csv(csvfile[0],names=[xname,yname],skiprows=1)
        dummdy_df[xname+'_norm'] = 2 * dummdy_df[xname] / divider_dicts[divider]
        dummdy_df[yname+'_norm'] = 2 * dummdy_df[yname] / rphic_sum
        dfs[keyname] = dummdy_df
    return dfs

def ls_handler(attributes, linstyles):
    """
    This function creates handles for matplotlib legend functions based on the line styles. It recieves a list  physical attributes and a list of linestyles (the two list have the same size) and returns a list of handles. 
    """
    handles = []
    for attr, linestyle in zip(attributes,linstyles):
        mline= mlines.Line2D([], [],  color='black', ls=linestyle, label=attr, lw=2)
        handles.append(mline)
    return handles

def lw_handler(attributes, linewidths):
    """
    This function creates handles for matplotlib legend functions based on the line widths. It recieves a list  physical attributes and a list of linewidths (the two list have the same size) and returns a list of handles. 
    """
    handles = []
    for attr, linewidth in zip(attributes,linewidths):
        mline= mlines.Line2D([], [],  color='black', lw=linewidth, label=attr)
        handles.append(mline)
    return handles

def marker_handler(attributes, markers, **kwargs):
    """
    This function creates handles for matplotlib legend functions based on the line widths. It recieves a list  physical attributes and a list of linewidths (the two list have the same size) and returns a list of handles. 
    """
    handles = []
    for attr, marker in zip(attributes,markers):
        mline= mlines.Line2D([], [],  color='black', marker=marker, label=attr,**kwargs)
        handles.append(mline)
    return handles

def all_in_one_properties(ens_evg_properties_files,df_name,geometry,filename=None, round_to=4,rename_cols=True):
    ens_avg_properties_df = df_of_properties(ens_evg_properties_files, df_name, geometry, to_file=False,index_col=0)
    if rename_cols:
            ens_avg_properties_df.rename(columns={"vfrc_m": "phi_m", "vfrc_crowd": "phi_c",'rho_crowd':'rho_c'},inplace=True)
    ens_avg_properties_df = ens_avg_properties_df.round(round_to)
    ens_avg_properties_df['fsd_normalized'] = 0.0
    ens_avg_properties_df['gyr_normalized'] = 0.0
    ens_avg_properties_df['rflory_normalized'] = 0.0
    ens_avg_properties_df['phi_c_normalized'] = (ens_avg_properties_df['dmon'] * ens_avg_properties_df['phi_c']) / ens_avg_properties_df['dcrowd'] 
    unique_simulations = [list(input_set) for input_set in list(ens_avg_properties_df.groupby(['nmon','dcyl','dcrowd']).size().index)]
    for nmon, dcyl, dcrowd in unique_simulations:
        condition = (ens_avg_properties_df['nmon'] == nmon) & (ens_avg_properties_df['dcyl'] == dcyl) & (ens_avg_properties_df['dcrowd'] == dcrowd) & (ens_avg_properties_df['phi_c'] == 0)
        fsd = ens_avg_properties_df[condition]['fsd'].values[0]
        gyr = ens_avg_properties_df[condition]['gyr'].values[0]
        rflory = ens_avg_properties_df[condition]['rflory'].values[0]
        cond_normal = (ens_avg_properties_df['nmon'] == nmon) & (ens_avg_properties_df['dcyl'] == dcyl) & (ens_avg_properties_df['dcrowd'] == dcrowd)
        ens_avg_properties_df.loc[cond_normal,'fsd_normalized'] = ens_avg_properties_df.loc[cond_normal,'fsd'] / fsd
        ens_avg_properties_df.loc[cond_normal,'gyr_normalized'] = ens_avg_properties_df.loc[cond_normal,'gyr'] / gyr
        ens_avg_properties_df.loc[cond_normal,'rflory_normalized'] = ens_avg_properties_df.loc[cond_normal,'rflory'] / rflory

    ens_avg_properties_df['phi_c_eff'] = (np.pi * ens_avg_properties_df['dcrowd'] ** 3 / 6) * ens_avg_properties_df['ncrowd'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dcrowd']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['phi_c_eff_normalized'] = (ens_avg_properties_df['dmon'] * ens_avg_properties_df['phi_c_eff']) / ens_avg_properties_df['dcrowd'] 
    ens_avg_properties_df['rho_c_normalized'] = ens_avg_properties_df['rho_c'] * ens_avg_properties_df['dcrowd'] ** 2
    
    ens_avg_properties_df['rho_c_eff'] = ens_avg_properties_df['ncrowd'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dcrowd']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['rho_c_eff_normalized'] = ens_avg_properties_df['rho_c_eff'] * ens_avg_properties_df['dcrowd'] ** 2
    
    ens_avg_properties_df['phi_m_eff'] = (np.pi * ens_avg_properties_df['dmon'] ** 3 / 6) * ens_avg_properties_df['nmon'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dmon']) ** 2) * ens_avg_properties_df['lcyl'])
    ens_avg_properties_df['rho_m_eff'] = ens_avg_properties_df['nmon'] / ((np.pi / 4 * (ens_avg_properties_df['dcyl']-ens_avg_properties_df['dmon']) ** 2) * ens_avg_properties_df['lcyl'])
    
    
    ens_avg_properties_df = ens_avg_properties_df.round(round_to)
    if filename != None:
        ens_avg_properties_df.to_csv(filename+'.csv',index=False)
    return ens_avg_properties_df

def organize_hists_ens_evg(pair_csvs,properties,xname,yname,xnorm,ynorm,filename=None,round_to=4):
    """
    This funcions receives a list of paths to histogram data (each contains two columns; one bins, other counts) in different directions and return a dict of paris of filenames and Pandas dataframes as pairs of key and values.
    """
    group_dfs = []
    properties[['groupname','garbage']] = properties.filename.str.split(pat='-',expand=True)
    for crd_csv, mon_csv in pair_csvs:
        
        crd_df = pd.read_csv(crd_csv,index_col=0)
        group_name = list(crd_df.columns)[0].split('-')[0]
            
        crd_df = pd.read_csv(crd_csv,names=[xname,yname+'_c'],skiprows=1)
        mon_df = pd.read_csv(mon_csv,names=[yname+'_m'],skiprows=1,usecols=[1])
        merged_df = pd.concat([crd_df,mon_df],axis=1)
        
        
        for col in properties.columns:
            cond = properties['groupname'] == group_name
            merged_df[col] = properties[cond][col].values[0]
            
        merged_df[xname+'_norm'] = 2 * merged_df[xname] / merged_df[xnorm]
        merged_df['sum_rule_rho'] = merged_df['rho_c'] * merged_df['dcrowd'] ** 2
        merged_df['sum_rule_rho_norm'] = merged_df['sum_rule_rho']/merged_df['sum_rule_rho'].sum()
        
        merged_df['sum_rule_phi'] = merged_df['phi_c'] / merged_df['dcrowd']
        merged_df['sum_rule_phi_norm'] = merged_df['sum_rule_phi']/merged_df['sum_rule_phi'].sum()
        
        merged_df['sum_rule_rho_eff'] = merged_df['rho_c_eff'] * merged_df['dcrowd'] ** 2
        merged_df['sum_rule_rho_eff_norm'] = merged_df['sum_rule_rho_eff']/merged_df['sum_rule_rho_eff'].sum()
        merged_df['sum_rule_phi_eff'] = merged_df['phi_c_eff'] / merged_df['dcrowd']
        merged_df['sum_rule_phi_eff_norm'] = merged_df['sum_rule_phi_eff']/merged_df['sum_rule_phi_eff'].sum()
        
        if ynorm != 'hist' : 
            merged_df[yname+'_norm'+'_c'] = merged_df[yname+'_c'] / merged_df[ynorm+'_c']
            merged_df[yname+'_norm'+'_c'] = merged_df[yname+'_norm'+'_c'] / merged_df[yname+'_norm'+'_c'].sum()
            merged_df[yname+'_norm'+'_m'] = merged_df[yname+'_m'] / merged_df[ynorm+'_m']
            merged_df[yname+'_norm'+'_m'] = merged_df[yname+'_norm'+'_m'] / merged_df[yname+'_norm'+'_m'].sum()
            
            merged_df[yname+'_sum'] = (merged_df[yname+'_c'] + merged_df[yname+'_m']) / (merged_df[ynorm+'_c'] + merged_df[ynorm+'_m'])
            merged_df[yname+'_norm_sum'] = merged_df[yname+'_sum'] / merged_df[yname+'_sum'].sum()
            
            merged_df[yname+'_sum_eff'] = (merged_df[yname+'_c'] + merged_df[yname+'_m']) / (merged_df[ynorm+'_c_eff'] + merged_df[ynorm+'_m'])
            merged_df[yname+'_norm_sum_eff'] = merged_df[yname+'_sum_eff'] / merged_df[yname+'_sum_eff'].sum()
    
        
        merged_df.drop(['filename','garbage'],axis=1,inplace=True)
        group_dfs.append(merged_df)
    
    df_of_hists = pd.concat(group_dfs)
    df_of_hists = df_of_hists.round(round_to)
    df_of_hists.reset_index(inplace=True,drop=True)
    if filename != None:
        outname = filename+'.csv'
        df_of_hists.to_csv(outname)
        
    return df_of_hists