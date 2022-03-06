import datetime
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.coordinates.LAMMPS import DumpReader
from MDAnalysis import transformations as mdatransform
from MDAnalysis.analysis import (diffusionmap, align, rms)

from polyphys.manage.parser import SumRule

def test_log_outputs(log_file, geometry):
    pass


def test_lammps_log_details(log_files, details_out , runtime_out):
    pass


def test_cylinder_write(output,cellAttrs):
    pass


def test_chain_stats(array, output):
    pass

def test_error_calc_block(data, filename):
    pass

def test_end_to_end(r):
    pass

def test_max_distance(r):
    pass

def test_fsd(coordinates, axis=2):
    pass

def test_bin_ceate(sim_name, bin_size, lmin, lmax, edge_name, save_to):
    pass

class AmirDumpReader(DumpReader):
    pass



def test_extract_trj_bug(simulation_pair, geometry, save_to="./"):
    pass
def test_rmsd_trj_bug(fname, geometry, save_to):
    pass


def test_extract_trj_all(all_data, all_trj, geometry, save_to):
    pass
