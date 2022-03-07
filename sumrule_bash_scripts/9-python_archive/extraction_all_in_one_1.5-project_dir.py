from glob import glob
from PipeLine import *
from pathlib import Path

home = str(Path.home())

fname = glob(home+'/amirhsi_rrg/cylinder_simulations/N*-trjs_bug/N*bug*')

geom = 'cylindrical'
fname = glob("./N*.bug.*")
fname = PipeLine.file_reader(fname) # This is a list with one member
PipeLine.extract_trj_bug(fname[0], geom) # A list with one member, the member is a tuple of a trj and data pair.
PipeLine.rmsd_trj_bug(fname[0], geom)

trj_files = glob("$N*.all.lammpstrj")
all_tuples = PipeLine.file_reader(trj_files,extensions=['all.lammpstrj'])
all_trjs = [all_tuple[0] for all_tuple in all_tuples]

data_files = glob("./N*.all.data")
all_data = PipeLine.file_reader(data_files,extensions=['all.data'])
all_data = all_data[0][0]
    
for all_trj in all_trjs:
    PipeLine.extract_trj_all(all_data, all_trj, geom)

