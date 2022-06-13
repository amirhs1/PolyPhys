from glob import glob
from PipeLine import *

fname = glob("./N*.bug.*")
fname = PipeLine.file_reader(fname) # This is a list with one member

geom = 'cylindrical'
print(fname)
PipeLine.extract_trj_bug(fname[0], geom) # A list with one member, the member is a tuple of a trj and data pair.
PipeLine.bug_trj_rmsd(fname[0], geom)

trj_files = glob("./N*.lammpstrj")
all_tuples = PipeLine.file_reader(trj_files,extensions=['lammpstrj'])
all_trjs = [all_tuple[0] for all_tuple in all_tuples]

data_file = glob("./N*.all.data")
all_data = PipeLine.file_reader(data_file,extensions=['all.data'])
all_data = all_data[0][0]

    
for all_trj in all_trjs:
    PipeLine.extract_all_trj(all_data, all_trj, geom)

