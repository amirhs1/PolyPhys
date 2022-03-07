from glob import glob
from PipeLine import *

trj_files = glob("./N*.lammpstrj")
all_tuples = PipeLine.file_reader(trj_files,extensions=['lammpstrj'])
all_trjs = [all_tuple[0] for all_tuple in all_tuples]

data_file = glob("./N*.all.data")
all_data = PipeLine.file_reader(data_file,extensions=['all.data'])
all_data = all_data[0][0]

geom = 'cylinder'
    
for all_trj in all_trjs:
    PipeLine.extract_all_mon_trj(all_data, all_trj, geom)