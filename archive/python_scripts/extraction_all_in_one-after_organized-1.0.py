# This script extract different bug's information from pairs (toplogy and trajectory) of bug simulation files in oen or more organized *trjs_bug* directories.

from pathlib import Path
import os
from glob import glob
from PipeLine import *
from dask.distributed import Client
from dask import delayed
from dask import compute

cores = 32
print(f"number of workers set to {cores}; is this the same requested cores on the cluster?")
client = Client(n_workers=cores)


home = str(Path.home())
cwdir = str(Path.cwd())
# information extraction from simulations
geom = 'cylindrical'
bug_pairs = glob(home+'/amirhsi_rrg/cylinder_simulations/N*-trjs_bug/N*bug*')
bug_pairs = PipeLine.file_reader(bug_pairs) # each bug_pair is a pair of trajectory and topopgy file.
trjs_computed = []
for bug_pair in bug_pairs:
    sim_name = bug_pair[0].split("/")[-1].split('bug')[0]
    Path(cwdir+sim_name).mkdir(parents=True, exist_ok=False)
    sim_save_to = cwdir+group_path+"/"
    trj_delayed = delayed(PipeLine.extract_trj_bug)(bug_pair, geom,sim_save_to)
    trjs_computed.append(trj_delayed)
results = compute(trjs_computed)


all_files = glob(home+'/amirhsi_rrg/cylinder_simulations/N*-trjs_all/N*all*')
all_files = PipeLine.file_reader(all_files) # each bug_pair is a pair of trajectory and topopgy file.




trj_files = glob("$N*.all.lammpstrj")
all_tuples = PipeLine.file_reader(trj_files,extensions=['all.lammpstrj'])
all_trjs = [all_tuple[0] for all_tuple in all_tuples]

data_files = glob("./N*.all.data")
all_data = PipeLine.file_reader(data_files,extensions=['all.data'])
all_data = all_data[0][0]
    
for all_trj in all_trjs:
    PipeLine.extract_trj_all(all_data, all_trj, geom)

