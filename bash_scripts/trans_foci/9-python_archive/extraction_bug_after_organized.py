# This script extract different bug's information from pairs (toplogy and trajectory) of bug simulation files in oen or more organized *trjs_bug* directories.

from glob import glob
from pathlib import Path
from PipeLine import *
import os
from dask.distributed import Client
from dask import delayed
from dask import compute
ncores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
client = Client(n_workers=ncores)

home = str(Path.home())
cwdir = str(Path.cwd())
geom = 'cylindrical'
bug_pairs = glob(home+'/amirhsi_rrg/cylinder_simulations/N*-trjs_bug/N*bug*')

bug_pairs = PipeLine.file_reader(bug_pairs) # each bug_pair is a pair of trajectory and topopgy file.
trjs_computed = []
bug_dir = '/extraction_bug/'
for bug_pair in bug_pairs:
    sim_name = bug_pair[0].split("/")[-1].split('.bug')[0]
    sim_dir = cwdir+bug_dir+sim_name
    print(bug_pair)
    Path(sim_dir).mkdir(parents=True, exist_ok=False)
    sim_save_to = sim_dir+"/"
    trj_delayed = delayed(PipeLine.extract_trj_bug)(bug_pair, geom, sim_save_to)
    trjs_computed.append(trj_delayed)
    
results = compute(trjs_computed)