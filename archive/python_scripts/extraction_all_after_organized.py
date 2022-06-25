# This script extract different bug's information from pairs (toplogy and trajectory) of bug simulation files in oen or more organized *trjs_bug* directories.
from pathlib import Path
import os
from glob import glob
from PipeLine import *
from dask.distributed import Client
from dask import delayed
from dask import compute

cores = ncores = int(os.environ.get('SLURM_CPUS_PER_TASK',default=1))
print(f"number of workers set to {cores}; is this the same requested cores on the cluster?")
#client = Client(n_workers=cores)
home = str(Path.home())
cwdir = str(Path.cwd())
sim_all_dirs = glob(home+'/amirhsi_rrg/cylinder_simulations/N*-trjs_all/N*/')
geom = 'cylindrical'

trjs_computed = []
all_extraction_dir = '/extraction_all/'
for sim_all_dir in sim_all_dirs:
    sim_name = sim_all_dir.split("/")[-2]
    print(sim_name)
    all_trjs = glob(sim_all_dir+"N*.lammpstrj")
    all_trjs = PipeLine.file_reader(all_trjs,extensions=['lammpstrj'])
    all_trjs = [all_trj[0] for all_trj in all_trjs]

    all_topology = glob(sim_all_dir+"N*.all.data")
    all_topology = PipeLine.file_reader(all_topology,extensions=['all.data'])
    all_topology = all_topology[0][0]
    sim_extract_dir = cwdir+all_extraction_dir+sim_name
    Path(sim_extract_dir).mkdir(parents=True, exist_ok=False)
    sim_save_to = sim_extract_dir+"/"
    for all_trj in all_trjs:
        print(all_trj)
        trj_delayed = delayed(PipeLine.extract_trj_all)(all_topology, all_trj, geom, sim_save_to)
        trjs_computed.append(trj_delayed)

results = compute(trjs_computed)

