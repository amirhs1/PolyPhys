from glob import glob
from PipeLine import *
from dask.distributed import Client
from dask import delayed
from dask import compute
if __name__ == '__main__':
    client = Client(n_workers=4, processes=False)
    #client
    
database = '/Users/amirhsi_mini/' # parent path
simulation_type = 'all' # bug or all
input_db_name = "analysis" # or extraction
input_segments_name = "whole_simulations" # whole_simulations or
input_sim_group_name = ("-").join([simulation_type,input_db_name,input_segments_name])
input_sim_groups = glob(database+input_db_name+'/N*-'+input_sim_group_name+'/')
output_db_name = "analysis"
geometry = 'cylindrical'

analysis_delayed = [] 
for input_sim_group in input_sim_groups:
    analysis = delayed(PipeLine.whole_group_analysis_all)(input_sim_group, input_db_name, output_db_name, simulation_type, geometry)
    analysis_delayed.append(analysis)

results = compute(analysis_delayed)