from glob import glob
from PipeLine import *
from dask.distributed import Client
from dask import delayed
from dask import compute
if __name__ == '__main__':
    client = Client(n_workers=4, processes=False)
    #client

# Bugs
# This notebook cannot convert to a python file since I do not know how to efficently use dask in a *.py* file. I try to solve this issue by configuring ask apprportiately there. Check the file.
# I get these errors:
#distributed.utils_perf - WARNING - full garbage collections took 19% CPU time recently (threshold: 10%)
#distributed.worker - WARNING - Memory use is high but worker has no data to store to disk.  Perhaps some other process is leaking memory?  Process memory: 6.47 GiB -- Worker memory limit: 4.00 GiB

# Path to bug extraction folders/files
# analyze all the bug groups at once:
database = '/Users/amirhsi_mini/' # parent path
simulation_type = 'bug' # bug or all
input_db_name = "extraction"
input_sim_groups = glob(database+input_db_name+"/N*-extraction/") # Path to bug extraction folders/files
output_db_name = "analysis"
geometry = 'cylindrical'
analysis_delayed = []
for input_sim_group in input_sim_groups:
    analysis = delayed(PipeLine.whole_group_analysis_bug)(input_sim_group, input_db_name, output_db_name, simulation_type, geometry)
    analysis_delayed.append(analysis)

# it takes less than 5 min for 12 simulation groups with 4 workers.
results = compute(analysis_delayed)

# create one dataframe of all the properties files:
database_path = database+output_db_name+"/"
properties_files = glob(database_path+"N*-"+simulation_type+"-"+output_db_name+"/N*-properties.csv")
properties_files = PipeLine.file_reader(properties_files,extensions=['-properties.csv'])
properties_all_in_one = PipeLine.all_in_one_properties(properties_files, save_to=database_path, round_to=4, index_col=0)
properties_files_ens_avg = glob(database_path+"N*-"+simulation_type+"-"+output_db_name+"-ens_avg/N*-properties-ens_avg.csv")
properties_files_ens_avg = PipeLine.file_reader(properties_files_ens_avg,extensions=['-properties-ens_avg.csv'])
properties_all_in_one_ens_avg = PipeLine.all_in_one_properties(properties_files_ens_avg, ens_avg=True, norm_func=PipeLine.cylindrical_norm, save_to=database_path, round_to=4, index_col=0)