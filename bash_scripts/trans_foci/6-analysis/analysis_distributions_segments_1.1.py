from glob import glob
import pandas as pd
from PipeLine import *
from dask.distributed import Client
from dask import delayed
from dask import compute
if __name__ == '__main__':
    client = Client(n_workers=4, processes=False)
    #client
database = '/Users/amirhsi_mini/'
simulation_type = 'all' # bug or all
input_db_name = 'extraction'
input_sim_groups = glob(database+input_db_name+"/N*-"+input_db_name+"/")
output_db_name = "analysis"
all_properties_file = database+output_db_name+"/properties-all_in_one.csv"
all_properties = pd.read_csv(all_properties_file,index_col=0)

# use this if there is no "properties-all_in_one.csv" file
"""
database = '/Users/amirhsi_mini/extraction/'
input_db_name = 'extraction'
csv_files = glob(database+input_db_name+"/N*-extraction/*.all*.csv") # the dot "." with"all" is crutial.
output_db_name = "analysis"
properties_csvs = glob( database+output_db_name+"/N*-bug-analysis/N*-properties.csv")
all_properties = []
for properties_csv in properties_csvs:
    df = pd.read_csv(properties_csv,index_col=0)
    all_properties.append(df) 
all_properties = pd.concat(all_properties)
all_properties.reset_index(inplace=True,drop=True)
"""

geometry = 'cylindrical'
analysis_delayed = []   
for input_sim_group in input_sim_groups: # run over all simulations in all the groups
    analysis = delayed(PipeLine.whole_group_analysis_segments)(input_sim_group, input_db_name, output_db_name, simulation_type, all_properties, geometry)
    analysis_delayed.append(analysis)

results = compute(analysis_delayed) # this takes around 20 minutes for 12 simulation groups with 4 workers.