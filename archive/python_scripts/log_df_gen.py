# Importing necessary packages:
from glob import glob
from PipeLine import *

log_files = glob('./N*.log')
geometry = 'cylinder'
log_files = PipeLine.file_reader(log_files,extensions=['.log'])
details_out , runtime_out = PipeLine.log_outputs(log_files[0][0], geometry=geometry)
PipeLine.lammps_log_details(log_files, details_out , runtime_out)    


