#!/bin/bash
# Script to run a lammps MD simjulation
#List of variables and their names defined by command line method. Their description and the step they will be used in are described in next sections:
#variable ${in_filename}
#variable ${run_date} #Current date of the system in this format: YYYYMMDD
#variable ${i} # i is the seed for random number generator. You may want to run many SAME simulations for ensemble average.
#variable ${r}
#variable ${lz}
#variable ${n_crowd} # Number of crowders
#variable ${sig2} #
#variable ${run_before} # Total number of timesteps/runs to reach the equilibirum.
#variable ${run_after} # Total number of timsteps/runs after equilibrium

CURRENT_DATE="$(date +%Y%m%d)" # Taking the current date in the YYYYMMDD format to use in the log file name 

mkdir images # Make a directory for saving images
mkdir archive # Make a directory for saving restart files

lmp_serial -v in_filename data.chain.50 -v run_date $CURRENT_DATE -v i 1 -v l 36 -v sig2 0.3 -v n_crowd 1000 -i  *.lmp
