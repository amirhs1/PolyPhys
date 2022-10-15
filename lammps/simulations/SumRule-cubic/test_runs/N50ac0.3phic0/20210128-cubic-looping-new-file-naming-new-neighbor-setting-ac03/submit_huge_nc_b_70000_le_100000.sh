#!/bin/bash

# Test run to estimate the resources needed for the final submissions
# Date:2019 March 13

#SBATCH --ntasks=32 						  	# number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --mem-per-cpu=3000M      						# memory; default unit is megabytes.
#SBATCH --time=02-10:00                 	# Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha 					# The user/account holder who use this computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=ahsadegh@uwaterloo.ca 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.
# Load the module:

module --force purge
module load nixpkgs/16.09   intel/2018.3    openmpi/3.1.4   lammps-omp/20190807

echo "Starting run at: `date`"

lmp_exec=lmp_icc_openmpi
lmp_input="input.lmp"
lmp_output="output.txt"

srun ${lmp_exec} < ${lmp_input} > ${lmp_output}

echo "Program finished with exit code $? at: `date`"
