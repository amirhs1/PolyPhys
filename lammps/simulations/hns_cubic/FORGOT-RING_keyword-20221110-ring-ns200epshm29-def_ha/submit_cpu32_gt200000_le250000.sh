#!/bin/bash

#SBATCH --nodes=1 					# number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --ntasks-per-node=32
#SBATCH --mem=96G                	# memory; default unit is megabytes.
#SBATCH --time=05-00:00             # Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha 		    # The user/account holder who use this computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.

# Load the module:

module --force purge
module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3   lammps-omp/20210929

echo "Starting run at: $(date)"

lmp_exec=lmp
lmp_input="input.lmp"
lmp_output="output.txt"

srun ${lmp_exec} < ${lmp_input} > ${lmp_output}

echo "Program finished with exit code $? at: $(date)"