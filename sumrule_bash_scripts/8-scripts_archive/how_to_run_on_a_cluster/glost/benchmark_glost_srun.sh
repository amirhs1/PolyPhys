#!/bin/bash

#SBATCH --ntasks=1 						  	# number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --mem=2500M							# memory for each node 
#SBATCH --time=00-08:00                 	# Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha 					# The user/account holder who use this computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=ahsadegh@uwaterloo.ca 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.

# Load the module:
# Load GLOST module along with the modules required to run your application:

module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps-omp/20170811 glost/0.3.1

echo "Starting run at: `date`"

lmp_exec=lmp_icc_openmpi

# Run GLOST with the argument: list_glost_tasks.txt

srun glost_launch glost_input.txt

echo "Program glost_launch finished with exit code $? at: `date`"
echo "Program finished with exit code $? at: `date`"
