#!/bin/bash

##SBATCH --nodes=1 --ntasks-per-node=32  	# number of nodes and number of task per node for MPI job

##In this case you could also say --mem-per-cpu=3G. 
##The advantage of --mem=45G is that the memory consumed by each individual process doesn't matter, 
##as long as all of them together don’t use more than 45GB. With --mem-per-cpu=3G, the job 
##will be canceled if any of the processes exceeds 3GB.

##SBATCH --mem=65G
#SBATCH --mem-per-cpu=2G                	# Default unit is megabytes. It should be less than 4G
#SBATCH --time=00-00:30                 	# Total time needed for this job. --time DD-HH:MM or -t 0-01:00;  Check the manual of sbatch
#SBATCH --ntasks=1
#SBATCH --output=%J.log                 	# output .log file
#SBATCH --account=def-byha 					# The user/account holder who use this computecanada project. you can use -a instead of --acount
#SBATCH --job-name=JOB_NAME					# Name of the input file.
#SBATCH --output=OUTPUT_%.out 				# Name of the output file(s). Usually, it is JOB_NAME_%.out
#SBATCH --mail-user=ahsadegh@uwaterloo.ca 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.

# Load the module:

module load nixpkgs/16.09  intel/2016.4  openmpi/2.1.1 lammps-omp/20170811

echo "Starting run at: `date`"

CURRENT_DATE="$(date +%Y%m%d)" # Taking the current date in the YYYYMMDD format to use in the log file name
lmp_exec=lmp_icc_openmpi
lmp_input="-v in_filename data.chain.80 -v run_date $CURRENT_DATE -v i 1367 -v r 4.0 -v sig2 0.3 -v n_crowd 1000 -v epsilon1 5.0 -v lz 14.0 -v run_before 50000 -v run_after 100000  -i chain_simple.lmp" 
lmp_output="lammps_lj_output.txt"

${lmp_exec}  ${lmp_input} > ${lmp_output}

echo "Program finished with exit code $? at: `date`"


#----- Description of the submitted task -----
# This script is used to run a lammps MD simulation on a computecanada HPC facility.
# List of variables:
# variable ${run_date} #Current date of the system in this format: YYYYMMDD
# variable ${in_filename}	# the input data file for our MD simulation; Usually, it contains initial configuration o the system.
# variable ${i} 			# i is the seed for random number generator. You may want to run many SAME simulations for ensemble average.
# variable ${r}				# Radius of the cylinderical region of the MD simulation	
# variable ${lz}			# Length of the cylinderical region of the MD simulation
# variable ${n_crowd} 		# Number of crowders
# variable ${sig2} 			# FENE potential constant
# variable ${run_before} 	# Total number of timesteps/runs to reach the user-defined equilibrium.
# variable ${run_after} 	# Total number of timsteps/runs after reaching the user-defined equilibrium
