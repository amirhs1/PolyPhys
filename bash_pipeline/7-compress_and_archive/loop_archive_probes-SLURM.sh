#!/bin/bash

#SBATCH --ntasks=1 						  	# number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --mem-per-cpu=3000M      						# memory; default unit is megabytes.
#SBATCH --time=00-01:00                	# Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha 					# The user/account holder who use th computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.
#SBATCH --job-name="archive_job"

> archive_report.txt  # Clear previous report contents
echo "Starting run at: $(date)"
for dir in N*/; do
    tar -zcvf "${dir::-1}".tar.gz "${dir}" 
    echo "Finished archiving ${dir}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"

