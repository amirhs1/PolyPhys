#!/bin/bash

#SBATCH --ntasks=1                                                      # number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --mem-per-cpu=3000M                                                     # memory; default unit is megabytes.
#SBATCH --time=00-05:00                 # Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha                                      # The user/account holder who use th computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com   # The mail address for notifications about submitted job
#SBATCH --mail-type=ALL                                         # Set the types of notifications which will be emailed.
#SBATCH --job-name="archive_job"

echo "Starting run at: $(date)"
for dir in N*; do
    echo $dir
    rsync -axvH --no-g --no-p  "${HOME}/scratch/HnsCub-probes-tar/${dir}" "${HOME}/amirhsi_rrg/"
    echo finished
done
echo "Run finsihed at: $(date)"

