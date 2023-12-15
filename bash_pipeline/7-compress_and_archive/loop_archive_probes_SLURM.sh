#!/bin/bash

#SBATCH --ntasks=1                                                     
#SBATCH --mem-per-cpu=3000M                   
#SBATCH --time=01-00:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com 
#SBATCH --mail-type=ALL
#SBATCH --job-name="archive_job"

> archive_report.txt  # Clear previous report contents
echo "Starting run at: $(date)"
for dir in N*/; do
    tar -zcvf "${dir::-1}".tar.gz "${dir}" 
    echo "Finished archiving ${dir}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"