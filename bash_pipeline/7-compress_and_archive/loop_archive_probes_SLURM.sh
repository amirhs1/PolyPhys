#!/bin/bash

#SBATCH --ntasks=1                                                     
#SBATCH --mem-per-cpu=3000M                   
#SBATCH --time=01-00:00
#SBATCH --account=def-group
#SBATCH --mail-user=54184877+amirhs1@users.noreply.github.com 
#SBATCH --mail-type=ALL
#SBATCH --job-name="archive_job"

> archive_report.txt  # Clear previous report contents
echo "Starting run at: $(date)"
for dir in N*/; do
    tar -zcvf "${dir::-1}".tar.gz "${dir}" 
    echo "Finished archiving ${dir}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"