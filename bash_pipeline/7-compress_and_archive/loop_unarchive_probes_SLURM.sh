#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=00-09:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="unarchive_job"

> archive_report.txt  # Clear previous report contents
echo "Starting run at: $(date)"
for tar_file in N*.tar; do
    # Determine the output directory name by replacing ".tar" with "-zip"
    output_dir=${tar_file%.tar}-zip
    
    # Extract only the specific files from the tar file
    tar -xf "${tar_file}" --wildcards --no-anchored 'N*.lammpstrj.gz' 'N*.data.gz'

    echo "Finished extracting specific files from ${tar_file}" >> archive_report.txt
    
    # Find and decompress specific .gz files within the output directory
    find "${output_dir}" -type f \( -name "N*.lammpstrj.gz" -o -name "N*.data.gz" \) -exec gunzip {} \;
    echo "Decompressed gz files in ${output_dir}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"
