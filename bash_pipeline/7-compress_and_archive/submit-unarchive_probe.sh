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
for tar_file in N*ens[1-8].tar.gz; do # SumRuleCyl
#for tar_file in N*ring.tar.gz; do # HnsCub HnsCyl
#for tar_file in al*ring.tar.gz; do # TransFociCub, SumRuleCubHeteroRing
#for tar_file in al*linear.tar.gz; do # TransFociCub, SumRuleCubHeteroLinear
#for tar_file in eps*.tar.gz; do # TransFociCyl
    # Determine the output directory name by replacing ".tar" with "-zip"
    output_dir=${tar_file%.tar.gz}-zip
    # Extract only the specific files from the tar file
    tar -zxf "${tar_file}"
    echo "Finished extracting specific files from ${tar_file}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"