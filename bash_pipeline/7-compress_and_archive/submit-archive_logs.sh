#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=00-02:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="archive_job"

> archive_report.txt  # Clear previous report contents
echo "Starting run at: $(date)"
for dir in N*ens[1-8]-logs-2ndRound/; do # SumRuleCyl-2ndRound
#for dir in N*ens[1-8]-logs/; do # SumRuleCyl
#for dir in N*ring-logs/; do # HnsCub HnsCyl
#for dir in al*ring-logs/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear-logs/; do # TransFociCub, SumRuleCubHeteroLinear
#for dir in eps*-logs/; do # TransFociCyl
    tar -zcvf "${dir::-1}".tar.gz "${dir}" 
    echo "Finished archiving ${dir}" >> archive_report.txt
done
echo "Program finished with exit code $? at: $(date)"
