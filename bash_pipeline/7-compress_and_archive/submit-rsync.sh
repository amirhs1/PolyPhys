#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=00-05:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="archive_job"

echo "Starting run at: $(date)"
for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear/; do # TransFociCub, SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
    echo $dir
    rsync -axvH --no-g --no-p  "${HOME}/scratch/${dir}" "${HOME}/amirhsi_rrg/"
    echo finished
done
echo "Run finsihed at: $(date)"

