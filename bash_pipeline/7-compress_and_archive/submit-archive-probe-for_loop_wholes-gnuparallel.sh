#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=3G
#SBATCH --time=0-09:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL

# record environment to unclutter gnu parallel run
parallel --record-env

# Create a function to execute your job
exe(){
parent=${1}
cd $parent
#for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
for dir in al*linear/; do # TransFociCub, SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
    file=$(echo "$dir" | cut -d / -f 1)
    tar -zcf ${file}.tar.gz ${file}
    done
cd ..
}

echo "Starting run at: $(date)"
# Run GNU Parallel
#export the function
export -f exe
# SumRuleCyl:
#parallel --will-cite --ungroup  --env _ exe ::: N*ens[1-8]-probe/
# HnsCub, HnsCyl:
#parallel --will-cite --ungroup  --env _ exe ::: N*ring-probe/
# TransFociCub, TransFociCyl, SumRuleCubHeteroRing, SumRuleCubHeteroLinear:
parallel --will-cite --ungroup  --env _ exe ::: ns*-probe/

echo "Program glost_launch finished with exit code $? at: $(date)"
