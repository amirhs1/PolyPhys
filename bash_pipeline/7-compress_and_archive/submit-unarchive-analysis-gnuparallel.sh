#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=3G
#SBATCH --time=0-00:50
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL

# record environment to unclutter gnu parallel run
parallel --record-env

# Create a function to execute your job
exe(){
dir=${1}
(tar -zxf ${dir})
}

echo "Starting run at: $(date)"

# Run GNU Parallel
#export the function
export -f exe

# run the loop in parallel
# SumRuleCyl:
#parallel --will-cite --ungroup  --env _ exe ::: N*ens[1-8].tar.gz 
# HnsCub, HnsCyl:
#parallel --will-cite --ungroup  --env _ exe ::: N*ring.tar.gz
# TransFociCub, TransFociCyl, SumRuleCubHeteroRing, SumRuleCubHeteroLinear:
parallel --will-cite --ungroup  --env _ exe ::: ns*.tar.gz

echo "Program glost_launch finished with exit code $? at: $(date)"
