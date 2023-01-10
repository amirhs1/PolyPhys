#!/bin/bash
# prints the last 5 lines of the SLURM scheduler's sbatch's slurm file.
#for dir in al*[1-8].ring/; do # TransFociCub
#for dir in N*.ring/; do #HnsCub
#for dir in N*[1-8]/; do #SumRuleCul
for dir in eps*[1-8].ring/; do #TransFociCyl
    echo "$dir"
    echo
    tail -n 5 "${dir}"slurm*
    echo
done