#!/bin/bash
# prints the last 5 lines of the SLURM scheduler's sbatch's slurm file.
#for dir in al*/; do # TrnasFociCub
for dir in eps*/; do #TransFociCyl
    echo "$dir"
    echo
    tail -n 5 "${dir}"slurm*
    echo
done