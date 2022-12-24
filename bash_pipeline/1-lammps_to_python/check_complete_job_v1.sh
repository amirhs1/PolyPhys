#!/bin/bash
# prints the last 5 lines of the SLURM scheduler's sbatch's output file.
#for dir in al*/; do # TransFociCub
for dir in eps*/; do # TransFociCyl
    echo "$dir"
    echo ""
    tail -n 5 "${dir}"output.txt
    echo ""
    tail -n 5 "${dir}"log.lammps
    echo ""
done