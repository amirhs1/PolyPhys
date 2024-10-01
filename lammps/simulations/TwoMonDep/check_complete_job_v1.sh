#!/bin/bash
# prints the last 5 lines of the SLURM scheduler's sbatch's output file.
for dir in am*ens[1-4]/; do
    echo "$dir"
    echo ""
    tail -n 5 "${dir}"output.txt
    echo ""
done