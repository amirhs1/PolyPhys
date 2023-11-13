#!/bin/bash
# prints the last 5 lines of the SLURM scheduler's sbatch's output file.
#for dir in al*[1-8].ring/; do # TransFociCub
#for dir in N*.ring/; do #HnsCub
for dir in N*[1-8]/; do #SumRuleCul
#for dir in eps*[1-8].ring/; do #TransFociCyl
    echo "$dir"
    # Check if log.lammps file exists in the folder
    if [[ -f "${dir}log.lammps" ]]; then
        # Read the last line of the file
        last_line=$(tail -n 1 "${dir}log.lammps")

        # Check if the last line starts with "Total wall time"
        if [[ $last_line != "Total wall time"* ]]; then
            # Copy the folder name to a desired location or file
            echo "$dir" >> folders_to_check.txt
        fi
    else
        echo "log.lammps not found in $dir"
    fi
done