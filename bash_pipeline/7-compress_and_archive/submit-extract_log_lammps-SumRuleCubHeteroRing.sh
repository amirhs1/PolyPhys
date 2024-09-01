#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem-per-cpu=3G
#SBATCH --time=0-05:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL

# Record environment to unclutter GNU Parallel run
parallel --record-env

# Create a function to execute your job
exe(){
    dir=${1}
    cd "$dir" || { echo "Failed to enter directory $dir"; exit 1; }
    find . -name "al*.tar" -exec tar --wildcards -xvf {} "al*-zip/log.lammps.gz" \;
    cd ..
}

echo "Starting run at: $(date)"

# Export the function so that GNU Parallel can use it
export -f exe

# Run the loop in parallel across directories matching ns*archived
parallel --will-cite --ungroup --env _ exe ::: ns*-all_simulations-archived

echo "Program finished with exit code $? at: $(date)"

