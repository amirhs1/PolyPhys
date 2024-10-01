#!/bin/bash

# Test run to estimate the resources needed for the final submissions
# Date:2019 March 13

#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=00-00:40
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL

# Load the module:

module --force purge
module load StdEnv/2023  intel/2023.2.1  openmpi/4.1.5 lammps-omp/20230802

echo "Starting run at: $(date)"

lmp_exec=lmp
lmp_input="minimize_initial_config.lmp"
lmp_output="minimize_initial_config-output.txt"

srun ${lmp_exec} -log log.minimize_initial_config < ${lmp_input} > ${lmp_output}

echo "Program finished with exit code $? at: $(date)"