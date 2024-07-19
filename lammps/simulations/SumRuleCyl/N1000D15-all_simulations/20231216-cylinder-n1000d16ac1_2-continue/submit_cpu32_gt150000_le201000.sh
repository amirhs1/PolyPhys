#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=07-00:00
#SBATCH --account=rrg-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL

# Load the module:
module --force purge
module load StdEnv/2020 intel/2020.1.217 openmpi/4.0.3 lammps-omp/20220623

echo "Starting run at: `date`"

lmp_exec=lmp
lmp_input="input.lmp"
lmp_output="output.txt"
srun ${lmp_exec} < ${lmp_input} > ${lmp_output}

echo "Program finished with exit code $? at: `date`"