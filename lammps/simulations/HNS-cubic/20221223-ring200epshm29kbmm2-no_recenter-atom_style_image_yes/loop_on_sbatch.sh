#!/bin/bash
for directory in N*/; do
   cd "${directory}" || exit
   eco "${directory}"
   mv log.lammps minimize_dna.lmp  
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out minimize_dna-"${slurm}".out
   sbatch submit.sh
   sleep 5
   cd ..
done
