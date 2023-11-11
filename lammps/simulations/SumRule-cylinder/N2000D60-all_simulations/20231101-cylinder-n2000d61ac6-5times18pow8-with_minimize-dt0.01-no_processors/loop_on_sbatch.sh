#!/bin/bash
for directory in N*/; do
   cd "${directory}" || exit
   echo "${directory}"
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out minimize_chain-"${slurm}".out
   sbatch submit.sh
   sleep 5
   cd ..
done
