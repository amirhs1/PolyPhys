#!/bin/bash
for dir in N*/; do
   cd "${dir}" || exit
   echo "${dir}"
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out minimize_dna-"${slurm}".out
   cd ..
done