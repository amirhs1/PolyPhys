#!/bin/bash
for dir in al*linear/; do
   cd "${dir}" || exit
   echo "${dir}"
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out min_init_config-"${slurm}".out
   cd ..
done