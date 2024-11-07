#!/bin/bash
echo "Start to submit"
for directory in al*linear/; do # TransFociCub
   cd "${directory}" || exit
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out minimize_initial_config-"${slurm}".out
   sbatch submit.sh
   sleep 5
   cd ..
done
echo "Finished!"
