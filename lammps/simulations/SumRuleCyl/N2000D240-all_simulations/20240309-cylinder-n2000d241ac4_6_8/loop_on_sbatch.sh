#!/bin/bash
echo "Start to submit"
for directory in al*/; do # TransFociCub
#for directory in eps*/; do # TransFociCyl
   cd "${directory}" || exit
   slurm=$(echo slurm*.out | cut -d . -f 1)
   mv "${slurm}".out minimize_system-"${slurm}".out
   sbatch submit.sh
   sleep 5
   cd ..
done
echo "Finished!"
