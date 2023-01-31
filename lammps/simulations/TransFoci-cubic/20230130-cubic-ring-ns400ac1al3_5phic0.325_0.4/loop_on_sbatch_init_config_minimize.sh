#!/bin/bash
for directory in al*/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit_init_config_minimize_dna.sh
   sleep 5
   cd ..
done
