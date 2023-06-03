#!/bin/bash
for directory in eps*/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit_minimize_initial_config.sh
   sleep 5
   cd ..
done
