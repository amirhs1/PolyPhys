#!/bin/bash
for directory in al*ring/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit_init_config_minimize.sh
   sleep 5
   cd ..
done