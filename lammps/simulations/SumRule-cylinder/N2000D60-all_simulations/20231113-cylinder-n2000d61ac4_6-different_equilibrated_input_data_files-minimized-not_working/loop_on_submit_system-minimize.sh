#!/bin/bash
for directory in N*/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit_system_minimize.sh
   sleep 5
   cd ..
done
