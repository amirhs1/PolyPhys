#!/bin/bash
for directory in N*/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit.sh
   sleep 5
   cd ..
done