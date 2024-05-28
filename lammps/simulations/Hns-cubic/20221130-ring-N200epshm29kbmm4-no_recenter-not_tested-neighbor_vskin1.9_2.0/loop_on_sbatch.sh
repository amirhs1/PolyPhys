#!/bin/bash
for dir in N*/; do
   cd "${dir}" || exit
   echo "${dir}"
   sbatch submit.sh
   sleep 5
   cd ..
done
