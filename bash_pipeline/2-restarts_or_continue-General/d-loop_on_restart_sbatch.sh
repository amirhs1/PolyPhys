#!/bin/bash
for directory in N*_res/; do
   cd "${directory}" || exit
   echo "${directory}"
   sbatch submit.sh
   sleep 5
   cd ..
done
