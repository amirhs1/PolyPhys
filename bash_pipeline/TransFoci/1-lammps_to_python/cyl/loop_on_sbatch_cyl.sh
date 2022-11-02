#!/bin/bash
for directory in eps*/; do
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
