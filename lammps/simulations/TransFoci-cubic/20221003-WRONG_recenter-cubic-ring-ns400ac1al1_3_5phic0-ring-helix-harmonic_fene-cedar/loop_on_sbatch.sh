#!/bin/bash
for directory in ns*/; do
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
