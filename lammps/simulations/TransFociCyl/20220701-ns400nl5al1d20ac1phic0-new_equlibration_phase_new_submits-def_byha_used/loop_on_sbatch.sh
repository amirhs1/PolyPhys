#!/bin/bash
for directory in epss*/; do
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
