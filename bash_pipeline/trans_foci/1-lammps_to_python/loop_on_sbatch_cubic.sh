#!/bin/bash
for directory in al*/; do
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
