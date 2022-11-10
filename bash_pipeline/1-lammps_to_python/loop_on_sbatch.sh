#!/bin/bash
for directory in al*/; do # TransFociCub
#for directory in eps*/; do # TransFociCyl
   cd "${directory}" || exit
   sbatch submit.sh
   sleep 5
   cd ..
done
