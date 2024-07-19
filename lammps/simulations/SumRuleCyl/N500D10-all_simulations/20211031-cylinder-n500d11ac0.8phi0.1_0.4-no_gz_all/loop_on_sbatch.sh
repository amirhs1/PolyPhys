#!/bin/bash
for directory in */; do
   cd ${directory}
   sbatch submit.sh
   sleep 5
   cd ..
done
