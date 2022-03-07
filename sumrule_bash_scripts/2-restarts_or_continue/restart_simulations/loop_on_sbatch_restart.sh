#!/bin/bash
for directory in N*_res/; do
   cd ${directory}
   sbatch submit.sh
   sleep 5
   cd ..
done