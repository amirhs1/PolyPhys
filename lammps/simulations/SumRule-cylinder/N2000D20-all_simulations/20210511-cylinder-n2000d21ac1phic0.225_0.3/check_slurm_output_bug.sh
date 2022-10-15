#!/bin/bash

# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
for dir in N*/; do
    echo $dir
    cd $dir
    tail -n 1 output.txt
    tail -n 100 bug.i*
    tail -n 20 slurm*
    cd ..
done
