#!/bin/bash

# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
for dir in N*/; do
    cd $dir
    lastline=$(tail -n 1 slurm*)
    if [[ $lastline == *"exit code 0"* ]]; then
        echo $dir
        ls .
        tail -n 50 slurm*
        echo
        tail -n 10 output.txt
    fi
    cd ..
done