#!/bin/bash

# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
for dir in N*/; do
    cd $dir
    lastline=$(tail -n 1 slurm*)
    if [[ $lastline == *"exit code 137"* ]]; then
        echo $dir
        rm log.lammps
        rm output.txt
        rm slurm*
        cp ../chain.1000.data .
    fi
    cd ..
done