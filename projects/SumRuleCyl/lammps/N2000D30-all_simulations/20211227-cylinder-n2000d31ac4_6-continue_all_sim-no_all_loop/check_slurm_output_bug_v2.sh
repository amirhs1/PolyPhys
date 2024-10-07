#!/bin/bash

# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
for dir in N*/; do
    echo $dir
    cd $dir
    echo
    tail -n 10 output.txt
    echo
    tail -n 10 log.lammps
    echo
    tail -n 20 slurm*
    echo
    cd ..
done