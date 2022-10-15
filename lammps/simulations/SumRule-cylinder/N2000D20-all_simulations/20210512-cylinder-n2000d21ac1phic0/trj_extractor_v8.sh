#!/bin/bash

# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
for dir in N*/; do
    echo $dir
    ncrowd=$(echo $dir | cut  -d'c' -f 2 | cut -d'e' -f 1)
    cd $dir
    fname=$(echo ${dir} | cut -d / -f 1)
    cp bug.i*.lammpstrj $fname.bug.lammpstrj
    mv $fname.bug.lammpstrj ..
    cp all.*.data $fname.all.data
    mv $fname.all.data ..
    cp log.lammps  $fname.log
    mv $fname.log ..
    cd ..
done
