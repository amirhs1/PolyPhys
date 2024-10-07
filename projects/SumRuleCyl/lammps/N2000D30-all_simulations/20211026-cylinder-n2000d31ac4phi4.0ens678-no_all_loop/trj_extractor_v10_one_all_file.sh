#!/bin/bash
# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
trjdir=${name}-extraction_bug
echo "Trajectory directory: $name"
mkdir ${trjdir}
for dir in N*[1-8]/; do
    echo $dir
    cd $dir
    fname=$(echo ${dir} | cut -d / -f 1)
    cp bug.i*.lammpstrj $fname.bug.lammpstrj
    mv $fname.bug.lammpstrj ../${trjdir}/
    cp all.*.data $fname.all.data
    mv $fname.all.data ../${trjdir}/
    cp log.lammps  $fname.log 
    mv $fname.log ../${trjdir}/
    gzip -dk  all.*.lammpstrj.gz
    cp all.*.lammpstrj ${fname}.all.lammpstrj 
    mv ${fname}.all.lammpstrj ../${trjdir}/
    cd ..
done