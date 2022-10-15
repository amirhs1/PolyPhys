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
    for gzfile in all.*.lammpstrj.gz;do 
            gzip -dk $gzfile
            allfile=$( echo $gzfile | cut -d . -f -4)
            j=$(echo $allfile | cut -d . -f 3)
            mv $allfile ${fname}.${j}.all.lammpstrj 
            mv ${fname}.${j}.all.lammpstrj ../${trjdir}/
    done
    cd ..
done
