#!/bin/bash
# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
pname=$(pwd | rev | cut -d / -f 1 | rev) # pattern of the parent name: D#al#nl#ns#ac#-SIMTYPE
name=$( echo "$pname" | cut -d - -f 1)
logdir=${name}-logs # name of trjactories directory
mkdir "$logdir"
for dir in al*[1-8].ring/; do
    echo "$dir"
    fname=$(echo "$dir" | cut -d / -f 1)
    cd "$dir" || exit
    cp log.lammps  "$fname".log 
    mv "$fname".log ../"$logdir"/
    cd ..
done
#move Lammps running files to run_files directory
rundir=run_files-${name}
mkdir "$rundir"
mv ./*.data ./*.lmp ./*.sh "$rundir"
cp ../PolyPhys/bash_scripts/4-probe/cubic/probe_runfiles-cubic.sh .