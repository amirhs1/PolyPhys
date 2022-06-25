#!/bin/bash
# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
pname=$(pwd | rev | cut -d / -f 1 | rev) #parent name
name=$( echo "$pname" | cut -d - -f 1)
simtype=$( echo "$pname" | cut -d - -f 2 | cut -d _ -f 1) # whether the simulations are "all" or "cont"; "all" means that the first time simulations are run; "cont" means this is the 2nd time they are run since the system has not yet reached equilibrium.
if [ "$simtype" = "all" ]; then
    trjdir=${name}-trjs
elif [ "$simtype" = "cont" ]; then
    trjdir=${name}-${simtype}_trjs
else
    echo "not a 'all_simulations' or 'cont_simulations' directory"
    exit 1
fi
echo "Trajectory directory: $trjdir"
mkdir "${trjdir}"
for dir in N*esp*[1-8]/; do
    echo "$dir"
    fname=$(echo "$dir" | cut -d / -f 1)
    mkdir "$trjdir"/"$fname"
    cd "$dir" || exit
    cp bug.i*.lammpstrj "${fname}".j02.bug.lammpstrj
    mv "$fname".j02.bug.lammpstrj ../"$trjdir"/"$fname"/
    cp log.lammps  "$fname".cont.log 
    mv "$fname".cont.log ../"${trjdir}"/"$fname"/
    gzip -dk  all.*.lammpstrj.gz
    cp all.*.lammpstrj "${fname}".j02.all.lammpstrj 
    mv "${fname}".j02.all.lammpstrj ../"${trjdir}"/"$fname"/
    cd ..
done

#move Lammps runnning files to run_files directory
rundir=run_files-${name}
mkdir "${rundir}"
mv ./*.data ./*.lmp ./*.sh "${rundir}"