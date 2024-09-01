#!/bin/bash
# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
echo "Start to extract ..."
pname=$(pwd | rev | cut -d / -f 1 | rev) # pattern of the parent name: D#al#nl#ns#ac#-SIMTYPE
name=$( echo "$pname" | cut -d - -f 1)
simtype=$( echo "$pname" | cut -d - -f 2 | cut -d _ -f 1) # whether the simulations are "all" or "cont"; "all" means that the first time simulations are run; "cont" means this is the 2nd time they are run since the system has not yet reached equilibrium.
if [ "$simtype" = "all" ]; then
    trjdir=${name}-trjs # name of trjactories directory
    logdir=${name}-logs # name of trjactories directory
elif [ "$simtype" = "cont" ]; then
    trjdir=${name}-trjs_${simtype}
    logdir=${name}-logs_${simtype} # name of trjactories directory
else
    echo "not a 'all_simulations' or 'cont_simulations' directory"
    exit 1
fi
echo "Trajectory directory: $trjdir"
mkdir "$trjdir"
mkdir "$logdir"
#for dir in al*[1-8].ring/; do # TransFociCub
#for dir in eps*[1-8].ring/; do # TransFociCyl
for dir in N*[1-8]/; do # SumRuleCyl
#for dir in N*[1-8].ring/; do # HnsCyl, HnsCub
    echo "$dir"
    fname=$(echo "$dir" | cut -d / -f 1)
    mkdir "$trjdir"/"$fname"
    cd "$dir" || exit
    cp "$fname".nucleoid.lammpstrj ../"$trjdir"/"$fname"/
    cp "$fname".all.data ../"$trjdir"/"$fname"/
    for gzfile in N*.all.lammpstrj.gz;do 
    #for gzfile in al*.all.lammpstrj.gz;do 
    #for gzfile in eps*.all.lammpstrj.gz;do 
            gzip -dk "$gzfile"
            allfile=${gzfile[*]:0: -3}
            mv "$allfile" ../"$trjdir"/"$fname"/
    done
    cp log.lammps  "$fname".log 
    mv "$fname".log ../"$logdir"/
    cd ..
done
# move Lammps running files to run_files directory
rundir=run_files-${name}
mkdir "$rundir"
mv ./*.data ./*.lmp ./*.sh ./*.mol "$rundir"
echo "Finished!"