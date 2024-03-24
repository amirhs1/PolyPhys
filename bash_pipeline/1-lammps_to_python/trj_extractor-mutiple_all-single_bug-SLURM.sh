#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=00-10:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="extract_job"

echo "Starting run at: $(date)"
report_name="$(basename "$(pwd)")-extract_report.txt"
echo "Starting run at: $(date)"
touch "$report_name"

echo "Starting run at: $(date)" >> $report_name 
# This script copies the dump bug/momomer file from a run folder to the simulation folder and rename it to the run folder name
echo "Start to extract ..." >> $report_name
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
    echo "not a 'all_simulations' or 'cont_simulations' directory" >> $report_name
    exit 1
fi
echo "Trajectory directory: $trjdir" >> $report_name
mkdir "$trjdir"
mkdir "$logdir"
#for dir in al*[1-8].ring/; do # TransFociCub
#for dir in eps*[1-8].ring/; do # TransFociCyl
for dir in N*[1-8]/; do # SumRuleCyl
#for dir in N*[1-8].ring/; do
    echo "$dir"
    fname=$(echo "$dir" | cut -d / -f 1)
    mkdir "$trjdir"/"$fname"
    cd "$dir" || exit
    cp "$fname".bug.lammpstrj ../"$trjdir"/"$fname"/
    cp "$fname".all.data ../"$trjdir"/"$fname"/
    for gzfile in N*.all.lammpstrj.gz;do 
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
mv ./*.data ./*.lmp ./*.sh "$rundir"
echo "Finished!" >> $report_name
echo "Program finished with exit code $? at: $(date)" >> $report_name
echo "Program finished with exit code $? at: $(date)"
