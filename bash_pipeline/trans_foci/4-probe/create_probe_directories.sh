#!/bin/bash
# create a 'probe' driectory for a given space and moves all the probe files to the directory.
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo "$currentname" | cut -d - -f 1)
probe=${name}-probe
mkdir "${probe}"
cp slurm*.out "${probe}"-slurm_report.out 
mv "${probe}"-slurm_report.out ./"${probe}"/
for dir in eps*/; do
    mkdir "${probe}/$dir"
    echo "$dir"
    cd "$dir" || exit
    mv ./eps*.csv "../${probe}/$dir"
    mv ./eps*.npy "../${probe}/$dir"
    mv ./*.txt "../${probe}/$dir"
    cd ..
done