#!/bin/bash
# create a 'probe' driectory for a given space and moves all the probe files to the directory.
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo "$currentname" | cut -d - -f 1)
probe=${name}-probe
mkdir "${probe}"
job=$(echo slurm*.out | cut -d . -f 1)
cp slurm*.out "${probe}-${job}_report.out" 
mv "${probe}"-"${job}"_report.out ./"${probe}"/
for dir in al*linear/; do
    mkdir "${probe}/$dir"
    echo "$dir"
    cd "$dir" || exit
    mv ./al*.csv "../${probe}/$dir"
    mv ./al*.npy "../${probe}/$dir"
    mv ./al*.txt "../${probe}/$dir"
    cd ..
done
echo "Finished!"