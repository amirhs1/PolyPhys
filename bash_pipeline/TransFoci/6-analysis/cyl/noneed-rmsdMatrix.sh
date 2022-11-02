#!/bin/bash

# script to generate folders with necessary files for a job array on the slurm schedular
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
rmsddir=${name}-rmsd
echo "Trajectory directory: $name"
mkdir ${rmsddir}

for dir in N*[1-8]/; do
        echo $dir

        cd $dir 
            cp *.npy ../${rmsddir}
        cd ..
        
done

mv $rmsddir ..