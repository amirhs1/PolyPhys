#!/bin/bash
# script to generate folders with necessary files for a job array on the slurm schedular
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
analyzedir=${name}-analyze_bug
echo "Trajectory directory: $name"
mkdir ${analyzedir}
mv N*.csv $analyzedir/
mv $analyzedir ..