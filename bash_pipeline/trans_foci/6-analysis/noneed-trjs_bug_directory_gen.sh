#!/bin/bash

# script to generate folders with necessary files for a job array on the slurm schedular
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
trjsdir=${name}-trjs_bug
logsdir=${name}-logs
echo "Trajectory directory: $name"
mkdir ${trjsdir}
mkdir ${logsdir}

for dir in N*[1-8]/; do
	echo $dir
	cd $dir 

	rm extraction*.py
	rm -r PipeLine

	mv *.lammpstrj ..
	mv *.all.data ..
	mv *.bug.data ..
	mv *.log ..
	cd ..
	
	mv *.lammpstrj $trjsdir/
	mv *.all.data $trjsdir/
	mv *.bug.data $trjsdir/
	mv *.log $logsdir/
done

mv $trjsdir ..
mv $logsdir ..