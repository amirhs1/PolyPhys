#!/bin/bash

# script to generate folders with necessary files for a job array on the slurm schedular

for file in N*.bug.data; do
	dir="$(echo $file | sed -n -e 's/\(^.*\)\(\(.bug.data\).*\)/\1/p')"
	echo $dir
	trj="$(echo $dir.bug.lammpstrj)"
	segments="$(echo ${dir}*.all.lammpstrj)"
	log="$(echo $dir.log)"
	all="$(echo $dir.all.data)"
	
	mkdir $dir
	
	cp extraction_all_in_one.py $dir
	cp __init__.py $dir
	cp PipeLine.py $dir
	
	mv $file $dir
	mv $trj $dir
	mv $segments $dir
	mv $log $dir
	mv $all $dir
	
	cd $dir
	mkdir PipeLine
	mv __init__.py PipeLine
	mv PipeLine.py PipeLine
	
	cd ..
done
