#!/bin/bash
# generate folders with necessary files for a job array on the slurm schedular
for file in N*.bug.data; do
	dir=$(echo "$file" | sed -n -e 's/\(^.*\)\(\(.bug.data\).*\)/\1/p')
	echo "${dir}"
	trj=$(echo "${dir}"*.bug.lammpstrj)
	segments=$(echo "${dir}"*.all.lammpstrj)
	all="${dir}".all.data

	mkdir "${dir}"
	
	cp extraction_all_in_one_*.py extraction_all_in_one.py
	mv extraction_all_in_one.py "${dir}"
	
	cp -R polyphys "${dir}"

	mv "${file}" "${dir}"
	mv "${trj}" "${dir}"
	mv "${segments}" "${dir}"
	mv "${all}" "${dir}"
done
