#!/bin/bash
# organize topology and data files base on their *whole* simulation names into directories.
for file in N*.bug.data; do
	dir=$(echo "$file" | sed -n -e 's/\(^.*\)\(\(.bug.data\).*\)/\1/p')
	echo "${dir}"
	trj=$(echo "${dir}"*.bug.lammpstrj)
	segments=$(echo "${dir}"*.all.lammpstrj)
	all=${dir}.all.data

	mkdir "${dir}"

	mv "${file}" "${dir}"
	mv "${trj}" "${dir}"
	mv "${segments}" "${dir}"
	mv "${all}" "${dir}"
done
