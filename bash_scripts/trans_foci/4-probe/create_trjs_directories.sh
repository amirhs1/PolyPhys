#!/bin/bash
# organize topology and data files base on their *whole* simulation names into directories.
for file in eps*.bug.data; do
	dir=$(echo "$file"| rev | cut -d . -f 3- | rev)
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
