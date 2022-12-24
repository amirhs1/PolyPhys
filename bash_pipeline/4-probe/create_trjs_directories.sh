#!/bin/bash
# organize topology and data files base on their *whole* simulation names into directories.
for file in N*.bug.data; do
	# trunk-ignore(shellcheck/SC2086)
	dir=$(echo $file | sed -n -e 's/\(^.*\)\(\(.bug.data\).*\)/\1/p')
	echo "${dir}"
	# trunk-ignore(shellcheck/SC2086)
	trj=$(echo ${dir}*.bug.lammpstrj)
	# trunk-ignore(shellcheck/SC2086)
	segments=$(echo ${dir}*.all.lammpstrj)
	all=${dir}.all.data

	mkdir "${dir}"

	mv "${file}" "${dir}"
	mv "${trj}" "${dir}"
	mv "${segments}" "${dir}"
	mv "${all}" "${dir}"
done
