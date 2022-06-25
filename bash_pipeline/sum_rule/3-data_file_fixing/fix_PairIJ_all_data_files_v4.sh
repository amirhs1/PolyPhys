#!/bin/bash
# fix the 'all' data files with coreecting the number of crowders.
# pattern of the currentname (parent directory): N*D*ac*-trjs
currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo "${currentname}"| cut -d - -f 1)
temp=$( echo "${name}" | grep -Eo '[a-zA-Z\-]+|[0-9\.]+')
# This part is different between different projects
ac=$(echo "${temp}" | sed '6q;d') # position of "ac" in the pattern is 6
echo "crowders size:" "${ac}"
for datafile in eps*.all.data; do
    echo "${datafile}"
    # the sed command below (s)ubstitute the line starts with
    # "PairIJ Coeffs" with the "Pair Coeffs # lj_cut"
    # Here we assumed that the Pair Coeffs are for the "lj_cut"
    # style.
    # mac osx sed:
    #sed -i "" "s/^PairIJ Coeffs.*/Pair Coeffs # lj_cut/" "${datafile}"
    # linux sed: 
    sed -i "s/^PairIJ Coeffs.*/Pair Coeffs # lj_cut/" "${datafile}"
done
echo "report:"
for datafile in eps*.all.data;do
    echo "${datafile}"
    # Check the Pair info line with new pattern is in the file
    grep -n -E "Pair Coeffs" "${datafile}"
done