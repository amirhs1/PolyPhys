#!/bin/bash
for dir in eps*_incomplete/;do
    comdir=$(echo "$dir" | cut -d _ -f 1)
    cd "${dir}" || exit
    cp all.*.lammpstrj.gz ../"${comdir}"/
    cp bug.*.lammpstrj ../"${comdir}"/
    cp log.lammps ../"${comdir}"/
    cd ..
done