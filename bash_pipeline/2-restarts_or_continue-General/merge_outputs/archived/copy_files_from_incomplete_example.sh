#!/bin/bash
for dir in N*_incomplete/; do # SumRuleCyl
#for dir in al*_incomplete/; do # TransFociCub
#for dir in eps*_encomplete/; do # TransFociCyl
    comdir=$(echo "$dir" | cut -d _ -f 1)
    cd "${dir}" || exit
    cp ./*.all.lammpstrj.gz ../"${comdir}"/
    cp ./*.bug.lammpstrj ../"${comdir}"/
    cp log.lammps ../"${comdir}"/
    cd ..
done