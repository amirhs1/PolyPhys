#!/bin/bash

for dir in N*_res/; do # SumRuleCyl
#for dir in al*_res/; do # TransFociCub
#for dir in eps*_res/; do # TransFociCyl
    comdir=$(echo "$dir" | cut -d _ -f 1)
    echo "$comdir"
    cd "$dir" || exit
    cp all*.data ../"${comdir}/"
    cp all.*.lammpstrj.gz ../"${comdir}"/
    cp bug.*.lammpstrj ../"${comdir}"/
    cp log.lammps restart.log.lammps
    mv restart.log.lammps ../"${comdir}"/
    cd ..
done