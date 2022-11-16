#!/bin/bash
for dir in N*.ring/; do # HnsCub
#for dir in N*[1-8]/; do # SumRuleCyl
#for dir in al*_incomplete/; do # TransFociCub
#for dir in eps*_encomplete/; do # TransFociCyl
    comdir=$(echo "$dir" | cut -d _ -f 1)
    cd ${dir}_incomplete || exit
        cp ./*.all.lammpstrj.gz ../"${comdir}"/
        cp ./*.nucleoid.lammpstrj ../"${comdir}"/
        cp log.lammps ../"${comdir}"/
    cd ..
    cd ${dir}_res || exit
        cp N*all.data ../"${comdir}/"
        cp N*all.lammpstrj.gz ../"${comdir}"/
        cp N*nucleoid.restart.lammpstrj ../"${comdir}"/
        cp log.lammps restart.log.lammps
        mv restart.log.lammps ../"${comdir}"/
    cd ..
done