#!/bin/bash
for dir in N*.ring_incomplete/; do # HnsCub
#for dir in N*[1-8]_incomplete/; do # SumRuleCyl
#for dir in al*_incomplete/; do # TransFociCub
#for dir in eps*_encomplete/; do # TransFociCyl
    comdir=$(echo "$dir" | cut -d _ -f 1)
    mkdir "${comdir}"
    cd "${comdir}"_incomplete || exit
        cp ./*.all.lammpstrj.gz ../"${comdir}"/
        #cp ./*.nucleoid.lammpstrj ../"${comdir}"/ # HnsCub 
        cp ./*.bug.lammpstrj ../"${comdir}"/ # Others
        cp log.lammps ../"${comdir}"/
    cd ..
    cd "${comdir}"_res || exit
        cp N*all.data ../"${comdir}/"
        cp N*all.lammpstrj.gz ../"${comdir}"/
        #cp N*nucleoid.restart.lammpstrj ../"${comdir}"/ # HnsCub
        cp N*bug.restart.lammpstrj ../"${comdir}"/ # others
        cp log.lammps restart.log.lammps
        mv restart.log.lammps ../"${comdir}"/
    cd ..
done