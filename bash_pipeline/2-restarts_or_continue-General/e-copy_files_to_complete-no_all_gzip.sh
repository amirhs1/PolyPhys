#!/bin/bash
echo "Start copying ..."
#for dir in al*ring_incomplete/; do # SumRuleHeteroRing
for dir in al*linear_incomplete/; do # SumRuleHeteroLinear
    comdir=$(echo "$dir" | cut -d _ -f 1)
    cd "${comdir}"_incomplete || exit
        cp ./*.all.lammpstrj ../"${comdir}"/ 
        cp ./*.bug.lammpstrj ../"${comdir}"/ # Others
        cp ./log.lammps "${comdir}".incomplete.log
        mv ./"${comdir}".incomplete.log ../"${comdir}"/
    cd ..
    cd "${comdir}"_res || exit
        cp ./*all.data ../"${comdir}/"
        cp ./*all.restart.lammpstrj ../"${comdir}"/
        cp ./*bug.restart.lammpstrj ../"${comdir}"/ # others
        cp ./log.lammps "${comdir}".restart.log
        mv ./"${comdir}".restart.log ../"${comdir}"/
    cd ..
done
echo "Finished!"
