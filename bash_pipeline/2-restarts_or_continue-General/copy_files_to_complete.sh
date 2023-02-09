#!/bin/bash
echo "Start copying ..."
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
        cp ./log.lammps "${comdir}".incomplete.log
        mv ./"${comdir}".incomplete.log ../"${comdir}"/
    cd ..
    cd "${comdir}"_res || exit
        cp ./*all.data ../"${comdir}/"
        cp ./*all.restart.lammpstrj.gz ../"${comdir}"/
        #cp ./*nucleoid.restart.lammpstrj ../"${comdir}"/ # HnsCub
        cp ./*bug.restart.lammpstrj ../"${comdir}"/ # others
        cp ./log.lammps "${comdir}".restart.log
        mv ./"${comdir}".restart.log ../"${comdir}"/
    cd ..
done
echo "Finished!"