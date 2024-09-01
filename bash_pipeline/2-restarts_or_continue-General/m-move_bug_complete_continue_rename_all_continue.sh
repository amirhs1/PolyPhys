#!/bin/bash
# CHECK continue name; this script is simple
for dir in N*ens[1-8]_cont/; do # SumRuleCyl
#for dir in N*ring_cont/; do # HnsCub HnsCyl
#for dir in al*ring_cont/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear_cont/; do # SumRuleCubHeteroLinear
#for dir in eps*_cont/; do # TransFociCyl
    cd ${dir::-5}
    rm ${dir::-5}.bug.lammpstrj
    mv test/${dir::-5}*lammpstrj .
    mv ${dir::-5}.j21.all.continue.lammpstrj.gz ${dir::-5}.j21.all.lammpstrj.gz
    rm -r test
    cd ..
done
