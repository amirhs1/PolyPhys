#!/bin/bash
for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear/; do # SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
    cd $dir
    rm ${dir}*bug*lammpstrj
    mv test/${dir}*lammpstrj .
    rm -r test
    cd ..
done