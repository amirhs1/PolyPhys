#!/bin/bash

for dir in N*_incomplete/; do # SumRuleCyl
#for dir in al*_incomplete/; do # TransFociCub
#for dir in eps*_incomplete/; do # TransFociCyl
    sim=$( echo "$dir" | cut -d _ -f 1)
    echo "$sim"
    mkdir rerun/"$sim"
    cd "$dir" || exit
    cp input.lmp ../rerun/"$sim"/
    cp submit.sh ../rerun/"$sim"/
    cp "${sim}".restart ../rerun/"$sim"/
    cd ..
done
