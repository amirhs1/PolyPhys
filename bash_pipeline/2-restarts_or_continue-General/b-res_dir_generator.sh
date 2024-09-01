#!/bin/bash
# create restart directories for all the directories in a parent directory.

for dir in N*ens[1-8]/; do # SumRuleCyl
#for dir in N*ring/; do # HnsCub HnsCyl
#for dir in al*ring/; do # TransFociCub, SumRuleCubHeteroRing
#for dir in al*linear/; do # SumRuleCubHeteroLinear
#for dir in eps*/; do # TransFociCyl
    file=$(echo "$dir" | cut -d / -f 1)
    echo "$file"
    resdir=${file}_res
    mkdir "$resdir"
    indir=${file}_incomplete
    mv "$dir" "$indir"
    mkdir "$file"
done