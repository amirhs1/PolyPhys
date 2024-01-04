#!/bin/bash
# create restart directories for all the directories in a parent directory.

for dir in N*/; do # SumRuleCyl
#for dir in al*/; do # TransFociCub
#for dir in eps*/; do # TransFociCyl
    file=$(echo "$dir" | cut -d / -f 1)
    echo "$file"
    resdir=${file}_res
    mkdir "$resdir"
    indir=${file}_incomplete
    mv "$dir" "$indir"
    mkdir "$file"
    cp "${indir}"/input.lmp "$resdir"
    cp "${indir}"/submit.sh "$resdir"
done