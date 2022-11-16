#!/bin/bash

for dir in N*[1-8]/; do # SumRuleCyl
#for dir in al*[1-8]/; do # TransFociCub
#for dir in eps*[1-8]/; do # TransFociCyl
    echo "$dir"
    results=$(echo "$dir" | cut -d . -f 1 |  grep -o -E '[0-9]+')
    nmon=$( echo "$results" | cut -d " " -f 1)
    results=$(echo "$dir" | cut -d . -f 5 |  grep -o -E '[0-9]+')
    ncrowd=$( echo "$results" | cut -d " " -f 2)
    natoms=$(( "$ncrowd" + "$nmon" ))
    echo "number of atoms: " "$natoms"
    cd "$dir" || exit
    lines=0
    j=0
    for file in N*.j*.all.lammpstrj; do
        line=$(wc -l "$file" | cut -d " " -f 1)
        lines=$((10#$lines + 10#$line))
        j=$(($j + 1))
    done
    echo "Totsl number of all trjs: " "$j"
    echo "Total number of lines: " "$lines"
    lpf=$(( "$natoms" + 9)) # lines per frame
    sampled=$(( "$lines" / "$lpf" ))
    echo "Number of sampled frames: " "$sampled"
    echo "Number of all frames: " "$(( ( ( "$sampled" - "$j" ) * 5000 + 71000000) ))"
    cd ..
done
