#!/bin/bash
# changes the padding of the 'all' trajectories by increasing j value of a segment by 1.

for allfile in al*all.lammpstrj; do # TransFociCub
#for allfile in eps*all.lammpstrj; do# TransFociCyl
    name=$( echo "$allfile" | cut -d . -f -6) # whole name without j
    echo "$name"
    j=$( echo "$allfile" | cut -d . -f 7) # j value
    segment=${j:1}
    newj=$(( 10#$segment + 1))
    jpad=$(printf "%02d" $newj)
    newname="${name}.j${jpad}.ring.all.lammpstrj.finished"
    mv "$allfile" "$newname"
done