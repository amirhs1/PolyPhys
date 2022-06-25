#!/bin/bash
# changes the padding of the 'all' trajectories by increasing j value of a segment by 1.
for allfile in N*all.lammpstrj; do
    name=$( echo "$allfile" | cut -d . -f -6)
    echo "$name"
    j=$( echo "$allfile" | cut -d . -f 7)
    segment=${j:1}
    newj=$(( 10#$segment + 1))
    jpad=$(printf "%02d" $newj)
    newname="${name}.j${jpad}.all.lammpstrj.finished"
    mv "$allfile" "$newname"
done