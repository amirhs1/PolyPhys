#!/bin/bash
# changes the padding of the 'all' restart trajectories by increasing j value of a segment by 1.
for dir in al*ring; do
    cd "$dir" || exit
    ls -lh al*.ring.all.lammpstrj.gz
    com=$(ls -lh al*.ring.all.lammpstrj.gz | wc -l)
    echo com "$com"
    rm "${dir::-4}"j"${com}".ring.all.lammpstrj.gz
    for allfile in al*restart.lammpstrj.gz; do # TransFociCub
    #for allfile in eps*all.lammpstrj; do# TransFociCyl
        name=$( echo "$allfile" | cut -d . -f -2) # whole name without j
        j=$( echo "$allfile" | cut -d . -f 3) # j value
        segment=${j:1}
        newj=$(( 10#$segment))
        jj=$((${com}+${newj}-1))
        jpad=$(printf "%02d" $jj)
        newname="${name}.j${jpad}.ring.all.lammpstrj.gz.finished"
        mv "$allfile" "$newname"
    done
    cd ..
done