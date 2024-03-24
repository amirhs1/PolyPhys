#!/bin/bash
# CHECK continue name; this script is simple
for dir in N*[1-8]_cont; do
    cd ${dir::-5}
    rm ${dir::-5}.bug.lammpstrj
    mv test/${dir::-5}*lammpstrj .
    mv ${dir::-5}.j21.all.continue.lammpstrj.gz ${dir::-5}.j21.all.lammpstrj.gz
    rm -r test
    cd ..
done
