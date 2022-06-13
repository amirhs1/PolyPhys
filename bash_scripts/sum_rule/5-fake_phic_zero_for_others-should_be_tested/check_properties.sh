#!/bin/bash
# check the 'stamp' files to see whether 'ncrowd' is defined correctly or not.
for dir in N*nc0*/;do
    echo "$dir"
    cd "$dir" || exit
    tail -n 2 N*-stamps.csv
    cd ..
done