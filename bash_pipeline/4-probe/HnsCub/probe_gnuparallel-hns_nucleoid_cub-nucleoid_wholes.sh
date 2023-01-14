#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-hns_nucleoid_cub-nucleoid_wholes.py ./"$dir"
    cp -R polyphys ./"${dir}"
done
