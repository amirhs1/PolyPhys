#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-hns_cub-all_in_one-nucleoid_wholes.py ./"$dir"
    cp -R polyphys ./"${dir}"
done
