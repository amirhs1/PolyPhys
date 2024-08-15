#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in al*ring/; do
    echo "$dir"
    cp probe-trans_foci_cub-all_in_one-bug_wholes.py  ./"$dir"
    cp -R polyphys ./"${dir}"
done
