#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in al*/; do
    echo "$dir"
    cp probe-trans_foci_bug_cub-bug_wholes.py ./"$dir"
    cp -R ./polyphys ./"${dir}"
done
