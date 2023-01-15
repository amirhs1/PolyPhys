#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-sum_rule_bug-trj_segments "$dir"
    cp -R polyphys "${dir}"
done
