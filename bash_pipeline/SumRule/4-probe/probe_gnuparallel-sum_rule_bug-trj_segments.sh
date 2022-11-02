#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-sum_rule_bug-trj_segments*.py probe-sum_rule_bug-trj_segments.py
    mv probe-sum_rule_bug-trj_segments.py "$dir"
    cp -R polyphys "${dir}"
done
