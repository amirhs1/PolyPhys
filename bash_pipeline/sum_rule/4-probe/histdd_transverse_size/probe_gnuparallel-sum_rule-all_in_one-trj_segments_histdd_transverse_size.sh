#!/bin/bash
# copies the python script asnd package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-sum_rule-all_in_one-trj_segments*.py probe-sum_rule-all_in_one-trj_segments_histdd_transverse_size.py
    mv probe-sum_rule-all_in_one-trj_segments_histdd_transverse_size.py "$dir"
    cp -R polyphys "${dir}"
done
