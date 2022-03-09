#!/bin/bash
# copies the python script asnd package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-all_trj_segments*.py probe-all_trj_segments.py
    mv probe-all_trj_segments.py "$dir"
    cp -R polyphys "${dir}"
done
