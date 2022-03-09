#!/bin/bash
# copies the python script asnd package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-bug_trj_segments*.py probe-bug_trj_segments.py
    mv probe-bug_trj_segments.py "$dir"
    cp -R polyphys "${dir}"
done
