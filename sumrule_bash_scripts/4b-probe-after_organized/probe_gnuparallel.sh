#!/bin/bash
# copies the python script asnd package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe_all_in_one_*.py probe_all_in_one.py
    mv probe_all_in_one.py "$dir"
    cp -R polyphys "${dir}"
done
