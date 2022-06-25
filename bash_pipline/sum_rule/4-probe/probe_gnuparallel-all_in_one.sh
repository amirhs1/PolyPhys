#!/bin/bash
# copies the python script asnd package which are then executed by gnuparallal.
for dir in N*/; do
    echo "$dir"
    cp probe-all_in_one_*.py probe-all_in_one.py
    mv probe-all_in_one.py "$dir"
    cp -R polyphys "${dir}"
done
