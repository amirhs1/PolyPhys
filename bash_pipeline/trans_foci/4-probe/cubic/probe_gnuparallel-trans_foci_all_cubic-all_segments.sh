#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in al*/; do
    echo "$dir"
    cp probe-trans_foci_all_cubic-all_segments-*.py probe-trans_foci_all_cubic-all_segments.py
    mv probe-trans_foci_all_cubic-all_segments.py "$dir"
    cp -R polyphys "${dir}"
done