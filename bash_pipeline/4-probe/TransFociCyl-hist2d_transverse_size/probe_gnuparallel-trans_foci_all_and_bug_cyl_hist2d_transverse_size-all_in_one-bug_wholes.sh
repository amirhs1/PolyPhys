#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in eps*/; do
    echo "$dir"
    cp probe-trans_foci_cyl_hist2d_transverse_size-all_in_one-bug_wholes*.py probe-trans_foci_cyl_hist2d_transverse_size-all_in_one-bug_wholes.py
    mv probe-trans_foci_bug-trj_wholes_histdd_transverse_size.py "$dir"
    cp -R polyphys "${dir}"
done
