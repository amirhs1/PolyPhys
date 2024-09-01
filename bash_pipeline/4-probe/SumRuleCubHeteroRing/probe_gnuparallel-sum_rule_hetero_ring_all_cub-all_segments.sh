#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in al*ring/; do
    echo "$dir"
    cp probe-sum_rule_hetero_ring_all_cub-all_segments.py  ./"$dir"
    cp -R ./polyphys ./"${dir}"
done
