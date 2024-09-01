#!/bin/bash
# copies the python script and package which are then executed by gnuparallal.
for dir in al*linear/; do
    echo "$dir"
    cp probe-sum_rule_hetero_linear-all_in_one-bug_wholes.py  ./"$dir"
    cp -R polyphys ./"${dir}"
done
