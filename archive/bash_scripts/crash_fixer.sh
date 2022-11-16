#!/bin/bash
# show the name of an incomplete simulation on the screen.

for dir in N*/; do # SumRuleCyl
#for dir in al*/; do # TransFociCub
#for dir in eps*/; do # TransFociCyl
    cd "$dir" || exit
    lastline=$(tail -n 1 slurm*)
    if [[ $lastline != *"exit code 0"* ]]; then
        echo "$dir"
    fi
    cd ..
done