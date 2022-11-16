#!/bin/bash
# This script copies some information from slurm and output files to the screen
# if the simulation finished with exti code 0.

for dir in N*/; do # SumRuleCyl
#for dir in al*/; do # TransFociCub
#for dir in eps*/; do # TransFociCyl
    cd "$dir" || exit
    lastline=$(tail -n 1 slurm*)
    if [[ $lastline == *"exit code 0"* ]]; then
        echo "$dir"
        ls .
        tail -n 50 slurm*
        echo
        tail -n 10 output.txt
    fi
    cd ..
done