#!/bin/bash
# script to generate folders with necessary files for a job array on the slurm schedular
for dir in N*/; do
        echo $dir

        cp extraction_segments_hists_mon.py $dir
        cp __init__.py $dir
        cp PipeLine.py $dir

        cd $dir
        mkdir PipeLine
        mv __init__.py PipeLine
        mv PipeLine.py PipeLine
        cd ..
done