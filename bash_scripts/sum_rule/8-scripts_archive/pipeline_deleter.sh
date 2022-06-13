#!/bin/bash

# script to generate folders with necessary files for a job array on the slurm schedular

for dir in N*/; do
        echo $dir
        
        cd $dir
        rm -r PipeLine
        rm *.py
        rm *.txt
        cd ..
done