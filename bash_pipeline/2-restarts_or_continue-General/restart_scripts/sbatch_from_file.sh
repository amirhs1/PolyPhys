#!/bin/bash
# This scripts read a list of simulation (which should be restart) line by line
# from a file named "unfinished_run.txt"
while IFS= read -r line; do
    echo "$line"
    cd "$line" || exit
    sbatch submit.sh
    sleep 5
    cd ..
done < unfinished_run.txt