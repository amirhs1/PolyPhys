#!/bin/bash
for dir in N*/; do
    cd ${dir}
    sbatch submit.sh
    sleep 5
    cd ..
done
