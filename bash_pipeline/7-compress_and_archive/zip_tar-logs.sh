#!/bin/bash
# This bash scripts compress all the files in all_simulations directories and
# them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move node on your cluter.
# For graham, it is gra-dtn1.sharcnet.ca you need to submit this scrip with
# nohup ... &  command.
report=$(pwd | rev | cut -d / -f 1 | rev)-archive_report.txt #report name
touch "$report"
echo "Start archiving..."
for dir in ns*-logs/; do
#for dir in N*-logs/; do
    sim=$(echo "$dir" | cut -d / -f 1)
    echo "$sim"
    echo "${dir}" >> "${report}"
    tar -cvkzf "${sim}.tar.gz" "${dir}"
done
echo "Finished!"