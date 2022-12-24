#!/bin/bash
# This bash scripts compress all the files in all_simulations directories and
# them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move node on your cluter.
# For graham, it is gra-dtn1.sharcnet.ca you need to submit this scrip with
# nohup ... &  command.
read -rp "Enter Project Name > " project
read -rp "Enter Part Number > " part
report="${project}-part${part}-logs-archive_report.txt" #report name
touch "$report"
for dir in D*-logs/; do
    echo "${dir}" >> "${report}"
    project=${dir::-1}.tar.gz
    tar -czvf "${project}" "${dir}" >> "${report}"
done