#!/bin/bash
# This bash scripts compress all the files in all_simulations directories and
# them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move node on your cluter.
# For graham, it is gra-dtn1.sharcnet.ca you need to submit this scrip with
# nohup ... &  command.
report=$(pwd | rev | cut -d / -f 1 | rev)-archive_report.txt #report name
touch "$report"
for dir in eps*/; do
    echo "$dir"
    echo "${dir}" >> "${report}"
#    sim=${dir::-1}
#    zipdir=${sim}"-zip"
#    mkdir "${zipdir}"
    gzip -vr "${dir}" >> "${report}"
#    mv ./${dir}*.gz ./${zipdir}
#    project=${zipdir}.tar.gz
#    tar -czvf ${project} ${zipdir} >> "${report}"
done