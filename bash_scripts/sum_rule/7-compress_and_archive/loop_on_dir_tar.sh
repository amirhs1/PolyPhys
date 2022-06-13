#!/bin/bash
# This scripts go over the directories in *___-extraction_bug* directory and archive all the directories  in that directory via tar.
# README: to use this file: connect to the Data Transferring Node on the cluter. For graham, it is gra-dtn1.sharcnet.ca
# you need to submit this scrip with nohup ... &  command.
report=$(pwd | rev | cut -d / -f 1 | rev)_archive_report.txt #report name
touch $report
for directory in */; do
    echo ${directory} >> ${report}
    dir=${directory::-1}
    echo $dir
    project=$dir".tar"
    tar -cvf ${project} ${directory} >> ${report}
done