#!/bin/bash
# This scripts go over all the *___-analyze_bug* directories and archive them via tar.
# README: to use this file: connect to the Data Transferring Node on the cluter. For graham, it is gra-dtn1.sharcnet.ca
# you need to submit this scrip with nohup ... &  command.
read -p "report name >" groupname
touch report-${groupname}
report=report-${groupname}.txt
echo $report
for directory in N*-${groupname}/; do
    cp loop_on_dir_tar.sh $directory
    echo ${directory} >> ${report}
    dir=${directory::-1}
    echo $dir
    cd $directory
    nohup bash loop_on_dir_tar.sh &
    cd ..
done