#!/bin/bash
# This bash scripts compress all the files in all simulations directories and them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move NODe on your cluter. For graham, it is gra-dtn1.sharcnet.ca
# you need to submit this scrip with nohup ... &  command.
for dir in N*[1-8]/; do
    file=${dir::-1}
    echo $file
    cd $dir
    mv restarts restarts_old
    cd restarts_old
    ls 
    mkdir restarts_zip
    mv *.gz restarts_zip
    mv restarts_zip ..
    cd ..
    zipdir=${file}-zip
    mv restarts_zip ../${zipdir}
    cd ..
    ls $zipdir
done