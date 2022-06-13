#!/bin/bash

currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
extraction=${name}-extraction
mkdir ${extractionall}
for dir in N*ep*/; do
    mkdir ${extraction}/$dir
    echo $dir
    cd $dir
    mv N*.csv ../${extraction}/$dir
    mv *.txt ../${extraction}/$dir
    rm *.py
    rm -r PipeLine
    cd ..
done