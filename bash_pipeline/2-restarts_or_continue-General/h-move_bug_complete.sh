#!/bin/bash

for dir in N*[1-8]; do
    cd $dir
    rm ${dir}*bug*lammpstrj
    mv test/${dir}*lammpstrj .
    rm -r test
    cd ..
done