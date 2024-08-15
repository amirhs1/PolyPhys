#!/bin/bash

for dir in N*[1-8]; do
for dir in al*ring; do
    cd $dir
    rm ${dir}*bug*lammpstrj
    mv test/${dir}*lammpstrj .
    rm -r test
    cd ..
done