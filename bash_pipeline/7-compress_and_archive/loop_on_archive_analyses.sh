#!/bin/bash
#for dir in N*-probe; do
#for dir in ns*probe; do
for dir in ns*/; do
    name=${dir:0:${#dir}-1}
    echo $name
    tar -zcvf  "${name}".tar.gz "${dir}" > "${name}-tar_report.txt"
    echo "finished!" > "${name}-tar_report.txt"
done
