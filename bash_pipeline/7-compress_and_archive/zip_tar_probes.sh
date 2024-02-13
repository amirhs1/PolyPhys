#!/bin/bash
for dir in N*ring; do
    tar -zcvf  "${dir}".tar.gz "${dir}"
    echo finished!
done