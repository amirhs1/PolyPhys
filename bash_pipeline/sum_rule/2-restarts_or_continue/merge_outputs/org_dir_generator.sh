#!/bin/bash
for dir in N*_res_after_cont/;do
    comdir=$(echo "$dir" | cut -d _ -f 1)
    echo "$comdir"
    echo mv "$comdir" "${comdir}_incomplete"
    echo mv "$dir" "${comdir}_res"
    echo mkdir "$comdir"
done