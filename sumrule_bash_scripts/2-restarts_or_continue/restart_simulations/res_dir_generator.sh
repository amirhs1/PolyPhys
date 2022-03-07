#!/bin/bash
for dir in N*/;do
    file=$(echo "$dir" | cut -d / -f 1)
    echo "$file"
    resdir=${file}_res
    mkdir "$resdir"
    indir=${file}_incomplete
    mv "$dir" "$indir"
    mkdir "$file"
    cp input.lmp "$resdir"
    cp submit.sh "$resdir"
done