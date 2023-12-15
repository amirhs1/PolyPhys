#!/bin/bash
# CAUTION: Set the directory pattern before using
# Define the base directory and pattern
for dirres in N*ens1_res; do
    dir="${dirres::-4}"
    cd $dir
    file_count=$(ls ${dir}*.all.lammpstrj.gz | wc -l)
    largest_file=$(ls ${dir}*.all.lammpstrj.gz | sort -t 'j' -k2,2n | tail -1)
    rm $largest_file
    mv ${dir}.all.continue.lammpstrj.gz ${dir}.j21.all.lammpstrj.gz

    restart_files=$(ls ${dir}*.all.restart.lammpstrj.gz)
    for file in $restart_files; do
        num=$(echo $file | sed -n 's/.*j\([0-9]*\).all.restart.lammpstrj.gz/\1/p')
        num=$((10#$num)) # Convert to decimal
        new_num=$((file_count - 1 + num))
        new_file_name="${dir}.j${new_num}.all.lammpstrj.gz"
        mv $file $new_file_name
    done
    cd ..
done