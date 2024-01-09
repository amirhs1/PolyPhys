#!/bin/bash

# Define the base directory and pattern
for dir in N*[1-8]; do
    cd "$dir"
    ls -l .
    file_count=$(ls ${dir}*.all.lammpstrj.gz | wc -l)
    largest_file=$(ls ${dir}*.all.lammpstrj.gz | sort -t 'j' -k2,2n | tail -1)
    rm "$largest_file"
    restart_files=$(ls ${dir}*.all.restart.lammpstrj.gz)
    n_restart_files=$(ls ${dir}*.all.restart.lammpstrj.gz | wc -l)
    echo "Number of restart files: ${n_restart_files}"

    if [ "$n_restart_files" -gt 1 ]; then
        for file in $restart_files; do
            num=$(echo $file | sed -n 's/.*j\([0-9]*\).all.restart.lammpstrj.gz/\1/p')
            num=$((10#$num)) # Convert to decimal
            new_num=$((14 + file_count - 1 + num))
            new_file_name="${dir}.j${new_num}.all.lammpstrj.gz"
            mv "$file" "$new_file_name"
        done
    else
        single_file=$(echo $restart_files | head -1)
        new_num=$((14 + file_count))
        new_file_name="${dir}.j${new_num}.all.lammpstrj.gz"
        mv "$single_file" "$new_file_name"
    fi
    cd ..
done