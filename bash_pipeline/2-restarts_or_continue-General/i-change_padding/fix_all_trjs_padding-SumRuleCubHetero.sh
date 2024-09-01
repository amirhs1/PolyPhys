#!/bin/bash

# Define the base directory and pattern
#for dir in N*[1-8]; do
#for dir in al*ring; do
for dir in al*linear; do
    cd "$dir"
    ls -l .
    file_count=$(ls ${dir::-5}*.all.lammpstrj | wc -l)
    largest_file=$(ls ${dir::-5}*.all.lammpstrj | sort -t 'j' -k2,2n | tail -1)
    echo rm "$largest_file"
    restart_files=$(ls ${dir::-5}*.all.restart.lammpstrj)
    n_restart_files=$(ls ${dir::-5}*.all.restart.lammpstrj | wc -l)
    echo "Number of restart files: ${n_restart_files}"

    if [ "$n_restart_files" -gt 1 ]; then
        for file in $restart_files; do
            num=$(echo $file | sed -n 's/.*j\([0-9]*\).ring.all.restart.lammpstrj/\1/p')
            echo $num
            num=$((10#$num)) # Convert to decimal
            new_num=$((file_count - 1 + num))
            new_file_name="${dir::-5}.j${new_num}.ring.all.lammpstrj"
            echo mv "$file" "$new_file_name"
        done
    else
        single_file=$(echo $restart_files | head -1)
        new_num=$((file_count))
        echo $new_num
        new_file_name="${dir::-5}.j${new_num}.ring.all.lammpstrj"
        echo mv "$single_file" "$new_file_name"
    fi
    cd ..
done
echo "Finished!"