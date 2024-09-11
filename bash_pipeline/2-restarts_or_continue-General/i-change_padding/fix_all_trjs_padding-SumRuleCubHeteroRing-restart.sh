#!/bin/bash
first_round_files=20 # CHECK THIS: Number of all files in the first round
# Define the base directory and pattern
for dir in al*ring; do
    cd "$dir"
    simname=${dir::-5}
    ls -l .
    file_count=$(ls ${simname}*ring.all.lammpstrj.gz | wc -l)
    largest_file=$(ls ${simname}*ring.all.lammpstrj.gz | sort -t 'j' -k2,2n | tail -1)
    rm "$largest_file"
    restart_files=$(ls ${simname}*ring.all.restart.lammpstrj.gz)
    n_restart_files=$(ls ${simname}*ring.all.restart.lammpstrj.gz | wc -l)
    echo "Number of restart files: ${n_restart_files}"

    if [ "$n_restart_files" -gt 1 ]; then
        for file in $restart_files; do
            num=$(echo $file | sed -n 's/.*j\([0-9]*\).ring.all.restart.lammpstrj.gz/\1/p')
            num=$((10#$num)) # Convert to decimal
            new_num=$((file_count - 1 + num))
            new_num=$(printf "%02d" $new_num)
            new_file_name="${simname}.j${new_num}.ring.all.lammpstrj.gz"
            mv "$file" "$new_file_name"
        done
    else
        single_file=$(echo $restart_files | head -1)
        new_num=$((file_count))
        new_num=$(printf "%02d" $new_num)
        new_file_name="${simname}.j${new_num}.ring.all.lammpstrj.gz"
        mv "$single_file" "$new_file_name"
    fi
    cd ..
done
echo "Finished!"
