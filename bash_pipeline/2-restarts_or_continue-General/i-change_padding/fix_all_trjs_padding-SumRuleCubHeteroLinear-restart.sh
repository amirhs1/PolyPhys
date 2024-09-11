#!/bin/bash
first_round_files=20 # CHECK THIS: Number of all files in the first round
# Define the base directory and pattern
for dir in al*linear; do
    cd "$dir"
    simname=${dir::-7}
    ls -l .
    file_count=$(ls ${simname}*linear.all.lammpstrj | wc -l)
    largest_file=$(ls ${simname}*linear.all.lammpstrj | sort -t 'j' -k2,2n | tail -1)
    echo rm "$largest_file"
    restart_files=$(ls ${simname}*linear.all.restart.lammpstrj)
    n_restart_files=$(ls ${simname}*linear.all.restart.lammpstrj | wc -l)
    echo "Number of restart files: ${n_restart_files}"

    if [ "$n_restart_files" -gt 1 ]; then
        for file in $restart_files; do
            num=$(echo $file | sed -n 's/.*j\([0-9]*\).linear.all.restart.lammpstrj/\1/p')
            num=$((10#$num)) # Convert to decimal
            new_num=$((file_count - 1 + num))
            new_num=$(printf "%02d" $new_num)
            new_file_name="${simname}.j${new_num}.linear.all.lammpstrj"
            echo mv "$file" "$new_file_name"
        done
    else
        single_file=$(echo $restart_files | head -1)
        new_num=$((file_count))
        new_num=$(printf "%02d" $new_num)
        new_file_name="${simname}.j${new_num}.linear.all.lammpstrj"
        echo mv "$single_file" "$new_file_name"
    fi
    cd ..
done
echo "Finished!"
