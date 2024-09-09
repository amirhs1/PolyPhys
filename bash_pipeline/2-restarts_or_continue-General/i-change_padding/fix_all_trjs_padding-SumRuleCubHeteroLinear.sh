#!/bin/bash
first_round_files=20 # CHECK THIS: Number of all files in the first round
# Define the base directory and pattern
#for dir in N*[1-8]; do
for dir in al*linear; do
    cd "$dir"
    simname=${dir::-7}
    ls -l .
    file_count=$(ls ${simname}*.linear.all.continue.lammpstrj | wc -l)
    largest_file=$(ls ${simname}*.linear.all.continue.lammpstrj | sort -t 'j' -k2,2n | tail -1)
    cont_files=$(ls ${simname}*.linear.all.continue.lammpstrj)
    echo "Number of restart files: ${n_restart_files}"

    if [ "$file_count" -gt 1 ]; then
        for file in $cont_files; do
            num=$(echo $file | sed -n 's/.*j\([0-9]*\).linear.all.continue.lammpstrj/\1/p')
            num=$((10#$num)) # Convert to decimal
            new_num=$((first_round_files + num))
            new_file_name="${simname}.j${new_num}.linear.all.lammpstrj"
            mv "$file" "$new_file_name"
        done
    else
        single_file=$(echo $restart_files | head -1)
        new_num=$((first_round_files + file_count))
        new_file_name="${simname}.j${new_num}.linear.all.lammpstrj"
        mv "$single_file" "$new_file_name"
    fi
    cd ..
done
echo "Finished!"
