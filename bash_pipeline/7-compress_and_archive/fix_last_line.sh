#!/bin/bash

# Get the name of the parent directory
parent_dir=$(basename "$(pwd)")

# Define the output CSV file name using the parent directory name
output_csv="${parent_dir}-runtime_slurm-summary.csv"

# Write CSV header
echo "folder,total_runtime_sec,n_cores" > "$output_csv"

# Iterate over slurm*out files in N*/ directories
for incomplete in N*_incomplete; do
    echo $incomplete
    dir=$(echo $incomplete | cut -d '_' -f 1)
    file=$(find "${dir}_incomplete" -maxdepth 1 -name "slurm*" -type f)
    start_time=$(grep -oP 'Starting run at: \K.*' "$file")
#    potential_end_time=$(grep -oP 'Program finished with exit code [0-9]+ at: \K.*|(?<=CANCELLED AT )\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$file" | tail -1)
    potential_end_time=$(grep -oP 'Program finished with exit code [0-9]+ at: \K.*|CANCELLED AT \K\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$file" | tail -1)
    if [ "$potential_end_time" != "$start_time" ]; then
        end_time="$potential_end_time"
    else
        echo "The end time was not found."
    fi

    # Convert times to seconds since epoch
    start_sec=$(date -d "$start_time" '+%s')
    end_sec=$(date -d "$end_time" '+%s')

    # Calculate duration
    duration=$((end_sec - start_sec))
    # Calculate hours, minutes, and seconds
    hours=$((duration / 3600))
    minutes=$(( (duration % 3600) / 60 ))
    seconds=$((duration % 60))
    # Format the duration
    formatted_duration=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
    sed -i '$d' "${dir}/${dir}".incomplete.log
    echo "Total wall time: $formatted_duration" >> ${dir}/${dir}.incomplete.log
    tail -n 2 ${dir}/${dir}.incomplete.log
    cp ${dir}/${dir}.incomplete.log ./${dir}.j01.log
done
