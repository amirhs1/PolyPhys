#!/bin/bash

# Get the name of the parent directory
parent_dir=$(basename "$(pwd)")

# Define the output CSV file name using the parent directory name
output_csv="${parent_dir}-runtime_slurm-summary.csv"

# Write CSV header
echo "Folder,Total Runtime (seconds),Number of Cores" > "$output_csv"

# Iterate over slurm*out files in N*/ directories
for file in N*/slurm*out; do
    # Extract folder name
    folder_name=$(dirname "$file")

    # Extract start and end times
    start_time=$(grep -oP '\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \w{3} \d{4}' "$file" | head -1)
    potential_end_time=$(grep -oP '\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \w{3} \d{4}' "$file" | tail -1)

    if [ "$potential_end_time" != "$start_time" ]; then
        end_time="$potential_end_time"
    else
        end_time=$(grep -oP '\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$file" | tail -1)
    fi

    # Convert times to seconds since epoch
    start_sec=$(date -d "$start_time" '+%s')
    end_sec=$(date -d "$end_time" '+%s')

    # Calculate duration
    duration=$((end_sec - start_sec))

    # Extract number of cores
    n_cores=$(grep -oP '#SBATCH --ntasks-per-node=\K\d+' "${folder_name}/submit.sh")

    # Append to CSV
    echo "$folder_name,$duration,$n_cores" >> "$output_csv"
done

echo "CSV file created: $output_csv"