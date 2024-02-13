#!/bin/bash

# Get the name of the parent directory
parent_dir=$(basename "$(pwd)")

# Define the output CSV file name using the parent directory name
output_csv="${parent_dir}-runtime_slurm-summary.csv"

# Write CSV header
echo "folder,total_runtime_sec,n_cores" > "$output_csv"

# Iterate over slurm*out files in N*/ directories
for file in N*/slurm*out; do
    # Extract folder name
    folder_name=$(dirname "$file")

    # Extract start and end times
    start_time=$(grep -oP 'Starting run at: \K.*' "$file")
    potential_end_time=$(grep -oP 'Program finished with exit code [0-9]+ at: \K.*|(?<=CANCELLED AT )\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$file" | tail -1)

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

    # Extract number of cores
    n_cores=$(grep -oP '#SBATCH --ntasks-per-node=\K\d+' "${folder_name}/submit.sh")

    # Append to CSV
    echo "$folder_name,$duration,$n_cores" >> "$output_csv"
done

echo "CSV file created: $output_csv"

