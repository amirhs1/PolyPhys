#!/bin/bash

# Get the name of the parent directory
parent_dir=$(basename "$(pwd)")

# Define the output CSV file name using the parent directory name
output_csv="allInOne-${parent_dir}-runtime_summary.csv"

# Write CSV header
echo "log_name,wall_time_sec,wall_time_hr,n_cores,n_atoms,avg_step_per_sec,bin_scheme,rskin" > "$output_csv"

# Iterate over log files in specified directories
for file in ns*logs/al*ring.log; do
    # Extract the necessary details using awk and grep
    wall_time=$(grep "Total wall time:" "$file" | awk '{split($4,t,":"); print t[1]*3600 + t[2]*60 + t[3]}')
    wall_time_hr=$(awk -v seconds="$wall_time" 'BEGIN {print seconds / 3600}')
    read n_cores n_atoms <<< $(grep -m 1 "Loop time of" "$file" | awk '{print $6, $12}')
    avg_tstep_per_sec=$(awk '/Performance/ {sum += $(NF-3); count++} END {if (count > 0) print sum / count}' "$file")
    
    # Find the correct line with "neighbor" that has an actual value and exactly 3 fields
    read rskin bin_scheme <<< $(awk '/neighbor/ && NF==3 && $2 !~ /^\$/ {print $2, $3; exit}' "$file")
    
    # Check if neighbor details were found
    if [ -z "$rskin" ] || [ -z "$bin_scheme" ]; then
        echo "Error: Neighbor details not found in the expected format (NF==3 with numeric value) in file $file" >&2
        rskin="N/A"
        bin_scheme="N/A"
    fi

    # Append to CSV
    echo "$(basename "$file"),$wall_time,$wall_time_hr,$n_cores,$n_atoms,$avg_tstep_per_sec,$bin_scheme,$rskin" >> "$output_csv"

done

echo "CSV file created: $output_csv"

