#!/bin/bash

# Get the name of the parent directory
parent_dir=$(basename "$(pwd)")

# Define the output CSV file name using the parent directory name
output_csv="allInOne-${parent_dir}-runtime_summary.csv"

# Write CSV header
echo "log_name,wall_time_slurm_sec,wall_time_sec,wall_time_hr,n_cores,n_atoms,avg_step_per_sec,bin_scheme,rskin" > "$output_csv"

# Iterate over log files in specified directories
#for slurm in al*[1-8].{ring,ring_res,ring_incomplete}/slurm*.out; do # TransFociCub, SumRuleCubHeteroRing
for slurm in al*[1-8].{linear,linear_res,linear_incomplete}/slurm*.out; do # SumRuleCubHeteroLinear
#for slurm in N*{ring,ring_res,ring_incomplete}/slurm*.out; do #HnsCub
#for slurm in N*{[1-8],[1-8]_res,[1-8]_incomplete}/slurm*.out; do #SumRuleCul
#for slurm in eps*{ring,ring_res,ring_incomplete}/slurm*.out; do #TransFociCyl
    # Extract folder name
    folder_name=$(dirname "$slurm")
    file="${dirname}/log.lammps"
    # SLURM out file
    start_time=$(grep -oP '\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \w{3} \d{4}' "$slurm" | head -1)
    potential_end_time=$(grep -oP '\w{3} \w{3} [ \d]\d \d\d:\d\d:\d\d \w{3} \d{4}' "$slurm" | tail -1)

    if [ "$potential_end_time" != "$start_time" ]; then
        end_time="$potential_end_time"
    else
        end_time=$(grep -oP '\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}' "$slurm" | tail -1)
    fi
    start_sec=$(date -d "$start_time" '+%s')
    end_sec=$(date -d "$end_time" '+%s')
    wall_time_slurm=$((end_sec - start_sec))
    # LAMMPS log file
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
    echo "$folder_name,$wall_time_slurm,$wall_time,$wall_time_hr,$n_cores,$n_atoms,$avg_tstep_per_sec,$bin_scheme,$rskin" >> "$output_csv"

done

echo "CSV file created: $output_csv"

