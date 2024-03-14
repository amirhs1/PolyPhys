#!/bin/bash

output_file="performance_report.txt"
> "$output_file" # Clear the file or create it if it doesn't exist

convert_to_hours() {
    local time_str=$1
    local days=$(echo $time_str | cut -d'-' -f1)
    local hours=$(echo $time_str | cut -d'-' -f2 | cut -d':' -f1)
    local minutes=$(echo $time_str | cut -d':' -f2)
    # local seconds=$(echo $time_str | cut -d':' -f3)
    echo $((days*24 + hours + minutes/60))
    # + seconds/3600))
}

for folder in N*ens1; do
    echo "Processing folder: $folder" >>  "$output_file"
    filename=${folder}/log.lammps
    submit_file=${folder}/submit.sh
    #echo $filename >> "$output_file"

    # Process submit.sh file
    if [ -f "$submit_file" ]; then
        sbatch_line=$(grep "^#SBATCH --time=" "$submit_file")
        echo "SBATCH Line: $sbatch_line" >> "$output_file"

        # Extract and convert runtime to hours
        runtime_str=$(echo $sbatch_line | awk -F"time=" '{print $2}')
        runtime_hours=$(convert_to_hours "$runtime_str")
        echo "SLURM runtime in hours: $runtime_hours" >> "$output_file"
    else
        echo "submit.sh not found in $folder" >> "$output_file"
    fi

    total=0
    count=0

    # Process every Performance line
    while IFS= read -r line; do
        value=$(echo $line | awk '{print $(NF-1)}')
        if [[ ! -z "$value" ]]; then
            #echo "Extracted value: $value" >> "$output_file"
            total=$(echo "$total + $value" | bc)
            count=$((count + 1))
        fi
    done < <(awk '/^Performance:/ {print}' "$filename")

    if [ $count -gt 0 ]; then
        average=$(echo "$total / $count" | bc -l)
        echo "Average for $folder: $average timesteps/s" >> "$output_file"
        result=$(echo "230000000 / ($average * 3600)" | bc -l)
        result_rounded=$(printf "%.0f" "$result")
        echo "Estimated time for $folder: $result_rounded hours" >> "$output_file"

        # Compare with SLURM runtime
        if [ ! -z "$runtime_hours" ] && [ $result_rounded -gt $runtime_hours ]; then
            echo "Warning: Estimated time ($result_rounded hours) is greater than SLURM runtime ($runtime_hours hours)" >> "$output_file"
        fi
    else
        echo "Insufficient performance data in $folder." >> "$output_file"
    fi

    echo >> "$output_file"
done

