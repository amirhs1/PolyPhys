#!/bin/bash

# Define the file name
FILE="slurm_jobs.txt"

# Read each line from the file
while IFS= read -r line; do
    # Check if the line contains 'Submitted batch job'
    if [[ $line == *"Submitted batch job"* ]]; then
        # Extract the job ID and store in variable
        job_id=$(echo $line | awk '{print $4}')
        
        # Execute the cancel command with the job ID
        scancel $job_id
    fi
done < "$FILE"