#!/bin/bash

# Directory where the incomplete folders will be moved
destination_dir="incomplete"

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Read each line from the file
while IFS= read -r folder; do
    # Check if the folder exists
    if [[ -d "$folder" ]]; then
        # Move the folder to the destination directory
        mv "$folder" "$destination_dir/"
    else
        echo "Folder not found: $folder"
    fi
done < folders_to_check.txt
