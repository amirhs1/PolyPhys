#!/bin/bash
# Initialize a counter
counter=1
# Loop through files and modify them
cd "data_files"
for file in N2000epsilon*r*; do
    # New file name
    new_file="file_${counter}.data"

    # Get the first line of the file
    first_line=$(head -n 1 "$file")

    # Add the original file name at the end of the first line
    # and then add the rest of the file
    (echo "$first_line $file"; tail -n +2 "$file") > "${new_file}"

    # Temporary file
    temp_file="temp_data_file.data"

    # Delete the line containing "Masses" and the next three lines
    sed '/^Masses$/{N;N;N;d;}' $new_file > $temp_file

    # Change "1 atom types" to "2 atom types"
    # We use awk to search for the line containing "atom types" and change the number preceding it
    awk '{if ($2 == "atom" && $3 == "types") $1="2"; print}' $temp_file > $new_file

    # Clean up the temporary file
    rm $temp_file

    # Increment the counter
    ((counter++))
done
cd ..