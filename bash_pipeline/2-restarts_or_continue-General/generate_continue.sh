#!/bin/bash
# CUATION: check folder pattern and namings

for incomplete_folder in N*_res; do
    # Extract N* pattern from folder name
    pattern=$(echo "$incomplete_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    res_folder="${pattern}_cont"
    # Copy restart.lmp to N*_res folder
    cp "./restart.lmp" "${res_folder}/"

    # Copy the first 10 lines of input.lmp to restart.lmp in res folder
    temp_file="${res_folder}/temp_restart.lmp"
    head -n 10 "${res_folder}/input.lmp" > "$temp_file"
    cat "${res_folder}/restart.lmp" >> "$temp_file"
    mv "$temp_file" "${res_folder}/restart.lmp"

    # Find restart.1.* file
    res_file=$(find "${res_folder}/" -name "restart.1.*" | cut -d "/" -f 2)
    echo $res_file
    # Add read_restart and variable nloop lines to the top of restart.lmp
    sed -i "1s/^/read_restart $res_file\n/" "${res_folder}/restart.lmp"
    sed -i "1s/^/variable nloop equal $jcount\n/" "${res_folder}/restart.lmp"
    mkdir ${res_folder}/restarts
done