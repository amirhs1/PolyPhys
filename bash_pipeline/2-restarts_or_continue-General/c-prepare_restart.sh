#!/bin/bash
init_tstep=70000000 # rounded to 5000000
res_tstep=5000000
initial_total_time_step=301000000  # Example: 3*10^8
for incomplete_folder in N*_incomplete; do
    # Extract N* pattern from folder name
    pattern=$(echo "$incomplete_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    res_folder="${pattern}_res"
    # Count the number of .all.lammpstrj.gz files
    jdone=$(ls -l ${incomplete_folder}/${pattern}*.gz | wc -l)
    file_tstep=$((init_tstep + 2*(jdone-1)*res_tstep))
    #ls -l ${incomplete_folder}/restarts/restart.i*.${file_tstep}
    cp ${incomplete_folder}/restarts/restart.i*.${file_tstep} ${res_folder}/
    cp ${incomplete_folder}/submit.sh ${res_folder}/
    jcount=$((24 - jdone))
    #echo $jcount

    # Copy the first 10 lines of input.lmp to restart.lmp in res folder
    temp_file="${res_folder}/input.lmp"
    head -n 11 "${incomplete_folder}/input.lmp" | tail -n 9 > "$temp_file"
    cat "./restart.lmp" >> "$temp_file"

    # Find restart.1.* file
    res_file=$(find "${res_folder}/" -name "restart.i*.${file_tstep}" | cut -d "/" -f 2)
    echo $res_file
    # Add read_restart and variable nloop lines to the top of restart.lmp
    sed -i "1s/^/read_restart $res_file\n/" "${res_folder}/input.lmp"
    sed -i "1s/^/variable nloop equal $jcount\n/" "${res_folder}/input.lmp"

    current_total_time_step=$(($initial_total_time_step - $file_tstep))

    initial_runtime=$(grep "#SBATCH --time=" ${res_folder}/submit.sh | cut -d= -f2)

    initial_days=$(echo $initial_runtime | cut -d- -f1)
    initial_hours=$(echo $initial_runtime | cut -d- -f2 | cut -d: -f1)

    total_initial_hours=$((initial_days * 24 + initial_hours))

    ratio=$(echo "$current_total_time_step / $initial_total_time_step" | bc -l)

    new_total_hours=$(echo "$total_initial_hours * $ratio * 1.30" | bc -l) # Adding 20% margin of error
    new_total_hours_rounded=$(printf "%.0f" $new_total_hours)

    new_days=$(echo "$new_total_hours_rounded / 24" | bc)
    new_hours=$(echo "$new_total_hours_rounded % 24" | bc)

    sed -i "s/#SBATCH --time=$initial_runtime/#SBATCH --time=$(printf "%02d-%02d:00" $new_days $new_hours)/" ${res_folder}/submit.sh

    mkdir ${res_folder}/restarts
done

