#!/bin/bash
# CHECK the lines have CHECK, they deffier from one project to another.
# CHECK neighbor setting in the restart template.
res_tstep=5000000 # CHECK 
adump_res_ratio=3 # CHECk 
initial_total_time_step=301000000  # CHECK
jtotal=20 # CHECK
for incomplete_folder in N*_incomplete; do # CHECK
    # Extract N* pattern from folder name
    pattern=$(echo "$incomplete_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    res_folder="${pattern}_res"
    # Count the number of .all.lammpstrj.gz files
    jdone=$(ls -l ${incomplete_folder}/${pattern}*.gz | wc -l)
    file_tstep=$((adump_res_ratio*(jdone-1)*res_tstep))
    ls -l ${incomplete_folder}/restarts/restart.*.${file_tstep}
    cp ${incomplete_folder}/restarts/restart.*.${file_tstep} ${res_folder}/
    cp ${incomplete_folder}/submit.sh ${res_folder}/
    jcount=$((jtotal + 1 - jdone))
    echo $jcount
    # Copy the first 10 lines of input.lmp to restart.lmp in res folder
    temp_file="${res_folder}/input.lmp"
    head -n 11 "${incomplete_folder}/input.lmp" | tail -n 9 > "$temp_file" # CHECK: 11 & 9? 
    cat "./restart.lmp" >> "$temp_file" # CHECK: restart file contents

    # Find restart.1.* or restart.i1.* file
    res_file=$(find "${res_folder}/" -name "restart.i*.${file_tstep}" | cut -d "/" -f 2)  # CHECK: the restart.i* or restart.*?
    # Check if res_file is empty
    if [[ -z "$res_file" ]]; then
    # If the first search didn't find anything, search with the alternate pattern
    res_file=$(find "${res_folder}/" -name "restart.*.${file_tstep}" | cut -d "/" -f 2)
    fi
    echo "restart file: $res_file"

    # Add read_restart and variable nloop lines to the top of restart.lmp
    sed -i "1s/^/read_restart $res_file\n/" "${res_folder}/input.lmp"
    sed -i "1s/^/variable nloop equal $jcount\n/" "${res_folder}/input.lmp"

    current_total_time_step=$(grep -B 3 "Performance" ${incomplete_folder}/log.lammps | tail -n 4 | head -n 1 | awk '{print $1}')
    echo "current time step: $current_total_time_step"

    avg_tstep_per_sec=$(awk '/Performance/ {if (seen) {sum += $(NF-1); count++} else seen = 1} END {if (count > 0) print sum / count}' ${incomplete_folder}/log.lammps)
    avg_tstep_per_sec=$(printf "%.0f" "$avg_tstep_per_sec")

    initial_runtime=$(grep "#SBATCH --time=" ${res_folder}/submit.sh | cut -d= -f2)
    initial_days=$(echo $initial_runtime | cut -d- -f1)
    initial_hours=$(echo $initial_runtime | cut -d- -f2 | cut -d: -f1)
    total_initial_hours=$((initial_days * (jtotal + 1) + initial_hours))
    echo "total time: $total_initial_hours"
    
    left_tsteps=$((initial_total_time_step - current_total_time_step))
    left_runtime=$(echo "$left_tsteps / $avg_tstep_per_sec" | bc -l)
    left_runtime=$(printf "%.0f" "$left_runtime")
    left_runtime_hr=$(echo "$left_runtime / 3600" | bc -l)
    echo "left time: $left_runtime_hr"

    ratio=1.2

    new_total_hours=$(echo "$left_runtime_hr * $ratio " | bc -l) # Adding 20% margin of error
    new_total_hours_rounded=$(printf "%.0f" $new_total_hours)

    echo "new left time: $new_total_hours_rounded"
    new_days=$(echo "$new_total_hours_rounded / 24" | bc)
    new_hours=$(echo "$new_total_hours_rounded % 24" | bc)

    sed -i "s/#SBATCH --time=$initial_runtime/#SBATCH --time=$(printf "%02d-%02d:00" $new_days $new_hours)/" ${res_folder}/submit.sh

    mkdir ${res_folder}/restarts
done
echo "Finished!"