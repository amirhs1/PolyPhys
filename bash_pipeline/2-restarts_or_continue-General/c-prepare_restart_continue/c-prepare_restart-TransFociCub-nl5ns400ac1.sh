#!/bin/bash
initial_tstep=3000000
res_tstep=1000000
adump_res_ratio=5
initial_total_time_step=103000000  # Example: 3*10^8
jtotal=20
for incomplete_folder in al*_incomplete; do
    # Extract N* pattern from folder name
    pattern=$(echo "$incomplete_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    res_folder="${pattern}_res"
    # Count the number of .all.lammpstrj.gz files
    jdone=$(ls -l ${incomplete_folder}/${pattern::-5}*.gz | wc -l)
    file_tstep=$((initial_tstep+adump_res_ratio*(jdone-1)*res_tstep))
    ls -l ${incomplete_folder}/restarts/restart.*.${file_tstep}
    cp ${incomplete_folder}/restarts/restart.*.${file_tstep} ${res_folder}/
    cp ${incomplete_folder}/submit.sh ${res_folder}/
    jcount=$((jtotal + 1 - jdone))
    echo $jcount
    # Copy the first 10 lines of input.lmp to restart.lmp in res folder
    temp_file="${res_folder}/input.lmp"
    head -n 13 "${incomplete_folder}/input.lmp" | tail -n 12 > "$temp_file"
    cat "./restart.lmp" >> "$temp_file"

    # Find restart.1.* or restart.i1.* file
    res_file=$(find "${res_folder}/" -name "restart.i*.${file_tstep}" | cut -d "/" -f 2)
    # Check if res_file is empty
    if [[ -z "$res_file" ]]; then
    # If the first search didn't find anything, search with the alternate pattern
    res_file=$(find "${res_folder}/" -name "restart.*.${file_tstep}" | cut -d "/" -f 2)
    fi
    echo "restart file: $res_file"

    # Add read_restart and variable nloop lines to the top of restart.lmp
    sed -i "1s/^/read_restart $res_file\n/" "${res_folder}/input.lmp"
    sed -i "1s/^/variable nloop equal $jcount\n/" "${res_folder}/input.lmp"

    current_total_time_step=$(grep -B 10 "Performance" ${incomplete_folder}/log.lammps | tail -n 8 | head -n 1 | grep "Step" | awk '{print $3}')
    echo "current time step: $current_total_time_step"

    avg_tstep_per_sec=$(awk '/Performance/ {if (seen) {sum += $(NF-3); count++} else seen = 1} END {if (count > 0) print sum / count}' ${incomplete_folder}/log.lammps)
    avg_tstep_per_sec=$(printf "%.0f" "$avg_tstep_per_sec")
    echo "avg time step per sec: $avg_tstep_per_sec"

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

    ratio=1

    new_total_hours=$(echo "$left_runtime_hr * $ratio * 1.30" | bc -l) # Adding 20% margin of error
    new_total_hours_rounded=$(printf "%.0f" $new_total_hours)

    echo "new left time: $new_total_hours_rounded"
    new_days=$(echo "$new_total_hours_rounded / 24" | bc)
    new_hours=$(echo "$new_total_hours_rounded % 24" | bc)

    sed -i "s/#SBATCH --time=$initial_runtime/#SBATCH --time=$(printf "%02d-%02d:00" $new_days $new_hours)/" ${res_folder}/submit.sh

    mkdir ${res_folder}/restarts
done

