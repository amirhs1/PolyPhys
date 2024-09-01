#!/bin/bash
res_tstep=1000000
adump_res_ratio=5
initial_total_time_step=103000000  # Example: 3*10^8
jcont=20
for cont_folder in al*ring_cont; do
    # Extract N* pattern from folder name
    pattern=$(echo "$cont_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    incomplete_folder="${pattern}_incomplete"
    res_folder="${pattern}_cont"
    cp ${pattern}_incomplete/submit.sh ${res_folder}/
    # Copy the first 10 lines of input.lmp to restart.lmp in res folder
    temp_file="${res_folder}/input.lmp"
    head -n 15 "${incomplete_folder}/input.lmp" | tail -n 14 > "$temp_file"
    cat "./continue.lmp" >> "$temp_file"

    # Find restart.1.* or restart.i1.* file
    res_file=$(find "${res_folder}/" -name "*restart_after*" | cut -d "/" -f 2)
    echo "restart file: $res_file"

    # Add read_restart and variable nloop lines to the top of restart.lmp
    sed -i "1s/^/read_restart $res_file\n/" "${res_folder}/input.lmp"
    sed -i "1s/^/variable nloop equal $jcont\n/" "${res_folder}/input.lmp"

    avg_tstep_per_sec=$(awk '/Performance/ {if (seen) {sum += $(NF-1); count++} else seen = 1} END {if (count > 0) print sum / count}' ${pattern}_res/log.lammps)
    avg_tstep_per_sec=$(printf "%.0f" "$avg_tstep_per_sec")
    echo "avg time step per sec: $avg_tstep_per_sec"

    initial_runtime=$(grep "#SBATCH --time=" ${pattern}_incomplete/submit.sh | cut -d= -f2)
    initial_days=$(echo $initial_runtime | cut -d- -f1)
    initial_hours=$(echo $initial_runtime | cut -d- -f2 | cut -d: -f1)
    total_initial_hours=$((initial_days * 24 + initial_hours))
    echo "total time: $total_initial_hours"
    total_initial_hours=$((3 * 24 ))
    echo "new left time: $total_initial_hours"
    new_days=$(echo "$total_initial_hours / 24" | bc)
    new_hours=$(echo "$total_initial_hours % 24" | bc)

    sed -i "s/#SBATCH --time=$initial_runtime/#SBATCH --time=$(printf "%02d-%02d:00" $new_days $new_hours)/" ${res_folder}/submit.sh

    mkdir ${res_folder}/restarts
done

