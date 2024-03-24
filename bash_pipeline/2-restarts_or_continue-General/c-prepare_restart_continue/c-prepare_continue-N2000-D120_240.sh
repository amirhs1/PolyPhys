#!/bin/bash
res_tstep=5000000
adump_res_ratio=3
initial_total_time_step=301000000  # Example: 3*10^8
jtotal=20
for incomplete_folder in N*_res; do
    # Extract N* pattern from folder name
    pattern=$(echo "$incomplete_folder" | cut -d "_" -f 1)
    echo $pattern
    # Corresponding res folder
    res_folder="${pattern}_cont"
    # Count the number of .all.lammpstrj.gz files
    cp ${pattern}_incomplete/submit.sh ${res_folder}/
    temp_file="${res_folder}/input.lmp"
    head -n 11 "${incomplete_folder}/input.lmp" | tail -n 9 > "$temp_file"
    cat "./continue.lmp" >> "$temp_file"
    sed -i "1s/^/read_restart restart_after\n/" "${res_folder}/input.lmp"
    avg_tstep_per_sec=$(awk '/Performance/ {if (NR > 1) sum += $(NF-1); count++} END {if (count > 0) print sum / count}' ${incomplete_folder}/log.lammps)
    avg_tstep_per_sec=$(printf "%.0f" "$avg_tstep_per_sec")

    initial_runtime=$(grep "#SBATCH --time=" ${pattern}_incomplete/submit.sh | cut -d= -f2)
    initial_days=$(echo $initial_runtime | cut -d- -f1)
    initial_hours=$(echo $initial_runtime | cut -d- -f2 | cut -d: -f1)
    total_initial_hours=$((initial_days * (jtotal + 1) + initial_hours))
    echo "total time: $total_initial_hours"
    
    left_tsteps=1000000
    left_runtime=$(echo "$left_tsteps / $avg_tstep_per_sec" | bc -l)
    left_runtime=$(printf "%.0f" "$left_runtime")
    left_runtime_hr=$(echo "$left_runtime / 3600" | bc -l)
    echo "left time: $left_runtime_hr"

    ratio=3

    new_total_hours=$(echo "$left_runtime_hr * $ratio" | bc -l) # Adding 20% margin of error
    new_total_hours_rounded=$(printf "%.0f" $new_total_hours)

    echo "new left time: $new_total_hours_rounded"
    new_days=$(echo "$new_total_hours_rounded / 24" | bc)
    echo $new_day
    new_hours=$(echo "$new_total_hours_rounded % 24" | bc)
    echo $new_hours
    sed -i "s/#SBATCH --time=$initial_runtime/#SBATCH --time=$(printf "%02d-%02d:00" $new_days $new_hours)/" ${res_folder}/submit.sh

#    mkdir ${res_folder}/restarts
done

