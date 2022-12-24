#!/bin/bash
# Fix the log file of a broken simulation with merging that log file with the
# log file of the restared simulation. The total wall time of simulation is
# either lost or equal to that of the restarted simulation.
thermoOne="Step Temp E_pair E_mol TotEng Press"
for dir in N*-logs-broken; do
    cd "$dir" || exit
    for res in N*_res.log; do
        echo "$res"
        log=$(echo "$res" | cut -d _ -f 1)
        dummy=$(awk -v thermoOne="$thermoOne" '$0 ~ thermoOne {print NR}' "${res}" | head -n 1)
        lineNum=$((dummy + 1))
        csplit "$res" "$lineNum"
        rm xx00
        mv xx01 "$log"_rest.log
        firstLine=$(head -n 1 "$log"_rest.log)
        dummy=$(awk -v firstLine="$firstLine" '$0 ~ firstLine {print NR}' "${log}.log" | tail -n 1)
        csplit "${log}.log" "$dummy"
        rm xx01
        cp xx00 "${log}_complete.log"
        cat "${log}_rest.log" >> $"${log}_complete.log"
        rm "${log}_rest.log" xx00 
        # check wall time
        lastline=$(tail -n 1 $"${log}_complete.log")
        totaltime="Total wall time"
        if [[ "${lastline:0:15}" == "${totaltime}" ]]; then
            mv $"${log}_complete.log"  $"${log}_complete_restart_walltime.log"
        else
            echo "Total wall time: 00:00:00" >> $"${log}_complete.log" 
            mv $"${log}_complete.log"  $"${log}_complete_missing_walltime.log" 
        fi
    done
    cd ..
done