#!/bin/bash
# incomfile (the main trj file) tstep (timestep ast which the merge is done) trjfile (the restart trj file)
thermoOne="Step Temp E_pair E_mol TotEng Press"
for dir in N*-logs-broken-different/; do
    cd "$dir" || exit
    for res in N*.j03.log; do
        echo "$res"
        # trunk-ignore(shellcheck/SC2086)
        log=$(echo $res | cut -d . -f 1-6)
        dummy=$(awk -v thermoOne="$thermoOne" '$0 ~ thermoOne {print NR}' "${res}" | head -n 1)
        lineNum=$((dummy + 1))
        csplit "$res" "$lineNum"
        rm xx00
        mv xx01 "$log"_rest.log
        firstLine=$(head -n 1 "$log"_rest.log)
        dummy=$(awk -v firstLine="$firstLine" '$0 ~ firstLine {print NR}' "${log}.j02.log" | tail -n 1)
        csplit "${log}.j02.log" "$dummy"
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