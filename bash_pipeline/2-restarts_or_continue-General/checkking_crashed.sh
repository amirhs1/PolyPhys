#!/bin/bash
while IFS= read -r line; do
    echo "$line"
    #mkdir "${line}_res"
    #mkdir "${line}_res/restarts"
    #cp "${line}/submit.sh" "${line}_res/"
    #cp "${line}/input.lmp" "${line}_res/"
    #cp ./template_loop.lmp "${line}_res"
    #ls "${line}_incomplete"
    #ls "${line}_res"
    #cd "${line}_res"
    ##head -n 14 input.lmp > res.lmp
    #tail -n +16 template_*.lmp >> res.lmp
    #rm input.lmp template_*.lmp
    #mv res.lmp input.lmp
    #head -n 15 "${line}_res/input.lmp"
    #tail -n 25 "${line}_res/input.lmp"
    #head -n 14 "${line}_res/input.lmp"
    #head -n 14 "${line}_res/submit.sh"
    #cd "$line" || exit
    #sbatch "${line}_res/submit.sh"
    #sleep 2
    mv "$line" "${line}"_incomplete
    #cd ..
done < last_dir.txt
