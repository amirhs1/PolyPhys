#!/bin/bash
# Use this command to copy between scratch and project folders:
read -rp "Enter the name of folder to rsync > " dir
echo "$dir"
#nohup rsync -axvH --no-g --no-p  "${HOME}/scratch/${dir}" "${HOME}/amirhsi_rrg/trans_foci_cylinder" > "${dir}-rsync_report.txt" 2>&1 &
nohup rsync -axvH --no-g --no-p  "${HOME}/scratch/${dir}" "${HOME}/amirhsi_rrg/sum_rule_cylinder" > "${dir}-rsync_report.txt" 2>&1 &
echo $! > "${dir}-rsync-save_pid.txt"