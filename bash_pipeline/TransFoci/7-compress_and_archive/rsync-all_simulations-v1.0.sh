#!/bin/bash
# Use this command to copy between scratch and project folders:
read -rp "Enter the name of folder to rsync > " dir
project=$(echo "$dir" | cut -d - -f 1)
echo "$project"
nohup rsync -axvH --no-g --no-p  "${HOME}/scratch/${dir}" "${HOME}/amirhsi_rrg/trans_foci_cylinder" > "${project}-rsync_report.txt" 2>&1 &
echo $! > "${project}-rsync-save_pid.txt"