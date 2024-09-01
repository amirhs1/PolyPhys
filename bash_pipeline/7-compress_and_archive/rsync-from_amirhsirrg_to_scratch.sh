#!/bin/bash
# Use this command to copy between scratch and project folders:
read -rp "Enter the name of folder to rsync > " dir
project=$(echo "$dir" | cut -d / -f 1)
echo "$project"
nohup rsync -axvH --no-g --no-p  "${HOME}/amirhsi_rrg/${dir}" "${HOME}/scratch" > "${project}-rsync_report.txt" 2>&1 &
echo $! > "${project}-rsync-save_pid.txt"