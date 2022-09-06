#!/bin/bash
# Use this command to copy between scratch and project folders:
read -rp "Enter the name of folder to rsync > " dir
echo "$dir"
nohup tar -cvkzf "${dir}.tar.gz" "${dir}" > "${dir}-tar_report.txt" 2>&1 &
echo $! > "${dir}-tar-save_pid.txt"
