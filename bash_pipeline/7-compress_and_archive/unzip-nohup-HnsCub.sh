#!/bin/bash
# Use this command to copy between scratch and project folders:
for dir in N*gz;do
    nohup tar -xvf "${dir}" > "${dir::-7}-rsync_report.txt" 2>&1 &
    echo $! > "${dir::-7}-rsync-save_pid.txt"
done