#!/bin/bash
# Use this command to copy between scratch and project folders:
read -rp "Enter the name of folder to rsync > " dir
nohup tar -zcvf  "${dir}".tar.gz "${dir}" > "${dir}-tar-report.txt" 2>&1 &
echo $! > "${dir}-tar-save_pid.txt"
echo Finished!