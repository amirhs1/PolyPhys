#!/bin/bash
# Use this command to copy between scratch and project folders:
nohup rsync -axvH --no-g --no-p  "${HOME}/scratch/some_directory" "${HOME}/projects/<project>/some_other_directory" &
