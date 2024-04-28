#!/bin/bash

# Constants for the simulation
freq=5000
init=1000000
total=301000000
line_per_tstep=2009

# Calculate expected number of frames
n_frames=$((((total - init) / freq) + 1))
echo "Total frames: $n_frames"

# Read directory names from "missing_runs.txt"
while read dir; do
  # Ensure directory ends with a slash
  dir=${dir%/}
  echo $dir
  # Extract the last line from the accompanying text file which has the FLAG timestep
  flag_tstep=$(tail -n 1 "${dir}-missing_timestep.txt")
  echo "Flag timestep: $flag_tstep"

  # Calculate frame number for the flag timestep
  flag_frame=$(((flag_tstep - init) / freq))
  cd $dir
  # Calculate the number of lines and frames in the broken dump file
  trj_broken_lines=$(wc -l < "${dir}.bug.lammpstrj")
  trj_broken_frames=$((trj_broken_lines / line_per_tstep))
  echo "Broken trj frames $trj_broken_frames"

  # Calculate the number of missing frames and lines
  n_missing_frames=$((n_frames - trj_broken_frames))
  n_missing_lines=$((n_missing_frames * line_per_tstep))
  echo "Number of missing frames: $n_missing_frames"
  # Find the exact line number in the dump file for the flag frame
  flag_line=$(((flag_frame - n_missing_frames + 1) * line_per_tstep ))
  # Check whether the flag timestep is correctly identified in the dump file
  flag_frame_in_dump=$(head -n $flag_line "${dir}.bug.lammpstrj" | tail -n $line_per_tstep | head -n 2 | tail -n 1)
  echo "Flag timestep in file: $flag_frame_in_dump"
  missing_tstep_end=$((flag_frame_in_dump - freq))
  # Calculate and find the timestep before the flag timeste
  tstep_before_flag=$(head -n "$flag_line" "${dir}.bug.lammpstrj" | tail -n $(( 2 * line_per_tstep)) | head -n 2 | tail -n 1)
  missing_tstep_start=$((tstep_before_flag + freq))
  echo "Timestep before flag timestep: $tstep_before_flag"
  cd ..
  touch ${dir}-missing_tsteps-${missing_tstep_end}_${missing_tstep_end}-inclusive.txt
done < "missing_runs.txt"