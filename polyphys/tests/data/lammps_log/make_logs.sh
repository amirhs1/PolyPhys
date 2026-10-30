#!/bin/sh
# Regenerate the LAMMPS log fixtures: log.<name> from in.<name>, serial run.
#
# The committed logs came from LAMMPS (22 Jul 2025 - Update 4), serial
# build, on 2026-10-30. Rerunning changes timings and CPU figures, so the
# tests read the committed logs and never call LAMMPS.
#
# Run from this directory, with the executable in LMP (default lmp_serial):
#   LMP=lmp_serial sh make_logs.sh
LMP=${LMP:-lmp_serial}
for name in one multi custom yaml; do
  "$LMP" -in "in.$name" -log "log.$name" -screen none || exit 1
done
# in.error stops LAMMPS with an error by design.
"$LMP" -in in.error -log log.error -screen none
# in.killed: stop the run after 5 s, as a walltime limit would.
"$LMP" -in in.killed -log log.killed -screen none &
pid=$!
sleep 5
kill "$pid"
wait "$pid"
# The trajectories, data files and restarts are by-products; keep only logs.
rm -f ./*.lammpstrj ./*.data restart.*
ls -l log.*
