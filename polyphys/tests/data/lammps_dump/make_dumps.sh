#!/bin/sh
# Regenerate the LAMMPS dump and data fixtures from the in.* scripts, serial
# run. Each script's header lists the files it writes.
#
# The committed files came from LAMMPS (22 Jul 2025 - Update 4), serial
# build with gzip support, on 2026-10-31. The tests read the committed files
# and never call LAMMPS.
#
# Run from this directory, with the executable in LMP (default lmp_serial):
#   LMP=lmp_serial sh make_dumps.sh
LMP=${LMP:-lmp_serial}
for style in bond angle molecular full; do
  "$LMP" -in in.ring -var style "$style" -log none -screen none || exit 1
done
"$LMP" -in in.triclinic -log none -screen none || exit 1
"$LMP" -in in.atomic -log none -screen none || exit 1
ls -l ./*.lammpstrj ./*.lammpstrj.gz ./*.data
