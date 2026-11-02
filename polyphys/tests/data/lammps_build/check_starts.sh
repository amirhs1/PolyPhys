#!/bin/sh
# Build the start files with make_starts.py and run each under the three
# integrating modes of in.kg; print the warnings and errors of each log,
# apart from the communication-cutoff heuristic that in.kg describes.
# The same matrix runs as the lammps-marked tests of test_build.py.
#
# Run from this directory, in the project environment, with the executable
# in LMP (default lmp_serial):
#   LMP=lmp_serial sh check_starts.sh [outdir]
LMP=${LMP:-lmp_serial}
OUT=${1:-./starts}
python make_starts.py "$OUT" || exit 1
status=0
for case in "helix bond 0" "ring bond 0" "ring_angle angle 0" \
            "helix bond 5.2"; do
  set -- $case
  name=$1; style=$2; dcyl=$3
  for run in "nve 0.005" "nvelang 0.01" "limitlang 0.005"; do
    set -- $run
    mode=$1; dt=$2
    log="$OUT/$name-$mode-dcyl$dcyl.log"
    "$LMP" -in in.kg -log "$log" -screen none -var data "$OUT/$name.data" \
      -var style "$style" -var mode "$mode" -var dt "$dt" -var dcyl "$dcyl"
    warnings=$(grep '^WARNING' "$log" | grep -vc 'Communication cutoff')
    errors=$(grep -c '^ERROR' "$log")
    finished=$(grep -c '^Loop time' "$log")
    printf '%-12s %-10s dcyl=%-4s warnings=%s errors=%s finished=%s\n' \
      "$name" "$mode" "$dcyl" "$warnings" "$errors" "$finished"
    if [ "$warnings" != 0 ] || [ "$errors" != 0 ] || [ "$finished" != 1 ]
    then
      status=1
    fi
  done
done
exit $status
