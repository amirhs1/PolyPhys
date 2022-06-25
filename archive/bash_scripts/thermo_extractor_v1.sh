#!/bin/bash
grep -m 1 "Step Temp" > thermo.txt
grep '^[[:blank:]]*[^[:blank:]A-z-]'  N80epsilon5.0r3.0lz35.5sig0.3nc36693ens1.txt | grep -v % | tail -n +15 >> thermo.txt
