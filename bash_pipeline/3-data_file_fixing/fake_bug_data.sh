#!/bin/bash

# Use this only for the TransFoci project.
# Becareful; this script is correct only if we want to create a "bug" data file
# from tje "bug" data file of a simulation with nc=0. In these two simulations,
# everything is similar EXCEPT the number of crowders.
read -rp "Enter the line number of the last line of the header section in 'all' data file  > " allLine

# For TransFociCub project, allLine=31
# For TransFociCyl project, allLine=31
# For SumRuleCyl project, allLine=26
# For HnsCub project, allLine=38
headLine=$((allLine - 3))
nextLine=$((allLine + 1))

#for dir in eps*.ring/; do # TransFociCyl
#for dir in N*.ring/; do # SumRuleCyl
#for dir in N*.ring/; do # HnsCub
for dir in al*.ring/; do # TransFociCub
    echo "${dir}"
    name=${dir:0:-1}
    # The sed command copies the first 10 lines of all.data to bug.data
    # files to ensure that the box dimensions are the correct in the
    # generated bug.data files. This is done becuase all the input
    # parameters for a given simulation in a space do not change, except
    # nc and lz or l.
    # nc is the number of crowders and is the only physically meaningful
    # parameter that is changed within a space. lz (in the cylindrical
    # geometry) or l (in the free space or cubic geometry) are the lengths
    # simulation box in a period direction and are set to make nc as small
    # as possible to have a give volume fraction of crowders in the system
    head -n 3 fake_nc0_all.data > "${dir}${name}".bug.data
    head -n "${allLine}" "${dir}${name}.all.data" | tail -n "${headLine}" >> "${dir}${name}.bug.data"
    tail -n "+${nextLine}" fake_nc0_all.data >> "${dir}${name}.bug.data"
done
