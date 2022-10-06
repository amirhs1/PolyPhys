#!/bin/bash
for dir in ns*.ring/; do
    echo "$dir"
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
    head -n 10 "${dir}${name}.all.data" | tail -n 7 >> "${dir}${name}".bug.data 
    tail -n +11 fake_nc0_all.data >> "${dir}${name}".bug.data
done

