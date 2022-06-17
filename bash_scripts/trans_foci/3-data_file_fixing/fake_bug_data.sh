#!/bin/bash
# fix the 'all' data files with coreecting the number of crowders.
for alldata in eps*.all.data; do
    sim=$(echo "${alldata}" | rev | cut -d . -f 3- | rev ) # -f -6)
    echo "${sim}.bug.data"
    cp fake_nc0_all.data "${sim}".bug.data
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
    # mac osx sed:
    # sed -i "" -e "10r ${alldata}" -e '1,10d' "${sim}".bug.data
    # linux sed:
    sed -i -e "10r ${alldata}" -e '1,10d' "${sim}".bug.data
done 