#!/bin/bash
for dir in eps*.ring/; do
    echo "$dir"
    cd "$dir" || exit
    name=${dir:0:-1}
    echo "$name"
    cp ../fake_nc0_all.data "${name}".bug.data
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
    # The following command does not work on Max OSX
    #sed -i "" -e "4,10R ${alldata}" -e '4,10d' "${name}".bug.data
    # linux sed: This does not work
    #sed -i -e "4,10R ${name}.all.data" -e '4,10d' "${name}".bug.data
    cd ..
done

