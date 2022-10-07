#!/bin/bash
# This bash scripts compress all the files in all_simulations directories and
# them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move node on your cluter.
# For graham, it is gra-dtn1.sharcnet.ca you need to submit this scrip with
# nohup ... &  command.
report=$(pwd | rev | cut -d / -f 1 | rev)-archive_report.txt #report name
touch "$report"
for directory in eps*/; do

    echo "${directory}" >> "${report}"
    dir=${directory::-1}
    zipdir=$dir"-zip" 
    mkdir "${zipdir}"
    echo "$dir"
    gzip -vrk "${dir}" >> "${report}"
    cd "$directory" || exit
        cd restarts || exit
            mkdir restarts_zip
            mv ./*.gz restarts_zip
            mv restarts_zip ..
        cd ..

        mv restarts_zip ../"${zipdir}"/
        mv ./*.gz ../"${zipdir}"/
    cd ..

    project=${dir}.tar
    tar -cvkf "${project}" "${zipdir}" >> "${report}"
done
