#!/bin/bash
# This bash scripts compress all the files in all_simulations directories and
# them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move node on your cluster.
# For graham, it is gra-dtn1.sharcnet.ca you need to submit this scrip with
# nohup ... &  command.
name=$(pwd | rev | cut -d / -f 1 | rev) #report name
report=${name}-archive_report.txt
touch "$report"
for directory in N*/; do

    echo "${directory}" >> "${report}"
    dir=${directory::-1}
    zipdir=$dir"-zip"
    mkdir "${zipdir}"
    echo "$dir"
    gzip -vrk "${dir}" >> "${report}"
    cd "$directory" || exit
        if [ -d restarts  ]; then
            cd restarts || exit
                mkdir restarts_zip
                mv ./*.gz restarts_zip
                mv restarts_zip ../../"${zipdir}"/
            cd ..
        fi
        mv ./*.gz ../"${zipdir}"/
    cd ..
    project=${dir}.tar
    tar -zcvkf "${project}" "${zipdir}" >> "${report}"
done
runname=$(echo $name | cut -d '-' -f 1)
rundir=run_files-${runname}
tar -zcvkf ${rundir}.tar.gz ${rundir}
echo "Finished!"