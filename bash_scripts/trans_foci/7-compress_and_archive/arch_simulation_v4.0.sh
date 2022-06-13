#!/bin/bash
# This bash scripts compress all the files in all simulations directories and them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move NODe on your cluter. For graham, it is gra-dtn1.sharcnet.ca
# you need to submit this scrip with nohup ... &  command.
report=$(pwd | rev | cut -d / -f 1 | rev)_archive_report.txt #report name
touch $report
for directory in N*/; do

    echo ${directory} >> ${report}
    dir=${directory::-1}
    zipdir=$dir"-zip" 
    mkdir ${zipdir}

    echo $dir

    gzip -vrk ${directory} >> ${report}

    cd $directory
        mv restarts restarts_old

        cd restarts_old
            mkdir restarts_zip
            mv *.gz restarts_zip
            mv restarts_zip ..
        cd ..

        mv restarts_zip ../${zipdir}/
        mv *.gz ../${zipdir}/
    cd ..

    project=$dir".tar"
    tar -cvf ${project} ${zipdir} >> ${report}
done