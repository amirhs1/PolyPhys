#!/bin/bash

#SBATCH --ntasks=1 						  	# number of mpi tasks (cpus); for a serial job, it is the default; 1.
#SBATCH --mem-per-cpu=3000M      						# memory; default unit is megabytes.
#SBATCH --time=01-00:00                	# Total time needed for this job. --time (DD-HH:MM) or -t 0-01:00;  Check the manual of sbatch
#SBATCH --account=def-byha 					# The user/account holder who use th computecanada project. you can use -a instead of --acount
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com 	# The mail address for notifications about submitted job
#SBATCH --mail-type=ALL 					# Set the types of notifications which will be emailed.
#SBATCH --job-name="archive_job"

echo "Starting run at: $(date)"
name=$(pwd | rev | cut -d / -f 1 | rev) #report name
report=${name}-archive_report.txt
touch "$report"
for directory in N*[1-8]/; do
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
for directory in N*[1-8]_res/; do
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
for directory in N*[1-8]_incomplete/; do
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
for directory in N*[1-8]_cont/; do
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
echo "Program finished with exit code $? at: $(date)"
