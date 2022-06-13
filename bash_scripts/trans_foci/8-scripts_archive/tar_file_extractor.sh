# This is an example on how to extract a file from a tar file.
for tfile in N*.tar; do
    echo ${tfile}
    file=$(echo ${tfile} | sed -r 's/(.*)\..*/\1/')
    tar -xvf ${tfile} ${file}/log.lammps.gz
done
