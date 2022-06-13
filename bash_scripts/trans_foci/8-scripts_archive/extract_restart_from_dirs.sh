for tar in N*/;do
    name=$( echo $tar | cut -d / -f 1)
    #restart_after.i1.gz
    cp ${tar}restart_after.* ./${name}.restart
done