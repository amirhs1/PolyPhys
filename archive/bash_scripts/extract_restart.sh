for dir in N*/; do
    fname=$(echo ${dir} | cut -d / -f 1)
    echo $fname
    if test -f ${dir}restart_after* ; then
    cp ${dir}restart_after* ${fname}.restart
    fi
done