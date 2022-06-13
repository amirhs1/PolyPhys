for dir in N*/;do
    echo $dir
    csvfile=$(echo $dir | cut -d / -f 1)
    echo ${csvfile}
    cd $dir
    file=${csvfile}_properties.csv
    echo $file
    fsdsedname=$(head -n 2 $file | tail -n 1| tr -d '[:space:]')
    sed -i "s/filename.*/${fsdsedname}/g" ${file}
    tail -n 3 $file
    cd ..
done