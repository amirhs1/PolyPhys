for datafile in N*.all.data;do
    echo $datafile
    head -n 25 $datafile
done