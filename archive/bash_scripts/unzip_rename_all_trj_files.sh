for dir in N*[1-8]/;do
    fname=$(echo ${dir} | cut -d / -f 1)
    echo $fname
    cd $dir
        for gzfile in all.*.lammpstrj.gz;do 
            gzip -dk $gzfile
            allfile=$( echo $gzfile | cut -d . -f -4)
            j=$(echo $allfile | cut -d . -f 3)
            mv $allfile ${fname}.${j}.all.lammpstrj
        done
    cd ..
done