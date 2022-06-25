for parentdir in N*-all_simu*/; do
    cd $parentdir
        for dir in N*[1-8]/;do
            fname=$(echo ${dir} | cut -d / -f 1)
            echo $fname
            cd $dir
                gzip -dk all*data.gz
                mv all*.data ${fname}.all.data
            cd ..
        done
    cd ..
done