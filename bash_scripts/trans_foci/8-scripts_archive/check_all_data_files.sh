for parentdir in N*-all_simu*/; do
    cd $parentdir
        for dir in N*[1-8]/;do
            fname=$(echo ${dir} | cut -d / -f 1)
            echo $fname >> ../all_file_header.txt
            cd $dir
                head -n 25 ${fname}.all.data >> ../../all_file_header.txt
                echo >> ../../all_file_header.txt
            cd ..
        done
    cd ..
done