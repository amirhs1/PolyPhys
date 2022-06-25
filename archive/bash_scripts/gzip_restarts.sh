for dir in N*-restarts/;do
    cd $dir
    for file in N*.restart.gz; do
        gzip -d $file
        done
    cd ..
done