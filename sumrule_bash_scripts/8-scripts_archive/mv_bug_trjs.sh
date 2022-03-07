for dir in N*/; do
    echo $dir
    mv ${dir::-1}*bug* $dir
done