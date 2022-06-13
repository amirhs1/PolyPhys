for dir in N*-zip/;do
    newname=$(echo $dir | cut -d - -f 1)
    echo $newname
    mv $dir $newname
done