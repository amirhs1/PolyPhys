for dir in N*-all_simulations/;do 
    echo $dir
    temp=$( echo $dir | cut -d . -f 2)
    ac=${temp:3:4}
    cd $dir
    for childir in N*/;do
        cd $childir
        echo $childir
        file=$( echo $childir | cut -d / -f 1)
        datafile=${file}.all.data
        sed -i "/PairIJ Coeffs.*/i Pair Coeffs # lj_cut\n\n1 1 1\n2 1 ${ac}\n" $datafile 
        sed -i "/PairIJ Coeffs.*/,+5d" $datafile 
        cd ..
    done
    cd ..
done