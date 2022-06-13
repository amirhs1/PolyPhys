currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
temp=$( echo $name | cut -d . -f 2)
ac=${temp:3:4}
echo $ac
for datafile in N*.all.data;do
    echo $datafile
    sed -i "/PairIJ Coeffs.*/i Pair Coeffs # lj_cut\n\n1 1 1\n2 1 ${ac}\n" $datafile 
    sed -i "/PairIJ Coeffs.*/,+5d" $datafile 
done