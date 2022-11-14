#!/bin/bash
# This part is different between different projects
#for dir in al*.ring/; do # TransFociCub
#for dir in eps*.ring/; do # TransFociCyl
#for dir in N*[1-8]/; do # HnsCub
for dir in N*[1-8]/; do # SumRuleCyl
    echo "$dir"
    cd "$dir" || exit
    name=${dir:0:-1}
    datafile="$name".all.data
    # the sed command below (s)ubstitute the line starts with
    # "PairIJ Coeffs" with the "Pair Coeffs # lj_cut"
    # Here we assumed that the Pair Coeffs are for the "lj_cut"
    # style.
    # mac osx sed:
    #sed -i "" "s/^PairIJ Coeffs.*/Pair Coeffs # lj_cut/" "${datafile}"
    # linux sed: 
    sed -i "s/^PairIJ Coeffs.*/Pair Coeffs # lj_cut/" "${datafile}"
    echo "report:"
    # Check the Pair info line with new pattern is in the file
    grep -n -E "Pair Coeffs" "${datafile}"
    cd ..
done