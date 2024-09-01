#!/bin/bash
# SumRuleCyl:
#for dd in N*[1-8]-probe/; do
# HnsCub, HnsCyl:
#for dd in N*ring-probe/; do
# TransFociCub, TransFociCyl, SumRuleCubHeteroRing, SumRuleCubHeteroLinear:
for dd in ns*-probe/; do
    cd $dd
    for dir in al*.linear/; do # SumRuleCubHeteroLinear
    #for dir in al*.ring/; do # TransFociCub, SumRuleCubHeteroRing
    #for dir in eps*.ring/; do # TransFociCyl
    #for dir in N*[1-8].ring/; do # HnsCub HnsCyl
    #for dir in N*[1-8]/; do # SumRuleCyl
        echo "$dir"
        cd "$dir" || exit
        name=${dir:0:-1}
        datafile="$name".all.data
        cp "${datafile}" "${name}".all.original.data
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
    cd ..
done
echo "Finished!"
