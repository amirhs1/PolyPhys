#!/bin/bash
# Double check the mass 
acs="2 3 4 5 6"
for d in ns*ac1*probe; do
    echo $d
    cd $d
    for sigNew in $acs; do
        echo $sigNew
        mkdir $sigNew
        cp -R al* ${sigNew}/
        cd $sigNew
        sigOld=1
        massOld=1
        massNew=$((${sigNew}*${sigNew}*${sigNew}))
        col=10
        row=2
        colMass=19
        rowMass=2
        rename "ac$sigOld" "ac$sigNew" al*
        rename "ac$sigOld" "ac$sigNew" al*/al*.csv
        rename "ac$sigOld" "ac$sigNew" al*/al*.npy
        rename "ac$sigOld" "ac$sigNew" al*/al*.txt
        find al*/ -type f -name "al*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
        find al*/ -type f -name "al*stamps.csv" -exec gawk -i inplace -v new="${sigNew}.0" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +
        find al*/ -type f -name "al*stamps.csv" -exec gawk -i inplace -v new="${massNew}.0" -v col="$colMass" -v row="$rowMass" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +
        cd ..
    done
    cd ..
done
