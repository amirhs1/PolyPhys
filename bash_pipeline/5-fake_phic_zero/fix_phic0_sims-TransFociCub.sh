#!/bin/bash
# Double check the mass 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
read -rp "Enter old mass of crowders > " massOld
read -rp "Enter new mass of crowders > " massNew

# The following part gives the col and row of the value of the "dcrowd" value.
read -rp "column to replace crowder size > " col # 10
read -rp "row to replace crowder size > " row # 2
# The following part gives the col and row of the value of the "mcrowd" value.
read -rp "column to replace crowder mass > " colMass # 19
read -rp "row to replace crowder mass > " rowMass # 2

# mac os x
#rename -v dir -s "sig$sigOld" "sig$sigNew" al*
#rename -s "sig$sigOld" "sig$sigNew" al*/al*.csv
#rename -s "sig$sigOld" "sig$sigNew" al*/al*.npy
#rename -s "sig$sigOld" "sig$sigNew" al*/al*.txt
# linux
rename "ac$sigOld" "ac$sigNew" al*
rename "ac$sigOld" "ac$sigNew" al*/al*.csv
rename "ac$sigOld" "ac$sigNew" al*/al*.npy
rename "ac$sigOld" "ac$sigNew" al*/al*.txt

#find al*/ -type f -name "al*stamps.csv" -exec sed -i "" "s/sig$sigOld/sig$sigNew/g" {} +
find al*/ -type f -name "al*stamps.csv" -exec sed -i "s/sig$sigOld/sig$sigNew/g" {} + # Linux
#find al*/ -type f -name "al*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
find al*/ -type f -name "al*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux

# The following part change the value of the "dcrowd" value.
find al*/ -type f -name "al*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +
# The following part change the value of the "mcrowd" value.
find al*/ -type f -name "al*stamps.csv" -exec gawk -i inplace -v new="${massNew}" -v col="$colMass" -v row="$rowMass" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +