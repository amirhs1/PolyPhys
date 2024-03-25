#!/bin/bash
# Double check the mass 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
read -rp "Enter old mass of crowders > " massOld
read -rp "Enter new mass of crowders > " massNew
# The following part gives the col and row of the value of the "dcrowd" value.
read -rp "column to replace > " col # 10
read -rp "row to replace > " row # 2
# The following part gives the col and row of the value of the "mcrowd" value.
read -rp "column to replace mass > " colMass # 19
read -rp "row to replace mass > " rowMass # 2
# mac os x
#rename -v dir -s "sig$sigOld" "sig$sigNew" N*
#rename -s "sig$sigOld" "sig$sigNew" N*/N*.csv
#rename -s "sig$sigOld" "sig$sigNew" N*/N*.npy
#rename -s "sig$sigOld" "sig$sigNew" N*/N*.txt
# linux
rename "sig$sigOld" "sig$sigNew" N*
rename "sig$sigOld" "sig$sigNew" N*/N*.csv
rename "sig$sigOld" "sig$sigNew" N*/N*.npy
rename "sig$sigOld" "sig$sigNew" N*/N*.txt

#find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/sig$sigOld/sig$sigNew/g" {} +
find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/sig$sigOld/sig$sigNew/g" {} + # Linux
#find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
# The following part change the value of the "dcrowd" value.
find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +
# The following part change the value of the "mcrowd" value.
find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${massNew}" -v col="$colMass" -v row="$rowMass" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +