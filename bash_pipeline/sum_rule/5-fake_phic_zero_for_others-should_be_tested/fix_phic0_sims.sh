#!/bin/bash
# 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
read -rp "column to replace > " col
read -rp "row to replace > " row

rename -v dir -s "sig$sigOld" "sig$sigNew" N*
rename -s "sig$sigOld" "sig$sigNew" N*/N*.csv
rename -s "sig$sigOld" "sig$sigNew" N*/N*.npy
find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/sig$sigOld/sig$sigNew/g" {} +
#find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/sig$sigOld/sig$sigNew/g" {} + # Linux
find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
#find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +