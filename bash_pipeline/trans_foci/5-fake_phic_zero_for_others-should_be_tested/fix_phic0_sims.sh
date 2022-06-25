#!/bin/bash
# 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
read -rp "column to replace > " col
read -rp "row to replace > " row

rename -v dir -s "ac$sigOld" "ac$sigNew" eps*
rename -s "ac$sigOld" "ac$sigNew" eps*/eps*.csv
rename -s "ac$sigOld" "ac$sigNew" eps*/eps*.npy
find eps*/ -type f -name "eps*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
#find eps*/ -type f -name "eps*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
find eps*/ -type f -name "eps*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +