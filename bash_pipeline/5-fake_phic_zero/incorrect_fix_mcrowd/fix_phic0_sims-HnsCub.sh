#!/bin/bash
# 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
# The following part change the value of the dcrowd or ac column
read -rp "column to replace > " col
read -rp "row to replace > " row

# mac os x
rename -v dir -s "ac$sigOld" "ac$sigNew" N*
rename -s "ac$sigOld" "ac$sigNew" N*/N*.csv
rename -s "ac$sigOld" "ac$sigNew" N*/N*.npy
rename -s "ac$sigOld" "ac$sigNew" N*/N*.txt
# linux
#rename "ac$sigOld" "ac$sigNew" N*
#rename "ac$sigOld" "ac$sigNew" N*/N*.csv
#rename "ac$sigOld" "ac$sigNew" N*/N*.npy
#rename "ac$sigOld" "ac$sigNew" N*/N*.txt

find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
#find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
# The following part change the value of the dcrowd or ac column
find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +