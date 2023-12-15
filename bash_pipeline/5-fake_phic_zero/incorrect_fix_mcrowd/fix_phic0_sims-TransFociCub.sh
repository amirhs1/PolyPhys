#!/bin/bash
# 
read -rp "Enter old size of crowders > " sigOld
read -rp "Enter new size of crowders > " sigNew
# The following part change the value of the dcrowd or ac column
read -rp "column to replace > " col
read -rp "row to replace > " row

# mac os x
rename -v dir -s "ac$sigOld" "ac$sigNew" al*
rename -s "ac$sigOld" "ac$sigNew" al*/al*.csv
rename -s "ac$sigOld" "ac$sigNew" al*/al*.npy
rename -s "ac$sigOld" "ac$sigNew" al*/al*.txt
# linux
#rename "ac$sigOld" "ac$sigNew" al*
#rename "ac$sigOld" "ac$sigNew" al*/al*.csv
#rename "ac$sigOld" "ac$sigNew" al*/al*.npy
#rename "ac$sigOld" "ac$sigNew" al*/al*.txt

find al*/ -type f -name "al*stamps.csv" -exec sed -i "" "s/ac$sigOld/ac$sigNew/g" {} +
#find al*/ -type f -name "al*stamps.csv" -exec sed -i "s/ac$sigOld/ac$sigNew/g" {} + # Linux
# The following part change the value of the dcrowd or ac column
find al*/ -type f -name "al*stamps.csv" -exec gawk -i inplace -v new="${sigNew}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +