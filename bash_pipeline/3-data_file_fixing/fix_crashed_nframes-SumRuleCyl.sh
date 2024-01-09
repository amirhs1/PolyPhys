#!/bin/bash
# 
read -rp "new nframe > " nframe
read -rp "column to replace > " col # 25
read -rp "row to replace > " row 
find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${nframe}" -v col="$col" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +

