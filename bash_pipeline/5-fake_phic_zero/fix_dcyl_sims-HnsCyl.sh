#!/bin/bash
# 
dcylOld=9.0
dcylNew=8.0
# The following part change the value of the dcrowd or D column
col=9

# mac os x: instaal rename and gawk
#find N*/ -type f -name "N*stamps*.csv" | rename -s "D$dcylOld" "D$dcylNew"
# linux
#rename "D$sigOld" "D$sigNew" N*
#rename "D$sigOld" "D$sigNew" N*/N*.csv
#rename "D$sigOld" "D$sigNew" N*/N*.npy
#rename "D$sigOld" "D$sigNew" N*/N*.txt

#find N*/ -type f -name "N*stamps.csv" -exec sed -i "" "s/D$dcylOld/D$dcylNew/g" {} \; # osx
#find N*/ -type f -name "N*stamps-ensAvg.csv" -exec sed -i "" "s/D$dcylOld/D$dcylNew/g" {} \; # osx
#find N*/ -type f -name "N*stamps.csv" -exec sed -i "s/D$dcylOld/D$dcylNew/g" {} + # Linux
#find N*/ -type f -name "N*stamps-ensAvg.csv" -exec sed -i "s/D$dcylOld/D$dcylNew/g" {} + # Linux
# The following part change the value of the D column


#for row in {2..41}; do
#    find N*/ -type f -name "N*stamps.csv" -exec gawk -i inplace -v new="${dcylNew}" -v col="7" -v row="$row" 'BEGIN{FS=OFS=","} NR == row  { $col = new } 1' {} +
#done

#find N*/ -type f -name "N*stamps.csv" -exec gawk -v new="8" -v col="9" 'BEGIN{FS=OFS=","} NR > 1  { col = new }' {} +

for file in N*/N*stamps.csv; do
    gawk -i inplace -v new="${dcylNew}" -v col="9" -F',' -v OFS=',' 'NR>=2 && NR<=41 && NR%1==0 { $col = new } 1' "$file"
done
#for d in N*-nucleoid-ensAvg/; do
#    for row in {2..11}; do
#        find ${d} -type f -name "N*stamps-ensAvg.csv" -exec gawk -i inplace -v new="${dcylNew}" -v col="7" -v row="$row" 'BEGIN{FS=OFS=","} NR % row == 0 { $col = new } 1' {} +
#    done
#done