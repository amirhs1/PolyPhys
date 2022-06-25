#!/bin/bash
# 
read -rp "Enter old size of crowders (regex friendly) > " sig_old
read -rp "Enter new size of crowders (regex friendly) > " sig_new

for dir in N*nc0*/; do
    echo "$dir"
    cd "$dir" || exit
    echo ""
    for file in N*stamps.csv; do
        echo "$file"
        fsdsedname=$(head -n 1 "$file" | tr -d '[:space:]')
        #sed -i "s/filename.*/${fsdsedname}/g" ${file} # linux
        sed -i '' "s/filename.*/${fsdsedname}/g" "${file}"
        dataline=$(tail -n 1 "$file")
        newline=$(echo "$dataline" | awk -v new="$sig_new" 'BEGIN{FS=OFS=","} {$11=new ; print}')
        #sed -i "s/N.*/${newline}/g" ${file} # linux
        sed -i "" "s/N.*/${newline}/g" "${file}"
        sed -i "" "s/.(sig$sig_old})./.(sig${sig_new})./g" "${file}"
        tail -n 1 "$file" 
        echo ""
    done
    echo ""
    cd ..
done