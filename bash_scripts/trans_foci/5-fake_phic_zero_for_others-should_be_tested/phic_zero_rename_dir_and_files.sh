#!/bin/bash
# change the value of crowder size in all the files.
read -p "Enter old size of crowders > " sig_old
read -p "Enter new size of crowders > " sig_new
for dir in N*/;do
    cd "$dir" || exit
    for file in N*.npy;do
        mv -i -- "$file" "${file//sig${sig_old}/sig${sig_new}}"
    done
    for file in N*.csv;do
        mv -i -- "$file" "${file//sig${sig_old}/sig${sig_new}}"
    done
    cd ..
    mv -i -- "$dir" "${dir//sig${sig_old}/sig${sig_new}}"
done