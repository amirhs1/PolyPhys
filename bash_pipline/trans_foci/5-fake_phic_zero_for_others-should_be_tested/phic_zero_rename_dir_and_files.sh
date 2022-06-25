#!/bin/bash
# change the value of crowder size in all the files.
read -rp "Enter old size of crowders > " sig_old
read -rp "Enter new size of crowders > " sig_new
for dir in eps*/;do
    cd "$dir" || exit
    for file in eps*.npy;do
        mv -i -- "$file" "${file//sig${sig_old}/sig${sig_new}}"
    done
    for file in eps*.csv;do
        mv -i -- "$file" "${file//sig${sig_old}/sig${sig_new}}"
    done
    cd ..
    mv -i -- "$dir" "${dir//sig${sig_old}/sig${sig_new}}"
done