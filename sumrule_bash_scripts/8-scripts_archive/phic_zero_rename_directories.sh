read -p "Enter old size of crowders > " sig_old
read -p "Enter new size of crowders > " sig_new
for dir in N*/;do
    mv -i -- "$dir" "${dir//sig${sig_old}/sig${sig_new}}"
done