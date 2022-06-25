read -p "Enter old size of crowders > " sig_old
read -p "Enter new size of crowders > " sig_new

for dir in N*nc0*/;do
    echo $dir
    csvfile=$(echo $dir | cut -d / -f 1)
    echo ${csvfile}
    cd $dir
    file=${csvfile}_properties.csv
    echo $file
    fsdsedname=$(head -n 2 $file | tail -n 1| tr -d '[:space:]')
    sed -i "s/filename.*/${fsdsedname}/g" ${file}
    dataline=$(tail -n 1 $file)
    newline=$(echo $dataline | awk -v new="$sig_new" 'BEGIN{FS=OFS=","} {$10=new ; print}')
    sed -i "s/N.*/${newline}/g" ${file}
    filename=$( tail -n 1 $file | cut -d , -f 1)
    newfilename=$(echo "${filename/sig${sig_old}/sig${sig_new}}")
    newline=$( tail -n 1 $file  | awk -v new="$newfilename" 'BEGIN{FS=OFS=","} {$1=new ; print}')
    sed -i "s/N.*/${newline}/g" ${file}
    tail -n 3 $file
    cd ..
done