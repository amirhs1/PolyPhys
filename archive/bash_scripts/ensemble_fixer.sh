for dir in N*/; do
	file=$(echo ${dir} | cut -d / -f 1)
    cd $dir
	mv ${file}_ensemble.txt ${file}_properties.csv
	for f in *.txt; do 
    mv -- "$f" "${f%.txt}.csv"
	done
	cd ..
done
