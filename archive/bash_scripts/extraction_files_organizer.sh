for dir in N2000*/;do
	echo $dir
	cd $dir
	for ens in N*/; do
		cd ${ens}
		file=$(echo ${ens} | cut -d "/" -f 1)
		#mkdir ${file}
		echo ${file}
		mv ${file}_ensemble.csv ${file}_properties.csv
		mv ${file}_rfloryEdges.csv ${file}_rFloryEdges.csv
		#for f in *.txt; do
		#	mv -- "$f" "${f%.txt}.csv"
		#done
		#mv ${file}*.csv ${file}/
		#mv ${file}*.pdf ${file}/
		cd ..
	done
	cd ..
done
