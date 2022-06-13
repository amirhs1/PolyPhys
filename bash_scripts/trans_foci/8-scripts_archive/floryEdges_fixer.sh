for dir in N*/; do
	file=$(echo ${dir} | cut -d / -f 1)
    cd $dir
	mv ${file}_rfloryEdges.txt ${file}_rFloryEdges.txt
	cd ..
done
