for dir in N*/;do
	echo $dir
	cd $dir
	for file in  *-ens_avg*.csv;do
		header=$(echo $file | cut -d - -f 1)
		footer=$(echo $file | cut -d - -f 3)
		property=$(echo $footer | cut -d . -f 1)
		final=${header}-${property}-ens_avg.csv
		mv $file $final
	done
	cd ..
done
