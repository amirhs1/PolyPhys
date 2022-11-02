currentname=$(pwd | rev | cut -d / -f 1 | rev) #report name
name=$( echo $currentname | cut -d - -f 1)
extractionall=${name}-extraction_all
mkdir ${extractionall}
for dir in N*ep*/; do
    mkdir ${extractionall}/$dir
    echo $dir
    cd $dir
    mv N*.csv ../${extractionall}/$dir
    rm *.py
    rm -r PipeLine
    cd ..
done