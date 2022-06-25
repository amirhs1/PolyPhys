for dir in N*_incomplete/; do
    sim=$( echo "$dir" | cut -d _ -f 1)
    echo "$sim"
    mkdir rerun/$sim
    cd "$dir" || exit
    cp input.lmp ../rerun/"$sim"/
    cp submit.sh ../rerun/"$sim"/
    cp "${sim}".restart ../rerun/"$sim"/
    cd ..
done
