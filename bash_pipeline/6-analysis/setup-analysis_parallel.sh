#!/bin/bash

# Check if a project name is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <PROJECT_NAME>"
    exit 1
fi

project=$1

# copies the python script and package which are then executed by gnuparallal.
case $project in
    HnsCub|HnsCyl|SumRuleCylWhole|SumRuleCylSegment)
        dir_pattern="N*-probe*/"
        ;;
    TransFociCyl)
        dir_pattern="D*-probe*/"
        ;;
    TransFociCub|SumRuleCubHeteroRing|SumRuleCubHeteroLinear)
        dir_pattern="ns*-probe*/"
        ;;
    *)
        echo "Invalid project name: $project"
        exit 1
        ;;
esac

for dir in $dir_pattern; do
    echo "$dir"
    cp analysis-${project}.py analysis_phase.py 
    mv analysis_phase.py ./"$dir"
    cp -R polyphys ./"${dir}"
done

