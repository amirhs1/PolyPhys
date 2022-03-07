#!/bin/bash
# This bash scripts compress all the files in all simulations directories and them tar the directories in a *____-all_simulations* directory.
# README: to use this file: connect to the Data Move NODe on your cluter. For graham, it is gra-dtn1.sharcnet.ca
# you need to submit this scrip with nohup ... &  command.
for dir in N*-all_simulations/;do
    cp arch_simulation_v1.0.sh $dir
    cd $dir
    nohup bash arch_simulation_v1.0.sh & 
    cd ..
done
