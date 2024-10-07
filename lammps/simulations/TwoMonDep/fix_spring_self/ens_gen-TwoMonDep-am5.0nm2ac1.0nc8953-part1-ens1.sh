#!/bin/bash
# Generated on 20240930
# Parameters: sig1, n1, sig2, n2, hl, s_d, rand_seed
P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.0	10000"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.1	10010"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.2	10020"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.3	10030"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.4	10040"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.5	10050"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.6	10060"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.7	10070"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.8	10080"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	0.9	10090"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.0	10100"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.1	10110"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.2	10120"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.3	10130"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.4	10140"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.5	10150"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.6	10160"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.7	10170"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.8	10180"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	1.9	10190"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.0	10200"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.1	10210"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.2	10220"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.3	10230"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.4	10240"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.5	10250"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.6	10260"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.7	10270"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.8	10280"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	2.9	10290"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	3.0	10300"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	3.5	10310"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	4.0	10320"
#P[${#P[@]}]="5.0	2	1.0	8953	12.5	4.5	10330"
P[${#P[@]}]="5.0	2	1.0	8953	12.5	5.0	10340"

dt=0.0005
adump=5000
bdump=2000
tdump=2000

# Initialize ensemble values
ens='1' # Global variable defining total ensembles for simulation.

# Loop through ensemble values
for i in ${ens}; do 
    # Loop through parameter sets
    for ((j=0; j < ${#P[@]}; j++ )); do
        # Extract parameters from array P
        read -r sig1 n1 sig2 n2 hl sd randseed_group <<< "${P[$j]}"

        # Define directory and simulation names
        dirname="am${sig1}nm${n1}ac${sig2}nc${n2}hl${hl}sd${sd}dt${dt}bdump${bdump}adump${adump}tdump${adump}ens${i}"
        simname="am${sig1}nm${n1}ac${sig2}nc${n2}hl${hl}sd${sd}dt${dt}bdump${bdump}adump${adump}tdump${adump}ens${i}"

        # Create simulation input file from template
        cp TwoMonDep_tmpl_input.lmp input.lmp
        {
            echo "# ---------- Input Parameters ----------"
            echo "variable        i equal ${i} # Simulation counter"
            echo "variable        sig1 equal ${sig1} # Diameter of monomers"
            echo "variable        n1 equal ${n1} # Number of monomers"
            echo "variable        sig2 equal ${sig2} # Diameter of crowders"
            echo "variable        n2 equal ${n2} # Number of crowders"
            echo "variable        sd equal ${sd} # Surface-to-surface distance between monomers"
            echo "variable        hl equal ${hl} # Box size half-length"
            echo "variable        randseed equal $((i+randseed_group)) # Random seed"
            echo "variable        simname string ${simname} # Simulation/run name"
        } | cat - input.lmp > temp && mv temp input.lmp

        # Create directory if it does not exist and move input file
        mkdir -p "${dirname}"
        mv input.lmp "${dirname}"
        
        # Select appropriate submission script based on total number of particles
        total_particles=$((n2 + n1))
        if [ "${total_particles}" -ne "${n1}" ]; then
            if [ "${total_particles}" -le 5000 ]; then
                cp submit_cpu4_le5000.sh submit.sh
            elif [ "${total_particles}" -le 15000 ]; then
                cp submit_cpu4_gt5000_le15000.sh submit.sh
            elif [ "${total_particles}" -le 30000 ]; then
                cp submit_cpu8_gt15000_le30000.sh submit.sh
            elif [ "${total_particles}" -le 75000 ]; then
                cp submit_cpu16_gt30000_le75000.sh submit.sh
            elif [ "${total_particles}" -le 150000 ]; then
                cp submit_cpu32_gt75000_le150000.sh submit.sh
            elif [ "${total_particles}" -le 200000 ]; then
                cp submit_cpu32_gt150000_le200000.sh submit.sh
            elif [ "${total_particles}" -le 250000 ]; then
                cp submit_cpu32_gt200000_le250000.sh submit.sh
            else
                cp submit_cpu32_gt250000_le300000.sh submit.sh
            fi
        else
            cp submit_nocrowd_8hrs.sh submit.sh
        fi
        
        # Move the submission script to the directory
        mv submit.sh "${dirname}"
        
        # Enter the directory and create necessary sub-directory
        cd "${dirname}" || exit
        mkdir -p restarts
        cd ..
    done
done
