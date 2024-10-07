#!/bin/bash

P=(
    "1000 24 0 2 0 8 628.5 100000"
    "1000 24 30 2 0 8 628.5 700000"
    "1000 24 30 2 4243 8 628.5 700010"
    "1000 24 30 2 21212 8 628.5 700090"
    "1000 24 30 1 33939 8 628.5 700100"
    "1000 24 30 1 169695 8 628.5 700180"
)

ens='1 2'
data_dir="equilibrated_data_files-N1000kbmm24r8nh0ac2lz628.5nc0ens1.j03.ring.dnaMinimize"
first_step=250200000
delta_step=200000

counter=0

for i in $ens; do
    for param in "${P[@]}"; do
        read -r N kBend n_hns sig4 n_crowd r lz randseed_group <<<"$param"

        dir_name="N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}.ring"
        sim_name="N${N}kbmm${kBend}r${r}nh${n_hns}ac${sig4}lz${lz}nc${n_crowd}ens${i}"

        if [[ ! -d "$dir_name" ]]; then
            mkdir "$dir_name"
            file_number=$((first_step + counter * delta_step))
            file_name="time_step-$file_number.data"
            cp "./$data_dir/$file_name" "$dir_name/minimized.dna.data"

            cp hns.mol "$dir_name/"

            if [[ $n_crowd -eq 0 ]]; then
                cp hns-cylindrical_v1-ac_equal_a.lmp input.lmp
            else
                if [[ $sig4 -ne 1 ]]; then
                    cp hns-cylindrical_v1.lmp input.lmp
                else
                    cp hns-cylindrical_v1-ac_equal_a.lmp input.lmp
                fi
            fi

            sed -e "1i\\
				variable simname string ${sim_name}" \
					-e "1i\\
				variable randseed equal $((i + randseed_group))" \
					-e "1i\\
				variable i equal $i" \
					-e "1i\\
				variable r equal $r" \
					-e "1i\\
				variable lz equal $lz" \
					-e "1i\\
				variable n_crowd equal $n_crowd" \
					-e "1i\\
				variable sig4 equal $sig4" \
					-e "1i\\
				variable n_hns equal $n_hns" \
					-e "1i\\
				variable kBend11 equal $kBend" \
					-e '1i\\
				# Defining input parameters:' input.lmp > "$dir_name/input.lmp"

            N_no_crowd=$((N + n_hns))

            if [[ $((n_crowd + N_no_crowd)) -ne $N_no_crowd ]]; then
                if [[ $((n_crowd + N)) -le 5000 ]]; then
                    cp submit_cpu4_le5000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 5000 && $((n_crowd + N)) -le 15000 ]]; then
                    cp submit_cpu4_gt5000_le15000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 15000 && $((n_crowd + N)) -le 30000 ]]; then
                    cp submit_cpu8_gt15000_le30000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 30000 && $((n_crowd + N)) -le 75000 ]]; then
                    cp submit_cpu16_gt30000_le75000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 75000 && $((n_crowd + N)) -le 150000 ]]; then
                    cp submit_cpu32_gt75000_le150000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 150000 && $((n_crowd + N)) -le 200000 ]]; then
                    cp submit_cpu32_gt150000_le200000.sh "$dir_name/submit.sh"
                elif [[ $((n_crowd + N)) -gt 200000 && $((n_crowd + N)) -le 250000 ]]; then
                    cp submit_cpu32_gt200000_le250000.sh "$dir_name/submit.sh"
                else
                    cp submit_cpu32_gt250000_le300000.sh "$dir_name/submit.sh"
                fi
            else
                cp submit_nocrowd_20hrs.sh "$dir_name/submit.sh"
            fi

            mkdir "$dir_name/restarts"
        fi

        counter=$((counter + 1))
    done
done
