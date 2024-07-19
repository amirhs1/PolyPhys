#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3000M
#SBATCH --time=01-00:00
#SBATCH --account=def-byha
#SBATCH --mail-user=mr.a.h.saadeghi@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="archive_job"

echo "Starting run at: $(date)"
report_name="$(basename "$(pwd)")-archive_report.txt"
touch "$report_name"

process_directory() {
    local pattern=$1
    local dirs=( $pattern )
    # Check if directories exist
    if [ ${#dirs[@]} -eq 0 ]; then
        echo "No directories found for pattern: $pattern"
        return
    fi
    
    for directory in "${dirs[@]}"; do
        echo "${directory}" >> "${report_name}"
        dir_name=${directory%/} # Remove the trailing slash
        zip_dir="${dir_name}-zip"
        mkdir "${zip_dir}"
        echo "$dir_name"
        gzip -vrk "${dir_name}" >> "${report_name}"
        if [ -d "${directory}restarts" ]; then
            mkdir "${zip_dir}/restarts_zip"
            gzip -vrk "${directory}restarts"/*.gz >> "${report_name}"
            mv "${directory}restarts"/*.gz "${zip_dir}/restarts_zip"
        fi
        mv "${directory}"*.gz "${zip_dir}/"
        tar_project="${dir_name}.tar"
        tar -zcvf "${tar_project}" "${zip_dir}" >> "${report_name}"
    done
}

# Directory patterns to process
# SumRuleCyl
#directory_patterns=("N*ens[1-8]/" "N*ens[1-8]_res/" "N*ens[1-8]_incomplete/" "N*ens[1-8]_cont/")
# TransFociCyl
#directory_patterns=("eps*ens[1-8]*ring/" "N*ens[1-8]*ring_res/" "N*ens[1-8]*ring_incomplete/" "N*ens[1-8]*ring_cont/")
# TransFociCub
directory_patterns=("al*ens[1-8]*ring/" "al*ens[1-8]*ring_res/" "al*ens[1-8]*ring_incomplete/" "al*ens[1-8]*ring_cont/")
# SumRuleCubHetero
#directory_patterns=("al*ens[1-8]*ring_res/" "al*ens[1-8]*ring_res/" "al*ens[1-8]*ring_incomplete/" "N*ens[1-8]_cont/")

for pattern in "${directory_patterns[@]}"; do
    # Expanding the pattern to check for existing directories
    dirs=( $pattern )
    if [ ${#dirs[@]} -gt 0 ]; then
        process_directory "$pattern"
    else
        echo "No directories found matching pattern: $pattern"
    fi
done

run_name="$(cut -d '-' -f 1 <<<"$report_name")"
run_dir="run_files-${run_name}"
tar -zcvf "${run_dir}.tar.gz" "${run_dir}"

echo "Finished!"
echo "Program finished with exit code $? at: $(date)"

