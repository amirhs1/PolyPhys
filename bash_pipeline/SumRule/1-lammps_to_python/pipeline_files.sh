#!/bin/bash
# copies the files needed for different stage of the pipeline from this directory to a 'space' directory - a 'space' directory is a directory in which all the topology and trajectory files of a 'space' group are located.
# bug data file generating step
read -p "Enter bug size > " nbug
read -p "Enter cylinder diameter > " dcyl
template=data_template-cylinder_bug_n${nbug}d${dcyl}.data
mkdir PipeLine
cp ../sumrule_pipeline/"${template}" .
cp ../sumrule_pipeline/datafile_generator.py .
cp ../sumrule_pipeline/PipeLine/*.py ./PipeLine
cp ../sumrule_pipeline/3-data_file_fixing/fix_PairIJ_all_data_files_v*.sh .

# extraction bug data using gnuparallel
cp ../sumrule_pipeline/4a-extraction-first_time/extraction*.sh .
cp ../sumrule_pipeline/4a-extraction-first_time/sbatch*.sh .
cp ../sumrule_pipeline/4a-extraction-first_time/submit*.sh .
cp ../sumrule_pipeline/extractionbug_no_crowd_*.py .
cp ../sumrule_pipeline/extraction_all_in_one_*.py .
cp ../sumrule_pipeline/PipeLine/*.py .

# add phic zero to extraction file (it is necessary for some sets of input parameters)
cp ../sumrule_pipeline/5-fake_phic_zero_for_others/check_properties.sh . 
cp ../sumrule_pipeline/5-fake_phic_zero_for_others/fix_properties_files_new_without_comment.sh . 
cp ../sumrule_pipeline/5-fake_phic_zero_for_others/phic_zero_rename_directories.sh .
cp ../sumrule_pipeline/5-fake_phic_zero_for_others/phic_zero_rename_files.sh . 

# analyze bug and create clean csv files, log file and create clean csv files
cp ../sumrule_pipeline/6-analysis/trjs_bug_directory_gen.sh .
cp ../sumrule_pipeline/6-analysis/analyze_bug_directory_gen.sh . 
cp ../sumrule_pipeline/analyze_*.py .
