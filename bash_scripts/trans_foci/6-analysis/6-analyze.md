# How to use analyze data for bug 
 1. Request a computing node with **salloc --ntasks=1 --cpus-per-task=4 --mem=8G --time=1:0:0 --account=def-someuser**.
 2. Active **dasnEnv**.
 3. Run the **python analyze_bug.py** in the **____-extraction** directory.
 4. Change the name of **nohup.out** to  **N___D___ac___phic___-analyze_bug-report.out**.
 5. Correct the name of **...-log_details.csv**, **...-log_runtime.csv**, **...-properties.csv**, and **...-properties-ens_avg.csv**, based on the **phic range**.
 6. Run the *python analyze_distributions_segments.py* in the *____-extraction_bug* directory.
 7. Change the name of **nohup.out** to  **N___D___ac___phic___-analyze_all-report.out**.
 8. Run **analyze_dirs_gen.sh** to organize all the files.
 9. Go to **N___D___ac___phic___-analyze_bug** directory and cp **....log...** report to the **cylinder_logs_csvs** folder.
 