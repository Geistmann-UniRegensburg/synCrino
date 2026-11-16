# synCrino-validation

This code was used to compare three different approaches regarding the synCrino python script. The outcome include the resulting ranks, RMSD and executing times.
This repository contains the validation and benchmarking codebase supporting the methodological evaluation of automated 13C NMR spectrum alignment algorithms.To demonstrate computational efficiency, mathematical robustness, and accuracy, three distinct algorithmic matching strategies are benchmarked against one another across full reference datasets.
The three matching strategies are: 3-For-Loops (hierarchical 3-step grid search), Brute Force (BF) and Monotonic Dynamic Programming (DP)

Requirements: 

python version: 3.12

pandas, numpy, matplotlib, scipy, tqdm

Input files: One database in the standardized database excel sheet, one experimental 13C NMR data of a sample in the standardized sample excel sheet


# How to Use:

from synCrino_validation import run_and_compare_benchmarks, load_database, load_sample

#Load database and test sample:

database = load_database("path/to/standardized_database.xlsx")

sample = load_sample("path/to/standardized_sample.xlsx")


#Run full comparative benchmark

df_dp, df_bf, df_3l, top10 = run_and_compare_benchmarks(
    database=database, 
    sample=sample, 
    exact_cc=False, 
    filename="NMR_Validation_Benchmark_Results.xlsx")





Raw NMR Datasets: Raw spectral data files (Bruker FIDs) used for algorithm validation are permanently archived on Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX.

Repository DOI: This repository version is archived on Zenodo under DOI https://doi.org/10.5281/zenodo.YYYYYYY.
This project is licensed under the MIT License – see the LICENSE file for details.

If you use this code or database in your research, please cite our paper:
