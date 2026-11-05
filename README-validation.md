synCrino-validation

This code was used to compare three different approaches regarding the synCrino python script. The outcome include the resulting ranks, RMSD and executing times.

Requirements: 

python version: 3.12,
pandas, numpy, matplotlib, scipy, tqdm;
One database in the standardized database excel sheet, one experimental 13C NMR data of a sample in the standardized sample excel sheet


how to call:
df_dp, df_bf, df_3l, top10 = run_and_compare_benchmarks(df_database, sample_peaks, exact_cc=False, filename =)


Raw NMR Datasets: Raw spectral data files (Bruker FIDs) used for algorithm validation are permanently archived on Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX.

Repository DOI: This repository version is archived on Zenodo under DOI https://doi.org/10.5281/zenodo.YYYYYYY.
This project is licensed under the MIT License – see the LICENSE file for details.

If you use this code or database in your research, please cite our paper:
