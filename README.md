# synCrino
A Python-based framework designed for automated matching, shift-correction, and visualization of experimental 13C NMR chemical shifts against reference databases.The algorithm accounts for systematic chemical shift offsets (e.g., solvent effects) and handles datasets with unequal signal counts by modeling signal alignment as a minimum-cost bipartite matching problem. So far, it is exemplified with polycyclic polyprenylated acylphloroglucinols (natural compounds)

Requirements: 

python version: 3.12

pandas, numpy, matplotlib, scipy

Input files: One database in the standardized database excel sheet, one experimental 13C NMR data of a sample in the standardized sample excel sheet


# How to Use:

from synCrino import comparison_similar_compound, load_database, load_sample

#Load database and test sample:

database = load_database("path/to/standardized_database.xlsx")

sample = load_sample("path/to/standardized_sample.xlsx")

#Run comparison against database (by default with exact carbon count matching):

results_df = comparison_similar_compound(
    database=database, 
    sample=sample, 
    exact_carbon_count=True
)

# Example and further Use:

As an testing example, there is the sample excel file "sample_hyperandroneA_CDCl3.xlsx" and its associated database excel file "Database_1.xlsx" provided.

As described in the paper, Database 1 contains caged bicyclo[3.3.1]nonane Type A acylphloroglucinols, and Database 2 compiles uncaged bicyclo[3.3.1]nonane Type A acylphloroglucinols.

For further use, there are the standardized excel files for sample ("sample_empty.xlsx") and database ("database_empty.xlsx") provided.



Raw NMR Datasets: Raw spectral data files (Bruker FIDs) used for algorithm validation are permanently archived on Zenodo at https://doi.org/10.5281/zenodo.XXXXXXX.

Repository DOI: This repository version is archived on Zenodo under DOI https://doi.org/10.5281/zenodo.YYYYYYY.
This project is licensed under the MIT License – see the LICENSE file for details.

If you use this code or database in your research, please cite our paper:
