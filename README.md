This repository contains code and data from **Fluoroquinolone resistance-conferring _gyrA_ variants alter the fitness cost and potentiate the resistance of the zoliflodacin resistance mutation _gyrB_<sup>D429N</sup> in _Neisseria gonorrhoeae_**\
Authors: Aditi Mukherjee, Sofia OP Blomqvist, David Helekal, Daniel HF Rubin, Bailey Bowcutt, Samantha G Palace, Yonatan H Grad\
Preprint: [https://www.biorxiv.org/content/10.64898/2026.05.04.722797v1](url)

## **Pairwise Data**
- Raw data_pairwise fitness_growth curve.xlsx: Excel workbook of raw data from pairwise fitness assays and growth curve experiments

## **Barcode Data**
- gyra_gyrb_barcode_counts.csv: Raw results of isogenic strain competition experiment. Read counts of each barcode over 3 timepoints (0, 2, and 4 hours) across 2 replicates (A/B).
- gyra+gyrb_effects.csv: median and 95% CI of model draws from bc_regression_model.R for each genotype. For genotypes with *gyrB* D429N, median and CI are calculated from the sum of the corresponding *gyrA* draws and *gyrB* given *gyrA* genotype draws.

## **Scripts**
- bc_count.py: python file for counting barcode sequences from amplicon sequencing output
- bc_regression_model.R: regression model for calculation of _gyrA_ and _gyrB_ fitness effects from barcode count data

## **Extras**

- bc_count.yml: can be used to create a conda environment for running bc_count.py

- convert_to_model_format.py: can be used to convert from prop_to_ref files from multiple replicates (bc_count.py outputs) to single input formatted for bc_regression_model.R. This script requires an extra text file containing the barcode names with associated genotypes (columns for strain, *gyrA* genotype, and 0/1 for *gyrB*)
  
- sample_reffile: example of correct format for barcode references fasta file for bc_count.py
  
- sample_genotypefile: example of genotype file for convert_to_model_format.py
