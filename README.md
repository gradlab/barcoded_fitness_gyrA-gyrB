# Barcoded fitness assays for isogenic _gyrA_ and _gyrB_ strains 

## **bc_count.py**

### **Inputs (edit at the top of the script before running)**

- reffile: fasta file of each barcode name/sequence

- fastqdir: directory of fastq files of sequencing reads for a given replicate. Each fastq should be named after its timepoint and the name of the fastq directory will be used as the replicate name (ex. fastqs/A/2.fastq contains the reads from timepoint 2 in replicate A, and the fastqdir in this situation would be 'fastqs/A')

- refstrain: name of the strain used as a reference (used for the prop_to_ref.csv output). If the reference strain has multiple barcode sequences, their names in the reffile should all start the same (ex. if barcodes are AM370.1 and AM370.2, refstrain is 'AM370').

### **Outputs**

bc_count_results directory containing the following files (with rep replaced with the name of the fastq directory):

- rep_rawcounts.csv: number of instances of each barcode at each timepoint

- rep_proportions.csv: proportion of each barcode to total barcode count for each timepoint

- rep_normalized.csv: proportion of each barcode at each timepoint normalized to the barcode's own time 0 value

- rep_prop_toref.csv: proportion of each barcode to sum of reference strain barcode counts at that timepoint. This data is used for bc_regression_model.R

## **bc_regression_model.R**

### **Input:** csv file containing a row for each barcode in each replicate at each timepoint with columns for:

- barcode: name of barcode/strain (should match those in the reffile fasta file used by bc_count.py)

- replicate: name of replicate

- genotype: genotype (not a necessary column but nice to have)

- gyra: *gyrA* pos91.pos95 (ex. S.D)

- gyrb: 0 if *gyrB* 429D, 1 if *gyrB* D429N

- unique_id: identifier for each barcode in each replicate (ex. AM390.1_A for all timepoints of AM390.1 in replicate A)

- time: time in generations

- prop_to_ref: proportion of barcode to sum of reference strain barcode counts at that timepoint (output of bc_count.py)

### **Output:**

- gyra+gyrb_effects.csv: median and 95% CI of all model draws. Contains rows for effect of *gyrA* (median difference in growth rate of each *gyrA* genotype compared to reference strain) and sum of *gyrA* effect and *gyrB* effect given *gyrA* genotype (median difference in growth rate of each *gyrA* genotype + *gyrB* D429N compared to reference strain).

## **Data**
- gyra_gyrb_barcode_counts.csv: Raw results of isogenic strain competition experiment. Read counts of each barcode over 3 timepoints (0, 2, and 4 hours) across 2 replicates (A/B).
- gyra+gyrb_effects.csv: median and 95% CI of model draws from bc_regression_model.R for each genotype. For genotypes with *gyrB* D429N, median and CI are calculated from the sum of the corresponding *gyrA* draws and *gyrB* given *gyrA* genotype draws.


## **Extras**

- bc_count.yml: can be used to create a conda environment for running bc_count.py

- convert_to_model_format.py: can be used to convert from prop_to_ref files from multiple replicates (bc_count.py outputs) to single input formatted for bc_regression_model.R. This script requires an extra text file containing the barcode names with associated genotypes (columns for strain, *gyrA* genotype, and 0/1 for *gyrB*)
  
- sample_reffile: example of correct format for barcode references fasta file for bc_count.py
  
- sample_genotypefile: example of genotype file for convert_to_model_format.py
