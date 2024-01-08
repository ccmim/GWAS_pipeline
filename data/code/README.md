# Genetic data pre-processing
## Imputed genotypes
The following instructions allow to filter the set of imputed genotype data (in BGEN format) for a set of SNPs and subjects, and to split these files into regions to ease the parallelization of the subsequent tasks.
1. Generate a preliminary set of SNPs to filter for, so that subsequent tasks are faster, based on MAF and imputation INFO score.
2. Filter for SNPs generated in step 1, for subjects and split the genome into regions.
3. Compute SNP-wise statistics (MAF, missingness) and further filter the BGEN files.

### Step 1
Generate a file with the SNPs to be kept in a preliminary filtering stage:
```
bash generate_include_snps_list.sh
```

### Step 2
Filter BGENs as per the list of SNPs generated in 1, a list of subjects, and split into small regions:
```
python split_bgen_by_region_qctool.py
```

### Step 3
Compute SNP-wise statistics and further filter the BGENs from step 2:
```
python filter_bgens_for_qc_criteria.py
```
