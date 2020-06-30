#!/usr/bin/env bash

ROOT="/links/groups/borgwardt/Projects/master_project_julia"

# The root to the original atwell data
DATA_ROOT="${ROOT}/data/athaliana/atwell_data"

# The root to the data created with the simulated phenotype
SIMULATED_DATA_ROOT="${ROOT}/data/permutation_data/permutation_data"

# Directory containing the original geneotype files
ORIGINAL_DATA="${DATA_ROOT}/genotype"

# The output directory for the processed files
OUT_DIR="${SIMULATED_DATA_ROOT}/plink"
mkdir -p "${OUT_DIR}"

# File containing the phenotype
SIMULATED_PHENO="${SIMULATED_DATA_ROOT}/simulated_phenotype_medium.txt"

# Output file for the filtered binary genotype files
FILTERED_BINARY_DATA="${OUT_DIR}/filtered_binary_data"

# Binary version of the original files
BINARY_DATA="${OUT_DIR}/binary_data"

# The path to run PLINK
PLINK_PATH="/links/groups/borgwardt/agkbshare/software/LINUX64/tools/plink2/plink"

# The resulting file with the PCAs for the simulated phenotype
COVAR="${OUT_DIR}/binary_data_covar.eigenvec"

# Convert the data into binary files
$PLINK_PATH --file $ORIGINAL_DATA --make-bed --out $BINARY_DATA

# Filter the data on the simulated phenotype
$PLINK_PATH --bfile $BINARY_DATA --keep $SIMULATED_PHENO --make-bed --out $FILTERED_BINARY_DATA

# Calculate the covariates for the simulated phenotype
$PLINK_PATH --bfile $FILTERED_BINARY_DATA --out $OUT_DIR/binary-genotype-covar --pca  header

# Perform linear association study for the filtered data using 10 PCs
$PLINK_PATH --bfile $FILTERED_BINARY_DATA --out $OUT_DIR/linear_medium --linear \
 --allow-no-sex --adjust --covar-number 1-10 --covar $COVAR --pheno $SIMULATED_PHENO


