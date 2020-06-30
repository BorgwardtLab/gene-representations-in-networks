#!/usr/bin/env bash

ROOT="/links/groups/borgwardt/Projects/master_project_julia"

# The output directory for the assoc studies
OUT_DIR="${ROOT}/data/permutation_data/permutation_data/permutation_gwas/medium_min"
mkdir $OUT_DIR

# The path to run PLINK
PLINK_PATH="/links/groups/borgwardt/agkbshare/software/LINUX64/tools/plink2/plink"

# The path to the file containing the PCs which to include in the model
COVAR="${ROOT}/data/permutation_data/permutation_data/plink/binary-genotype-covar.eigenvec"

# The path to the files containing the binary data files (filtered with PLINK's --keep on the phenotype, see run_plink.sh)
BINARY_FILES="${ROOT}/data/permutation_data/permutation_data/plink/binary-filtered_permuted"

# The path to the directory containing the phenotype permutations (i.e pheno1.txt, pheno2.txt ...)
PHENOTYPE="${ROOT}/data/permutation_data/permutation_data/phenotype_permutations_medium_min"

# Run a GWAS for each of the 500 permutations
for n_perm in {1..500}
do

$PLINK_PATH --bfile $BINARY_FILES --out $OUT_DIR/linear_$n_perm --linear \
 --allow-no-sex --adjust --covar-number 1-10 --covar $COVAR --pheno $PHENOTYPE/pheno$n_perm.txt

done
