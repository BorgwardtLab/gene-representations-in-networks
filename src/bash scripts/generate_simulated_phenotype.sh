#!/usr/bin/env bash
ROOT="/links/groups/borgwardt/Projects/master_project_julia"

DATA_ROOT="${ROOT}/data/athaliana/atwell_data"

OUT_DIR="${DATA_ROOT}/athaliana/atwell_data/FRI"
mkdir -p "${OUT_DIR}"

# Directory containing the original geneotype files
ORIGINAL_DATA="${DATA_ROOT}/genotype"

# The root to the data created with the simulated phenotype
SIMULATED_DATA_ROOT="${ROOT}/data/simulation_data"

EXEC_PYTHON="${ROOT}/code"

PLINK_PATH="/links/groups/borgwardt/agkbshare/software/LINUX64/tools/plink2/plink"

# Find 5 genes that are connected in the network

python $EXEC_PYTHON/find_gene_subset.py --data_filepath ${DATA_ROOT} --network_filepath ${ROOT}/data/athaliana/networks

# Merge the snp mappings
python $EXEC_PYTHON/merge_snp_mappings.py

# Extract the snp subset for the 5 genes
python $EXEC_PYTHON/extract_snp_subset_from_5genes.py --data_filepath ${DATA_ROOT} --allow-no-sex

#Extract the overalapping SNP genotypes
$PLINK_PATH --file $ORIGINAL_DATA  --extract $DATA_ROOT/snp_set.txt --make-bed --out $DATA_ROOT/simulation_data --allow-no-sex

# Recombine simulation data
$PLINK_PATH --bfile $OUT_DIR/simulation_data  --recode12 --out $DATA_ROOT/simulation_data_recombined

# Generate a simulated phenotype of weak association strength
python $EXEC_PYTHON/generate_simulated_phenotype.py --data_filepath $DATA_ROOT

