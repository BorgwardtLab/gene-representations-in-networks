
ROOT="/links/groups/borgwardt/Projects/master_project_julia"

EXEC_PYTHON="${ROOT}/code"

DATA_ROOT="${ROOT}/data/athaliana/atwell_data"

OUT_DIR="${ROOT}/data/athaliana/processed_FRI/FRI"
mkdir -p "${OUT_DIR}"

# Directory containing the original geneotype files
ORIGINAL_DATA="${DATA_ROOT}/genotype"

# File containing the filtered phenotypes
FILTERED_PHENO="${ROOT}/data/athaliana/processed_FRI/FRI/FRI.pheno"

# Output file for the filtered binary genotype files
FILTERED_GENO="${OUT_DIR}/filtered_FRI"

# Output file for the filtered binary genotype files
BINARY_GENO="${OUT_DIR}/binary_FRI"

PLINK_PATH="/links/groups/borgwardt/agkbshare/software/LINUX64/tools/plink2/plink"

PVALUE_RESULTS="${OUT_DIR}/FRI-linear.assoc.linear"

SKATO_RESULTS="${OUT_DIR}/FRI-skato"

SKATO_PATH="/links/groups/borgwardt/Projects/master_project_julia/code_AG/appl/skato_genes.sh"

# Convert the data into binary files
#$PLINK_PATH --file $ORIGINAL_DATA --make-bed --out $BINARY_GENO --allow-no-sex

# Filter the data
#$PLINK_PATH --bfile $BINARY_GENO --keep $FILTERED_PHENO --make-bed --out $FILTERED_GENO --allow-no-sex

# Include covariates in the data 
#$PLINK_PATH --bfile $FILTERED_GENO --out ${OUT_DIR}/binary_FRI --pca  header

# Perform linear association study
#$PLINK_PATH --bfile $FILTERED_GENO --out $OUT_DIR/FRI-linear --allow-no-sex --linear --assoc --pheno $FILTERED_PHENO

# Map each gene to the corresponding SNP:s
#python $EXEC_PYTHON/map_snps_to_genes.py --filepath $ROOT/data

# From the gene-snp mapping and the linear assoc studies, map gene to p-values
python $EXEC_PYTHON/find_gene_pval_representation.py --data_filepath $DATA_ROOT --pvalue_results $PVALUE_RESULTS 

# If this is run on a simulated phenotype add the phenotype to the data in order to run SKATO 
#$PLINK_PATH --bfile $BINARY_GENO  --make-bed --out ${OUT_DIR}/FRI_skato --pheno $FILTERED_PHENO --allow-no-sex

# Create covariate file for the SKAT-O
#$PLINK_PATH --bfile ${OUT_DIR}/FRI_skato --out ${OUT_DIR}/FRI_skato --pca  header

# Run the SKAT-O script to obtain the p-value representations for SKAT-O
$SKATO_PATH

# Create the set file
python $EXEC_PYTHON/create_set_file.py --data_filepath $DATA_ROOT 

# Compute the set-tests representation using PLINK
$PLINK_PATH --bfile $FILTERED_GENO --set-test --set $DATA_ROOT/gene_snp_set_file.txt --mperm 1000 --set-p 1 --set-r2 0.0 --set-max 1000 --assoc --linear --out $OUT_DIR/set_FRI --pheno $FILTERED_PHENO --allow-no-sex



