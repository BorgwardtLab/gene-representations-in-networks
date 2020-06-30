
ROOT="/links/groups/borgwardt/Projects/master_project_julia"

# File path to the SKAT-O R wrapper
EXEC="${ROOT}/code_AG/methods/skato_wrapper.R"

# The directory containing the simulated data
DATA_ROOT="${ROOT}/data/simulated_data"

OUT_DIR="${ROOT}/output/athaliana/atwell_data/SKATO/simulated_data"
mkdir -p "${OUT_DIR}"

# File containing the sets (e.g. genes).
SET_FN="SET_FN=${ROOT}/data/athaliana/processed_SKATO/TAIR_geneset.txt"

# Takes the binary file, the covariate file, the gene set file, 
# the output file and the number of covariates to be used as input
Rscript ${EXEC} \
    --bfile "${DATA_ROOT}/permutation_data_SKATO" \
    --cfile "${DATA_ROOT}/permutation_data_SKATO.eigenvec" \
    --mfile "${SET_FN}" \
    --opref "${OUT_DIR}/permuted_SKATO_medium" \
    --ncov "10"
