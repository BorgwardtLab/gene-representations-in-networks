import numpy as np
import pandas as pd


phenotype_df = pd.read_csv('./permutation_data/simulated_phenotype_strong.txt',
                           names=["FID", "ID", "score"], sep=" ")
# Extract the phenotype scores to be permuted
phenotype = phenotype_df["score"]

# Set the number of permutations
number_of_permutations = 500
# Set the start index
start_index = 2

if __name__ == '__main__':
    phenotype_df.to_csv('./permutation_data/phenotype_permutations_strong_min/pheno'+ str(1) + '.txt',
                          sep=" ", index=False, header=None)
    for i in range(start_index, start_index + number_of_permutations):
       # Create a random permutation of the phenotype score for 500 phenotype permutations
       permuted_pheno = np.random.permutation(phenotype)
       permuted_df = phenotype_df.copy()
       permuted_df["score"] = permuted_pheno
       permuted_df.to_csv('./permutation_data/phenotype_permutations_strong_min/pheno'+ str(i) + '.txt',
                          sep=" ", index=False, header=None)