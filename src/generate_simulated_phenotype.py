import argparse

import pandas as pd
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_filepath')
    args = parser.parse_args()
    filepath = args.data_filepath

    # The file path leading to the data file containing the recombined 1-2 genotypes for the phenotype
    df_X = pd.read_csv(filepath + '/simulation_data_recombined.ped', delim_whitespace=True, header=None)

    # Exclude the six first columns since they include IID, FID etc.
    X = df_X.iloc[:, 6:].to_numpy()
    X_id = df_X.iloc[:, 0:2]

    offset = 0
    association_strength = np.random.uniform(0.01, 0.05, (72, 1))
    residual_error = np.random.normal(0, 1, (1307, 1))

    # The linear regression formula used to generate the phenotype
    y = offset + np.dot(X, association_strength) + residual_error

    # Keep the ID of the phenotype from the original file and combine with the new phenotype
    phenotype = pd.DataFrame(np.hstack([X_id, y]))

    phenotype[0] = [np.int(pheno) for pheno in phenotype[0]]
    phenotype[1] = [np.int(pheno) for pheno in phenotype[1]]

    phenotype.to_csv(filepath + '/simulated_phenotype_weak.txt',
                     sep=" ", header=None, index=False)