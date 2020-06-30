import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from compute_precision_recall import read_hotnet

def convert_pvalues_to_log10(df):
    df_pvalue_log10_results = df.copy()
    df_pvalue_log10_results['P'] = -np.log10(df['P'])
    return df_pvalue_log10_results

if __name__ == '__main__':
    plt.style.use('seaborn-darkgrid')
    plt.rcParams.update({'font.size': 15})

    # Read in the results from the HotNet clustering
    results_fn = './my_hierarchical_hotnet_data/strong_own_permutations/results_simulated_permutation_data_min_medium/clusters_network_1_scores_1.tsv'
    # Extract the genes found in the largest cluster
    hotnet_genes = read_hotnet(results_fn)[0]
    # Read in the genes in the orginial input score file for HotNet
    original_score_filepath = './my_hierarchical_hotnet_data/strong_own_permutations/data_permutation_min/score_1.tsv'
    df_original_scores = pd.read_csv(original_score_filepath, delim_whitespace=True, names=["gene", "p-value"])

    # Plot the genes HotNet found in the largest cluster instead of the 5 induced genes
    df_genes_in_largest_cluster = pd.DataFrame(hotnet_genes, columns=["gene"])
    df_genes_in_largest_cluster_scores = pd.merge(df_genes_in_largest_cluster, df_original_scores, on=['gene'])

    for i in range(1, 500):
        permutation_pval_filepath = './my_hierarchical_hotnet_data/strong_own_permutations/permutation_min/score_' + str(
            i) + '.txt'
        df_permutation_scores = pd.read_csv(permutation_pval_filepath, delim_whitespace=True, names=["gene", "p-value"])
        df_permutation_scores_for_largest_cluster = pd.merge(df_genes_in_largest_cluster, df_permutation_scores, on=['gene'])
        plt.scatter(df_permutation_scores_for_largest_cluster['gene'], df_permutation_scores_for_largest_cluster['p-value'], alpha=0.7, s=0.8)

    plt.title('Distribution of p-values for strong phenotype using Minimum')
    plt.xlabel('Gene')
    plt.ylabel('-log10(p-value)')
    plt.plot(df_genes_in_largest_cluster_scores['gene'], df_genes_in_largest_cluster_scores['p-value'], 'x', c='r')
    plt.show()