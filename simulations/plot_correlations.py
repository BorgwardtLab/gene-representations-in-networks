import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
from matplotlib import rc

def read_in_json_file(filepath):
    df_snp_mapping = pd.DataFrame(columns=['gene_id', 'nr_snp'])
    with open(filepath) as json_file:
        data = json.load(json_file)
        for p in data:
            df_snp_mapping = df_snp_mapping.append({'gene_id': p, 'nr_snp': data[p]}, ignore_index=True)
    return df_snp_mapping

def plot_gene_length_correlation(x, y, method, correlation):
    #fig = plt.figure()
    plt.rcParams.update({'font.size': 15})
    plt.scatter(x, y, alpha=0.7, s=0.9)
    plt.title('$R^2$: %.3f' %correlation)
    plt.ylabel('Number of SNPs')
    plt.xlabel(method + ' -log10($p$-value)')
    plt.show()
    #fig.savefig('../plots/new_correlation_plots/nsnp' + name + '.png', dpi=fig.dpi)


def plot_method_correlation(x, y, method1, method2, correlation, xval, yval):
    #fig = plt.figure()
    plt.rcParams.update({'font.size': 15})
    plt.scatter(x, y, alpha=0.7, s=0.8, label="")
    plt.plot(range(xval), linestyle='dotted', color='r', label='bisecting line')

    plt.title('$R^2$: %.3f' %correlation)
    plt.xlabel(method1 + ' -log10($p$-value)')
    plt.ylabel(method2 + ' -log10($p$-value)')
    plt.ylim(0, yval)
    plt.legend()
    plt.show()
    #fig.savefig('../plots/new_correlation_plots/'+ s1 + s2 + simulated + '.png',  dpi=fig.dpi)

def calculate_correlation(x_values, y_values):
    correlation_matrix = np.corrcoef(x_values, y_values)
    correlation_xy = correlation_matrix[0, 1]
    r_squared = correlation_xy ** 2
    return r_squared


def plot_all_gene_length_correlations(df_gene_count, df_pvalues):
    merged_df = df_gene_count.merge(df_pvalues, on='gene_id')
    x_data = merged_df['min_p_val']
    y_data = merged_df['nr_snp']
    r2 = merged_df.corr()**2

    plot_gene_length_correlation(x_data, y_data, r'\texttt{Minimum}', r2['nr_snp']['min_p_val'])

    x_data = merged_df['av_p_val']
    y_data = merged_df['nr_snp']


    plot_gene_length_correlation(x_data, y_data, r'\texttt{Fisher}', r2['nr_snp']['av_p_val'])

    x_data = merged_df['P.value']
    y_data = merged_df['nr_snp']

    plot_gene_length_correlation(x_data, y_data, r'\texttt{SKAT-O}', r2['nr_snp']['P.value'])

    x_data = merged_df['p-value']
    y_data = merged_df['nr_snp']
    r2 = merged_df.corr() ** 2

    plot_gene_length_correlation(x_data, y_data, r'\texttt{Set-based}', r2['nr_snp']['p-value'])


def transform_pvalues(df):
    df_pvalue_log10_results = df.copy()

    df_pvalue_log10_results['av_p_val'] = -np.log10(df['av_p_val'])
    df_pvalue_log10_results['min_p_val'] = -np.log10(df['min_p_val'])
    df_pvalue_log10_results['SKATO_p_val'] = -np.log10(df['SKATO_p_val'])
    return df_pvalue_log10_results


def read_and_plot_for_chromosome(current_chr):
    gene_count_filepath = './snp_count_per_gene_' + str(
        current_chr) + '.txt'
    df_gene_count = read_in_json_file(gene_count_filepath)

    pvalues_filepath = '/Users/juliaortheden/PycharmProjects/master_project/MLCB_project/gene_pvalue_representations/p-val_representation_chr' + str(
        current_chr) + '.txt'
    df_pvalue_results = pd.read_csv(pvalues_filepath)
    df_pvalue_results_log10 = transform_pvalues(df_pvalue_results)
    plot_all_gene_length_correlations(df_gene_count, df_pvalue_results_log10)

def plot_QQ_vs_PCs(x_covar, y_lambda, data):
    fig = plt.figure()
    plt.style.use('seaborn-darkgrid')
    plt.rcParams.update({'font.size': 15})
    plt.plot(x_covar, y_lambda, '.', alpha=0.8, markersize=3)
    plt.xlabel('number of leading PCs used as covariates')
    plt.xticks(np.arange(0, 21, step=1))
    plt.ylabel('Genomic inflation factor \lambda')
    plt.axes().set_yscale('log')
    plt.show()
    fig.savefig('../plots/QQ' + data + 'new.png', dpi=fig.dpi)




if __name__ == '__main__':
    plt.style.use('seaborn-darkgrid')
    rc('text', usetex=True)

    # Set filepaths
    gene_count_filepath = './nr_snp_total.txt'
    pvalue_filepath = './p-val_permutation_total.txt'

    # Read in the data into dataframes
    df_gene_count = pd.read_csv(gene_count_filepath)
    df_pvalue_results = pd.read_csv(pvalue_filepath)
    df_SKATO_pvalue_results = pd.read_csv("./permutation_data/permuted_SKATO_weak_results.txt",
                                          delim_whitespace=True, usecols=["SetID", "P.value"])
    df_set_pvalues = pd.read_csv('./my_hierarchical_hotnet_data/set_weak/data/scores_1.tsv',
                                 delim_whitespace=True, names=["gene_id", "p-value"])

    # Rename column in order to merge dataframes
    df_SKATO_pvalue_results = df_SKATO_pvalue_results.rename(columns={"SetID": "gene_id"})
    merged_pvalues = pd.merge(pd.merge(df_SKATO_pvalue_results, df_pvalue_results, on=["gene_id"]), df_set_pvalues, on=['gene_id'])
    df_pvalue_log10_results = transform_pvalues(merged_pvalues)
    names = [r'\texttt{Minimum}', r'\texttt{Fisher}', r'\texttt{SKAT-O}', r'\texttt{Set-based}']
    methods = ['min_p_val', 'av_p_val', 'P.value', 'p-value']

    # Plot gene length correlation for all methods
    plot_all_gene_length_correlations(df_gene_count, df_pvalue_log10_results)

    # Plot method correlation between all pairs of methods

    # Min vs Set
    x_data = df_pvalue_log10_results[methods[0]]
    y_data = df_pvalue_log10_results[methods[3]]
    r = df_pvalue_log10_results.corr()**2
    r_min_av = r['p-value']['min_p_val']

    plot_method_correlation(x_data, y_data, names[0], names[3], r_min_av, 30, 15)

    # Plot fisher vs min
    x_data = df_pvalue_log10_results[methods[0]]
    y_data = df_pvalue_log10_results[methods[1]]
    r = df_pvalue_log10_results.corr() ** 2
    r_min_av = r['av_p_val']['min_p_val']

    plot_method_correlation(x_data, y_data, names[0], names[1], r_min_av, 30, 30)

    # Min vs SKATO
    x_data = df_pvalue_log10_results[methods[0]]
    y_data = df_pvalue_log10_results[methods[2]]
    r_min_skat = r['P.value']['min_p_val']

    plot_method_correlation(x_data, y_data, names[0], names[2], r_min_skat, 30, 10)

    # Fisher vs SKATO
    x_data = df_pvalue_log10_results[methods[1]]
    y_data = df_pvalue_log10_results[methods[2]]
    r_min_skat = r['P.value']['av_p_val']

    plot_method_correlation(x_data, y_data, names[1], names[2], r_min_skat, 20, 15)

    # Set vs SKATO
    x_data = df_pvalue_log10_results[methods[2]]
    y_data = df_pvalue_log10_results[methods[3]]
    r_min_skat = r['P.value']['p-value']

    plot_method_correlation(x_data, y_data, names[2], names[3], r_min_skat, 15, 15)

    # Fisher vs Set
    x_data = df_pvalue_log10_results[methods[1]]
    y_data = df_pvalue_log10_results[methods[3]]
    r_min_skat = r['p-value']['av_p_val']
    plot_method_correlation(x_data, y_data, names[1], names[3], r_min_skat, 30, 15)


