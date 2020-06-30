import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rc

def convert_pvalues_to_log10(df):
    df_pvalue_log10_results = df.copy()
    df_pvalue_log10_results['P'] = -np.log10(df['P'])
    return df_pvalue_log10_results

if __name__ == '__main__':
    pval_filepath = '/Users/juliaortheden/PycharmProjects/master_project/MLCB_project/permutation_data/filtered_linear.assoc.linear'
    snp_locations_filepath = '/Users/juliaortheden/ResearchProject/data/atwell_data/genotype.map'
    df_snp_subset = pd.read_csv('/Users/juliaortheden/PycharmProjects/master_project/MLCB_project/permutation_data/gene_snp_subset_5genes.txt',
                                delim_whitespace=True, names=["gene", "SNP"])
    df_snp_locations = pd.read_csv(snp_locations_filepath, delim_whitespace=True, names=["CHR", "SNP", "START", "END"])
    df_transformed_pvalues = convert_pvalues_to_log10(pd.read_csv(pval_filepath, delim_whitespace=True))

    # Merge the dataframes containing all the log 10 transformed p-values with the snp locations
    df_pvalue_snp_locations = pd.merge(df_transformed_pvalues, df_snp_locations, on=['SNP'])

    # Merge the dataframes containing only the log 10 transformed p-values with the snp locations for the induced SNPs
    df_induced_snps_locations = pd.merge(pd.merge(df_transformed_pvalues, df_snp_subset, on=['SNP']),
                                         df_snp_locations, on=['SNP'])

    plt.style.use('seaborn-darkgrid')
    fig = plt.figure()
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 15})
    offset = 0

    for chromosome in range(1, 6):
        # Plot all the SNP p-values and locations for the current chromosome
        df_current_chromosome = df_pvalue_snp_locations.loc[df_pvalue_snp_locations['CHR_x'] == chromosome]
        chr_data = df_current_chromosome['END'].astype(int)
        chr_data += offset
        plt.scatter(chr_data, df_current_chromosome['P'], alpha=0.7, s=0.8, label="Chromosome " + str(chromosome))

        # Mark the induced SNPs for the current chromosome
        snps = df_induced_snps_locations.loc[df_induced_snps_locations['CHR_x'] == chromosome]
        chr_data = snps['END'].astype(int)
        chr_data += offset
        # Plot the label for the induced SNPs only once,for the last chromosome
        if chromosome < 5:
            plt.plot(chr_data, snps['P'], 'x', c='r', label="_nolegend_")
        else:
            plt.plot(chr_data, snps['P'], 'x', c='r', label="Induced 36 SNPs")
        offset = offset + max(df_current_chromosome['END'].astype(int))


    # Set parameters for the plots, i.e. size of legend and control axes ticks
    legend = plt.legend(loc="best", markerscale=5., scatterpoints=1, fontsize=13)
    legend.legendHandles[0]._legmarker.set_markersize(6)
    plt.title('Manhattan plot of strong association')
    plt.xlabel('SNP position')
    plt.ylabel('-log10($p$-value)')
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False
        )
    plt.show()
    # Set the path where to save the plot
    fig.savefig('../plots/manhattan_strong.png', dpi=fig.dpi)