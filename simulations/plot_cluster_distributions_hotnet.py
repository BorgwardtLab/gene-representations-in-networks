import numpy as np
from itertools import groupby
import matplotlib.pyplot as plt
from matplotlib import rc

# The root to the HotNet data locally
hierarchical_hotnet_filepath = "/Users/juliaortheden/PycharmProjects/master_project/MLCB_project/generated_data/my_hierarchical_hotnet_data"
# The file containing the indices mapped to the genes of the network (same for all methods)
df_index_genes = hierarchical_hotnet_filepath + '/medium_own_permutations/data_permutation_fisher/network_1_index_gene.tsv'
association_strength = 'medium'

# File paths to the clustering outputs for the different methods
fisher_fn = hierarchical_hotnet_filepath + '/hotnet_' + association_strength + '/results_fisher/clusters_network_1_scores_1.tsv'
min_fn = hierarchical_hotnet_filepath + '/hotnet_' + association_strength + '/results_min/clusters_network_1_scores_1.tsv'
skato_fn = hierarchical_hotnet_filepath + '/hotnet_' + association_strength + '/results_SKATO/clusters_network_1_scores_1.tsv'
set_fn = hierarchical_hotnet_filepath + '/set_'+ association_strength + '/results/clusters_network_score_1.tsv'

def read_hotnet(fn):
    clusters = []
    with open(fn, 'r') as fin:
        for line in fin:
            if line.startswith('#'): continue
            tmp = line.strip().split()
            clusters.append(tmp)
    return clusters

def plot_clusters(clusters, title, name):
    clusters_sizes = np.sort([len(x) for x in clusters])
    x = []
    y = []

    # Group by cluster sizes and count the occurences in order to plot the cluster distributions
    for key, group in groupby(clusters_sizes):
        x.append(key)
        y.append(len(list(group)))

    # Set plot parameters
    plt.style.use('seaborn-darkgrid')
    fig = plt.figure()
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 15})
    plt.bar(x, y)
    plt.title(title)
    plt.ylabel('Occurences')
    plt.xlabel('Cluster size')
    plt.yscale('log')
    plt.show()
    fig.savefig('../plots/new_correlation_plots/histogram_' + name + '_' + association_strength + '.png')


if __name__ == '__main__':

    # Read the HotNet clusters for the respective method
    skato_clusters = read_hotnet(skato_fn)
    fisher_clusters = read_hotnet(fisher_fn)
    min_clusters = read_hotnet(min_fn)
    set_clusters = read_hotnet(set_fn)

    # Plot the cluster size distribution of the respective method
    plot_clusters(skato_clusters, r'\texttt{SKAT-O}', 'skato')
    plot_clusters(fisher_clusters, r'\texttt{Fisher}', 'fisher')
    plot_clusters(min_clusters, r'\texttt{Minimum}', 'min')
    plot_clusters(set_clusters, r'\texttt{Set-based}', 'set')