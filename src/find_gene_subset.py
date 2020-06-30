import argparse
import pandas as pd
import random



def filter_network(genes):
    df_new = pd.DataFrame(columns=['gene0', 'gene1'])
    for gene in df_gene_edges.iterrows():
        gene0 = gene[1].values[0]
        gene1 = gene[1].values[1]

        if gene0 in genes and gene1 in genes:
            df_new = df_new.append({'gene0': gene0, 'gene1': gene1}, ignore_index=True)
    return df_new


def find_genes_with_correct_snp_mappings():
    genes = []
    for gene in df_snp_genes.iterrows():
	#print(gene)
        number_of_snp = gene[1].values[1]
        gene_id = gene[1].values[0]
	print(number_of_snp)
	print(gene_id)
        # Add all genes that have 5-10 overlapping SNPs to a list
        if 4 < number_of_snp and number_of_snp < 11:
            genes.append(gene_id)
    return genes


def get_neighbourhood(gene_set, df_new):
    df_neighbours = pd.DataFrame()
    # Find neighbors to the genes in the gene set
    for gene in gene_set:
        df_neighbours = df_neighbours.append(df_new.loc[df_new['gene1'] == gene])
    return df_neighbours

def random_draw(neighbourhood):
    #random_index = random.randint(0, len(neighbourhood))
    return neighbourhood['gene0'][0]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--data_filepath')
    parser.add_argument('--network_filepath')
    args = parser.parse_args()
    data_filepath = args.data_filepath
    network_filepath = args.network_filepath

    df_gene_edges = pd.read_csv(str(network_filepath) + "/TairProteinInteraction.20090527_edges_sym.txt", sep="\t", names=['gene1', 'gene2'])
    df_snp_genes = pd.read_csv(str(data_filepath) + "/number_of_snps_per_gene.txt", delim_whitespace=True)
    
    gene_set= []
    # Find genes that have 5-10 overlapping SNPs
    genes = find_genes_with_correct_snp_mappings()
    # Filter so only genes that are in the TAIR network are included
    df_new = filter_network(genes)
    # Initialize the search for 5 genes by grabbing the first gene
    gene_set.append(df_new.iloc[0, 0])

    # Create a gene set of length 5, that are connected in the network
    while len(gene_set) < 5:
        # Find all genes that are connected to any gene in gene_set
        df_tmp_neigh = get_neighbourhood(gene_set, df_new)
        for gene in df_tmp_neigh['gene0']:
            if gene not in gene_set:
                next = gene
        gene_set.append(next)

    with open(data_filepath + "/gene_set.txt", 'w') as f:
        for item in gene_set:
            f.write("%s\n" % item)
