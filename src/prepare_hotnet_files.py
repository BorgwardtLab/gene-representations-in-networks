import pandas as pd
import numpy as np

# Create dataframe with the edges present in the network
df_gene_edges = pd.read_csv("./TairProteinInteraction.20090527_edges_sym.txt", delim_whitespace=True, header=None)
# Create dataframe with the SKAT-O results
df_SKATO_pvalues = pd.read_csv("./SKATO_test_results.txt", delim_whitespace=True, usecols=["SetID", 'P.value'])
# Create dataframe with the p-value representations of Fisher and Minimum
df_pvalues = pd.read_csv("./p-val_representations.txt")
# Create dataframe with the genes to be used
df_genes = df_SKATO_pvalues["SetID"]

# Create dataframes for the three files that will be generated as HotNet inputs
df_gene_edges_indices = pd.DataFrame(columns=['index0', 'index1'])
df_indices_genes = pd.DataFrame(columns=['index', 'gene'])
df_gene_scores = pd.DataFrame(columns=['gene', 'score'])

if __name__ == '__main__':
    # Choose here which p-value representation to use for the score file
    method = 'min_p_val'
    # -log10 transform the p-values before using them as scores in the network
    df_pvalues[method] = -np.log10(df_pvalues[method])
    # If running for SKAT-O method use the commented lines instead
    #df_SKATO_pvalues['P.value'] = -np.log10(df_SKATO_pvalues['P.value'])

    # Create the gene-index mapping from the genes present in SKAT-O
    gene_index_map = {}
    index = 1
    for gene in df_genes:
        gene_index_map[gene] = index
        index += 1

    new_gene_index_map = {}
    i = 1

    # Go through each edge in the network
    for gene in df_gene_edges.iterrows():
        gene0 = gene[1].values[0]
        gene1 = gene[1].values[1]
        # Extract the gene indices if the genes are present in the gene_index_map
        index0 = gene_index_map.get(gene0)
        index1 = gene_index_map.get(gene1)

        # Only add the edge if both genes exist in the gene-index map
        if index0 is None or index1 is None:
            print('Removes edge, one of the genes are missing')
        # Both genes exists, re-index the gene such that values goes from 0 and upwards by creating a new map
        else:
            i1 = new_gene_index_map.get(gene0)
            i2 = new_gene_index_map.get(gene1)

            # Only add the gene the first time it is missing
            if i1 is None:
                i1 = i
                new_gene_index_map[gene0] = i1
                df_indices_genes = df_indices_genes.append({'index': i1, 'gene': gene0},
                                                           ignore_index=True)
                df_gene_scores = df_gene_scores.append({'gene': gene0, 'score':
                    #df_SKATO_pvalues.loc[df_SKATO_pvalues['SetID'] == gene0, ['P.value']].values[0].item()},
                    #                                   ignore_index=True)
                    df_pvalues.loc[df_pvalues['SET'] == gene0, [method]].values[0].item()},
                                                       ignore_index=True)
                i += 1
            if i2 is None:
                i2 = i
                new_gene_index_map[gene1] = i2
                df_indices_genes = df_indices_genes.append({'index': i2, 'gene': gene1},
                                                           ignore_index=True)
                df_gene_scores = df_gene_scores.append({'gene': gene1, 'score':
                    #df_SKATO_pvalues.loc[df_SKATO_pvalues['SetID'] == gene1, ['P.value']].values[0].item()},
                    #                                   ignore_index=True)
                    df_pvalues.loc[df_pvalues['SET'] == gene0, [method]].values[0].item()},
                                                       ignore_index=True)
                i += 1

            df_gene_edges_indices = df_gene_edges_indices.append({'index0': i1, 'index1': i2}, ignore_index=True)


    # Write the output files
    df_indices_genes.sort_values(by=['index'])
    # Change so the columns do not contain the weirdly formatted 0 values
    df_gene_scores =  df_gene_scores.mask(df_gene_scores == -0.0000000, 0)
    df_indices_genes.to_csv('./network_1_index_gene.tsv', header=None, sep=" ", index=False)
    df_gene_edges_indices.to_csv('./network_1_edge_list.tsv', header=None, index=False, sep=" ")
    df_gene_scores.to_csv('./score_1.tsv', header=None, index=False, sep=" ")

