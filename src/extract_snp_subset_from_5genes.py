import pandas as pd
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--data_filepath')
    args = parser.parse_args()
    filepath = args.data_filepath
	
    df_gene_subset = pd.read_csv(filepath + '/gene_set.txt', delim_whitespace=True, header=None)
    df_gene_snp = pd.read_csv('./gene_snp_mapping.txt', delim_whitespace=True, names=['gene_id', 'snp_id'])
    df_new = pd.DataFrame()    
    #df = pd.DataFrame(columns=['gene_id', 'number'])
    for gene in df_gene_subset.iterrows():
        gene_id = gene[1].values[0]
        print(gene_id)
        snps = df_gene_snp.loc[df_gene_snp['gene_id'] == gene_id, :]
        #nsnps = len(snps	)
	print(snps)
        #df = df.append({'gene_id': gene_id, 'number': nsnps}, ignore_index=True)
        df_new = df_new.append(snps)
    df_snps = df_new.iloc[:,1]
    #df.to_csv('./number_of_snps_per_gene.txt')
    #df_new.to_csv('../permutation_data/gene_snp_subset_5genes.txt', index=False, header=None, sep=" ")
    print(df_snps)
    df_snps.to_csv(filepath + '/snp_set.txt', index=False, header=None, sep=" ")


