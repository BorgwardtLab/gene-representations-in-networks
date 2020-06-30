import argparse
import pandas as pd
import numpy as np
import json
from scipy.stats import combine_pvalues

def extract_snp_pval_from_id(id, df):
    series = df.loc[df['SNP'] == id, ['P']]
    print(id)
    print(series)
    if series['P'].empty:
        return " "
    print(series['P'])
    return series['P'].values[0]


def extract_df_value(row, column, df):
    # returns a Pandas Series
    series = df.loc[[row], [column]]
    # convert series to an actual value that can be used to extract matches
    return series[column].values[0]


def calculate_fisher_pval(list):
    return combine_pvalues(list, 'fisher')


def calculate_min_p_val(list):
    return np.min(list)


def extract_int_from_string(chromosome):
    return ''.join(filter(str.isdigit, chromosome))


def read_in_json_file(filepath):
    df_snp_mapping = pd.DataFrame(columns=['snp_id', 'gene_id'])
    with open(filepath) as json_file:
        data = json.load(json_file)
        for p in data:
            df_snp_mapping = df_snp_mapping.append({'snp_id': p, 'gene_id': data[p]}, ignore_index=True)
    return df_snp_mapping


def map_genes_with_snp_count(gene, df):
    map = {}
    number_of_snps = len(df)
    map[gene] = number_of_snps
    return map


def find_pval_representations(df_snp_mapping, df_snp_pval):
    # The column in which to find the gene_id and snp_id
    gene_id_column = 1
    snp_id_column = 0

    # Create a df for the p-value representation results
    df_mapping = pd.DataFrame(columns=['gene_id', 'av_p_val', 'min_p_val'])
    unique_genes = df_snp_mapping.iloc[:, gene_id_column].unique()
    gene_with_snp_count_dict = {}

    # Iterate over the unique genes and find the possible p-value representations
    for gene in unique_genes:
        # Find all SNP mappings for the unique gene
        snps_mappings_of_gene = df_snp_mapping[df_snp_mapping.iloc[:, gene_id_column].str.contains(gene)]
        gene_with_snp_count_dict[gene] = len(snps_mappings_of_gene)
        pval_list = []

        # Create a list with the p-values for the overlapping SNPs
        for index, snp in snps_mappings_of_gene.iterrows():
            snp_id = snp[snp_id_column]
            p_val = extract_snp_pval_from_id(snp_id, df_snp_pval)
            if p_val != " ":
                pval_list.append(p_val)

        # Calculate the three different representations for the gene
        if len(pval_list) != 0:
            stat, average_pval = calculate_fisher_pval(pval_list)
            min_pval = calculate_min_p_val(pval_list)

            # Add the representations to the dataframe
            df_mapping = df_mapping.append({'gene_id': gene, 'av_p_val': average_pval,
                                            'min_p_val': min_pval}, ignore_index=True)
    return df_mapping, gene_with_snp_count_dict


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--data_filepath')
    parser.add_argument('--pvalue_results')
    args = parser.parse_args()
    filepath = args.data_filepath
    pvalue_results = args.pvalue_results

    # Create a dataframe containing the GWAS results
    df_snp_pval = pd.read_csv(pvalue_results, delim_whitespace=True)
    df_mapping_total = pd.DataFrame(columns=['gene_id', 'av_p_val', 'min_p_val'])
    # Find the p-value representations for each chromosome and add the results
    for current_chromosome in range(1,6):
        # In order for the script to work this filepath needs to be set to the files created by mapping the genes to SNPs earlier
        current_filepath = './snp_mappings_chr' + str(current_chromosome) + '.txt'
        # Output files
        output_filepath = filepath + '/df_snp_mapping.txt'
        output_count_dictionary_path = filepath + '/number_of_snps_per_gene.txt'

        # Convert the gene-SNPs json file to a dataframe
        df_snp_mapping = read_in_json_file(current_filepath)
        #df_snp_mapping.to_csv('./df_snp_mappings_chr' + str(current_chromosome) + '.txt')

        # Get the rows for the current chromosome
        df_snp_pval_chromosome = df_snp_pval.loc[df_snp_pval['CHR'] == current_chromosome]

        # Find the p-value representations for the current chromosome and add to the total dataframe
	    df_mapping, count_dict = find_pval_representations(df_snp_mapping, df_snp_pval_chromosome)
	    df_mapping_total = df_mapping_total.append(df_mapping)
	    print(df_mapping_total)
	
        # Write the SNP count per chromosome to a file  
        with open(output_count_dictionary_path, 'r+') as f:
	      [f.write('{0} {1}\n'.format(key, value)) for key, value in count_dict.items()]
    # Write the total results from the p-value representations to a text file
    df_mapping_total.to_csv(output_filepath, index=False, sep=" ")



