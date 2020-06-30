import argparse
import pandas as pd
import json

# List with chromosome names as written in the data file
chromosome_list = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
snp_chromosome_dfs = []
gene_chromosome_dfs = []


def create_chromosome_dataframe(chromosome, df):
    return df[df['chromosome'] == chromosome]


def extract_int_char_from_string(chromosome):
    return ''.join(filter(str.isdigit, chromosome))


def sort_df_by_column(df, column):
    df.astype({column: 'int32'})
    return df.sort_values(by=[column])


def map_gene_and_snps_for_chromosome(df_gene, df_snp):
    # Create a local dictionary for the gene-snp mapping
    snp_gene_dict = {}

    # Iterate over each row in the gene df for the current chromosome
    for gene_index, gene in df_gene.iterrows():
        gene_start_pos = gene['start_pos']
        gene_end_pos = gene['end_pos']

        # Iterate over each snp row in the snp df for the current chromosome
        for snp_index, snp in df_snp.iterrows():
            snp_pos = snp['position']
            if snp_pos > gene_end_pos:
                break
            elif gene_start_pos <= snp_pos <= gene_end_pos:
                snp_id = snp['snp_id']
                gene_id = gene['gene_id']

                # Add the match to the dictionary with the mappings
                snp_gene_dict[snp_id] = gene_id
                # Delete the row of the found snp match from the snp dataframe
                df_snp = df_snp.drop([snp_index])
        print(gene_index)
    return snp_gene_dict


if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--filepath', help='Filepath to the athaliana data on the server')
    args = parser.parse_args()
    filepath = args.filepath

    # Dataframe with gene id:s, start and end positions
    df_total_gene = pd.read_csv(filepath + "/athaliana/gene_annotations/Araport11_GFF3_genes_transposons.201606.geneloc",
                                sep="\t", names=["gene_id", "chromosome", "start_pos", "end_pos"])
    # Dataframe with SNP id:s, chromosomes and SNP positions
    df_total_snp = pd.read_csv(filepath + "/athaliana/atwell_data/genotype.map", sep=" ",
                               names=["chromosome", "snp_id", "", "position"])

    for chromosome in chromosome_list:
        snp_chromosome_df = create_chromosome_dataframe(chromosome, df_total_snp)
        # Sort by the SNP positions
        snp_chromosome_dfs.append(sort_df_by_column(snp_chromosome_df, 'position'))

        gene_chromosome_df = create_chromosome_dataframe(extract_int_char_from_string(chromosome), df_total_gene)
        # Sort by the start position of the gene
        gene_chromosome_dfs.append(sort_df_by_column(gene_chromosome_df, 'start_pos'))

    dict_chr1 = map_gene_and_snps_for_chromosome(gene_chromosome_dfs[0], snp_chromosome_dfs[0])
    f = open("./snp_mappings_chr1.txt", "w")
    f.write(json.dumps(dict_chr1))
    f.close()

    dict_chr2 = map_gene_and_snps_for_chromosome(gene_chromosome_dfs[1], snp_chromosome_dfs[1])
    f = open('./snp_mappings_chr2.txt', "w")
    f.write(json.dumps(dict_chr2))
    f.close()

    dict_chr3 = map_gene_and_snps_for_chromosome(gene_chromosome_dfs[2], snp_chromosome_dfs[2])
    f = open("./snp_mappings_chr3.txt", "w")
    f.write(json.dumps(dict_chr3))
    f.close()

    dict_chr4 = map_gene_and_snps_for_chromosome(gene_chromosome_dfs[3], snp_chromosome_dfs[3])
    f = open("./snp_mappings_chr4.txt", "w")
    f.write(json.dumps(dict_chr4))
    f.close()

    dict_chr5 = map_gene_and_snps_for_chromosome(gene_chromosome_dfs[4], snp_chromosome_dfs[4])
    f = open("./snp_mappings_chr5.txt", "w")
    f.write(json.dumps(dict_chr5))
    f.close()
