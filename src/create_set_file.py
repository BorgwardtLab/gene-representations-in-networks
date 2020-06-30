# Create a set file for PLINK's set-based method where each gene creates a set with the overlapping SNPs in the format:
#
# gene-id1
# SNP-id 1
# SNP-id 2
# END
#
# gene-id2
# SNP-id 1
# SNP-id 2
# END
#
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_filepath')
    args = parser.parse_args()
    filepath = args.data_filepath

    df_gene_snps = pd.read_csv(filepath + '/df_snp_mapping.txt', sep=" ", names=['gene', 'snp'])
    with open(filepath + '/gene_snp_set_file.txt', 'w') as outfile:
        for gene in df_gene_snps.iterrows():
            gene_id = gene[1].values[0]
            snps = df_gene_snps.loc[df_gene_snps['gene'] == gene_id]
            outfile.write(gene_id + "\n")
            for snp in snps.iterrows():
                outfile.write(str(snp[1].values[1]) + "\n")
            outfile.write("END\n")
            outfile.write("\n")

