# A systematic evaluation of gene representations in network based genetic analysis

The purpose of the semester project is to explore four different p-value representations for a gene in a hierarchical HotNet analysis. 
The p-value representations used are: Minimum, Fisher, SKAT-O and PLINK's set-based tests. The data used in this project comes from the *Athaliana* dataset created by Atwell et al., and the protein-protein network used is from the resource TAIR. On the server this data is located under "/links/groups/borgwardt/Projects/master_project_julia/data/athaliana/atwell_data/"

In order to analyze the importance of the p-value representation the following three steps are executed:
1. Map the SNPs to genes with overlapping locations. 
2. Obtain the p-value representations for the genes using the four different methods.
3. Perform a hierarchical HotNet analysis.

Step one and two are combined in the bash script src/bash scripts/get_gene_pval_representations.sh,
which can be executed on the server after loading the correct version of python by running the command: "module load local/python_bleeding_edge".
The varaibles that need to be set are the filepath to the phenotype to be used in the association study, 
"$FILTERED_PHENO", which in the example file is set to be the flowering phenotype FRI.

The script outputs three files:
1. Containing the p-value representations of Minimum and Fisher's method for the data.
2. Containing the p-value representations for SKAT-O. 
3. Containing the p-value representations for the Set-based method.

Thereafter three hierarchical HotNet input files needs to be created (more information are available in the readme for the source code https://github.com/raphael-group/hierarchical-hotnet/). 
This can be done by using the src/prepare_hotnet_files.py. The p-value representation to be used are set in the script, therefore in order to obtain
one for each method, the script will need to be run four times. When the input is created a HotNet analysis can be executed as in the example file src/bash scripts/run_hotnet.sh. 
Make sure to set the input data directory to the input files created for the hierarchical HotNet analysis. 

Another example on how to execute the hierarchical HotNet analysis using our own phenotype permutations (created via permute_phenotype.py and running association studies for each phenotype permutation via run_gwas_on_all_permutations.sh) can be found in run_hotnet_own_permutations.sh. 

In the project a simulation study was executed by extracting 5 genes that were connected in the network and had 5-10 SNPs each. This was done by using the script: src/bash script/generate_simulated_phenotype.sh. 
(In the Python file generate_simulated_phenotype.py the strength of the association can be specified, the example is set to a weak correlation). 

In the directory simulations a number of scripts are available in order to explore the data (in particular for the simulated dataset). 
For example src/simulations/plot_cluster_distributions_hotnet.py plots the distribution of cluster sizes for the different methods.