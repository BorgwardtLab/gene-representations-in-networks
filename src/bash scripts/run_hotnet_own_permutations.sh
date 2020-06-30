#!/usr/bin/env bash
# Set the paths to the directories, data contains the input and the others will store the HotNet data
data=/links/groups/borgwardt/Projects/master_project_julia/data/permutation_data/permutation_data/permutation_gwas/strong_min/data_permutation_min
intermediate=/links/groups/borgwardt/Projects/master_project_julia/data/permutation_data/permutation_data/permutation_gwas/strong_min
results=/links/groups/borgwardt/Projects/master_project_julia/data/permutation_data/permutation_data/permutation_gwas/strong_min/results_simulated_permutation_data_min

# Number of permutations performed on the data
num_permutations=500

################################################################################
#
#   Prepare data.
#
################################################################################

# Create the intermediate and results directory
mkdir -p $results

for network in network_1
do
    mkdir -p $intermediate/"$network"
done

for network in network_1
do
    for score in scores_1
    do
        mkdir -p $intermediate/"$network"_"$score"
    done
done

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

echo "Construct similarity matrices..."

for network in network_1
do
    # Set python paths to the src files from the hierarchical HotNet repository
    python /links/groups/borgwardt/Projects/master_project_julia/code/hierarchical_hotnet/src/construct_similarity_matrix.py \
        -i   $data/"$network"_edge_list.tsv \
        -o   $intermediate/similarity_matrix.h5 \
        -bof $intermediate/beta.txt
done


##################################################################################
#
#   Construct hierarchies.
#
##################################################################################

echo "Constructing hierarchies..."

for network in network_1
do
    for score in score_1
    do
        cp $data/score_1.tsv $intermediate/score_0.txt

        for i in `seq 0 $num_permutations`
        do
            python /links/groups/borgwardt/Projects/master_project_julia/code/hierarchical_hotnet/src/construct_hierarchy.py \
                -smf  $intermediate/similarity_matrix.h5 \
                -igf  $data/"$network"_index_gene.tsv \
                -gsf  $intermediate/score_"$i".txt \
                -helf $intermediate/h_edge_list_"$i".tsv \
                -higf $intermediate/h_index_gene_"$i".tsv
        done
    done
done

#################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 10 (because default is 10) for larger graphs.

for network in network_1
do
    for score in scores_1
    do
        python /links/groups/borgwardt/Projects/master_project_julia/code/hierarchical_hotnet/src/process_hierarchies.py \
            -oelf $data/"$network"_edge_list.tsv \
            -oigf $data/"$network"_index_gene.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/h_edge_list_"$i".tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/h_index_gene_"$i".tsv "; done) \
            -lsb  10 \
            -cf   $results/clusters_"$network"_"$score".tsv \
            -pl   $network $score \
            -pf   $results/sizes_"$network"_"$score".pdf
    done
done