#!/usr/bin/env bash

# Example of setting the files on the local computer
data=$PWD/data_permutation_fisher
intermediate=$PWD/intermediate_permutation_data_fisher
results=$PWD/results_simulated_permutation_data_fisher

# Compile Fortran module.
cd ../src
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

# Create data, intermediate data and results and results directories.
mkdir -p $intermediate
mkdir -p $intermediate/network_1
mkdir -p $intermediate/network_1_score_1
mkdir -p $results

# Number of permutations the HotNet algorithm will create
num_permutations=100

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

echo "Construct similarity matrices..."

for network in network_1
do
    python /Users/juliaortheden/PycharmProjects/master_project/MLCB_project/hierarchical-hotnet/src/construct_similarity_matrix.py \
        -i   $data/"$network"_edge_list.tsv \
        -o   $intermediate/similarity_matrix.h5 \
        -bof $intermediate/beta.txt
done


###############################################################################
#
# Permute data
#
###############################################################################

echo "Permuting scores..."

for network in network_1
do
    for score in score_1
    do
    cp $data/score_1.tsv $intermediate/network_1_score_1/scores_0.tsv

        python /Users/juliaortheden/PycharmProjects/master_project/MLCB_project/hierarchical-hotnet/src/find_permutation_bins.py \
            -gsf $intermediate/"$network"_"$score"/scores_0.tsv \
            -igf $data/"$network"_index_gene.tsv \
            -elf $data/"$network"_edge_list.tsv \
            -ms  100 \
            -o   $intermediate/"$network"_"$score"/score_bins.tsv

        for i in `seq $num_permutations`
        do
            python /Users/juliaortheden/PycharmProjects/master_project/MLCB_project/hierarchical-hotnet/src/permute_scores.py \
                -i  $intermediate/"$network"_"$score"/scores_0.tsv \
                -bf $intermediate/"$network"_"$score"/score_bins.tsv \
                -s  "$i" \
                -o  $intermediate/"$network"_"$score"/scores_"$i".tsv
        done
    done
done

###############################################################################
#
# Construct hierarchies.
#
###############################################################################

echo "Constructing hierarchies..."

for network in network_1
do
    for score in score_1
    do
        for i in `seq 0 $num_permutations`
        do
            python /Users/juliaortheden/PycharmProjects/master_project/MLCB_project/hierarchical-hotnet/src/construct_hierarchy.py \
                -smf  $intermediate/similarity_matrix.h5 \
                -igf  $data/"$network"_index_gene.tsv \
                -gsf  $intermediate/"$network"_"$score"/scores_"$i".tsv \
                -helf $intermediate/"$network"_"$score"/h_edge_list_"$i".tsv \
                -higf $intermediate/"$network"_"$score"/h_index_gene_"$i".tsv
        done
    done
done

################################################################################
#
# Process hierarchies.
#
###############################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 10 (because default is 10) for larger graphs.

        python /Users/juliaortheden/PycharmProjects/master_project/MLCB_project/hierarchical-hotnet/src/process_hierarchies.py \
            -oelf $intermediate/network_1_score_1/h_edge_list_0.tsv \
            -oigf $intermediate/network_1_score_1/h_index_gene_0.tsv \
            -pelf $(for i in `seq $num_permutations`; do echo " $intermediate/network_1_score_1/h_edge_list_"$i".tsv "; done) \
            -pigf $(for i in `seq $num_permutations`; do echo " $intermediate/network_1_score_1/h_index_gene_"$i".tsv "; done) \
            -lsb  10 \
            -cf   $results/clusters_network_score_1.tsv \
            -pl   $network $score \
            -pf   $results/sizes_network_score.pdf
