#!/usr/bin/env	python3

"""
PanTree: Mapping the Occurrence of Unique and Shared wgMLST loci Across a species phylogenetic tree 

By Veronica Mixao
@INSA
"""

import sys
import argparse
import datetime as datetime
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

version = "1.0.0"
last_updated = "2025-05-20"

def print_log(message, log):
	""" print messages in the terminal and in the log file """

	print(message)
	print(message, file = log)

def generate_binary(input_file, output_file):
    """ This function loads the allele matrix and returns a binary matrix """
    
    allelic_matrix = pd.read_csv(input_file, sep="\t")
    allelic_matrix[allelic_matrix.columns[0]] = allelic_matrix[allelic_matrix.columns[0]].astype(str)
    allelic_matrix = allelic_matrix.set_index(allelic_matrix.columns[0])

    # Convert all non-zero values to 1 (present) and keep 0 as absent
    binary_matrix = allelic_matrix.copy()
    binary_matrix[binary_matrix != 0] = 1

    binary_matrix.to_csv(output_file, index = True, header=True, sep ="\t")

    return binary_matrix, allelic_matrix

def order_loci(binary_matrix, allelic_matrix, tag_output, order, thr):
    """ This function orders loci in an allele matrix according to their prevalence """
    mx_percentage = pd.DataFrame()
    mx_percentage["Loci"] = binary_matrix.columns
    presence_percentage = binary_matrix.mean(axis=0)
    mx_percentage["Prevalence"] = presence_percentage.values
    mx_percentage.to_csv(f"{tag_output}_loci_prevalence.tsv", index = False, header=True, sep ="\t")
    if thr != "":
        mx_percentage = mx_percentage[mx_percentage["Prevalence"] < float(thr)]
        final_loci = mx_percentage["Loci"].values.tolist()
        binary_matrix = binary_matrix[final_loci]
        allelic_matrix = allelic_matrix[final_loci]
    if order == "yes":
        presence_percentage_flt = presence_percentage[presence_percentage.index.isin(binary_matrix.columns)]
        binary_matrix = binary_matrix.loc[:, presence_percentage_flt.sort_values(ascending=False).index]
        allelic_matrix = allelic_matrix.loc[:, presence_percentage_flt.sort_values(ascending=False).index]

    return binary_matrix, allelic_matrix, mx_percentage

def order_samples(binary_matrix, allelic_matrix, tree_file, tag_output):
    """ This function orders samples in an allele matrix according to their presence in a tree """
    
    tree = Phylo.read(tree_file, "newick")
    sample_order = [str(clade) for clade in tree.get_terminals()]
    
    # Ensure only common samples between the tree and the input matrix are considered
    common_samples = []
    for sample in sample_order:
        if sample in binary_matrix.index:
            common_samples.append(sample)
        
    if not common_samples:
        sys.exit("\nERROR!! Warning: No common samples between the tree and the input matrix.")

    binary_matrix = binary_matrix.reindex(common_samples)
    binary_matrix.to_csv(f"{tag_output}_binary_ordered.tsv", index = True, header=True, sep ="\t")
    allelic_matrix = allelic_matrix.reindex(common_samples)
    allelic_matrix.to_csv(f"{tag_output}_alleles_ordered.tsv", index = True, header=True, sep ="\t")

    return binary_matrix, allelic_matrix

def prevalence_heatmap(binary_matrix, output, order):
    """ This function generates a heatmap with the tree and the loci prevalence """

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(binary_matrix, cmap=["white", "blue"], cbar=False, xticklabels=False, yticklabels=False, ax=ax)

    y_labels_fontsize = max(4, 12 - len(binary_matrix) // 10)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=y_labels_fontsize)

    blue_patch = mpatches.Patch(color='blue', label='Present')
    white_patch = mpatches.Patch(color='white', label='Absent')
    plt.legend(handles=[blue_patch, white_patch], loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    # Add titles and labels
    plt.title("Heatmap of Loci Presence/Absence")
    xlabel = "Loci (Sorted by Presence Percentage)" if order == "yes" else "Loci"
    plt.xlabel(xlabel)
    plt.ylabel("Samples")

    # Save the figure to the specified path
    plt.tight_layout()
    plt.savefig(f"{output}_loci_prevalence.png", dpi=300, bbox_inches='tight')

def get_loci_prevalence_threshold(threshold, allelic_matrix, partitions, presence_percentage, tag, columns2check):
    """ This function summarizes the loci prevalence and diversity within a cgMLST threshold """

    cluster_column = f"cluster_{threshold}"
    allele_div = {cluster_column: []}
    loci_prev = {cluster_column: []}
    summary = presence_percentage
    loci_lst = summary["Loci"].values.tolist()
    clusters_series = partitions[str(threshold)]
    clusters = clusters_series.unique()

    cluster_to_samples = partitions.groupby(str(threshold))[partitions.columns[0]].apply(set).to_dict()

    n_alleles_lst = []
    n_clusters = []
    n_clusters_present = []
    n_clusters_present_pct = []

    for loci in loci_lst:
        if loci not in allele_div.keys():
            allele_div[loci] = []
            loci_prev[loci] = []
        
        locus_mx = allelic_matrix[allelic_matrix[loci].astype(str) != "0"]
        locus_samples = set(locus_mx.index)
        n_alleles = allelic_matrix[loci].nunique()
        n_alleles_lst.append(n_alleles)
        count_clusters_threshold = 0
        count_clusters_locus = 0

        for cluster, samples in cluster_to_samples.items():
            if len(samples) <= 1:
                continue
            else:
                if cluster not in allele_div[cluster_column]:
                    allele_div[cluster_column].append(cluster)
                    loci_prev[cluster_column].append(cluster)
                count_clusters_threshold += 1
                overlap = samples & locus_samples
                if overlap:
                    count_clusters_locus += 1
                    locus_mx_cluster = allelic_matrix.loc[list(overlap), loci]
                    n_alleles_locus_cluster = locus_mx_cluster.nunique()
                    prevalence_locus_cluster = len(overlap) / len(samples)
                else:
                    n_alleles_locus_cluster = 0
                    prevalence_locus_cluster = 0

                allele_div[loci].append(n_alleles_locus_cluster)
                loci_prev[loci].append(prevalence_locus_cluster)

        n_clusters.append(count_clusters_threshold)
        n_clusters_present.append(count_clusters_locus)
        n_clusters_present_pct.append(count_clusters_locus / count_clusters_threshold if count_clusters_threshold else 0)
        
    summary[f"N_alleles"] = n_alleles_lst
    summary[f"N_clusters_{threshold}"] = n_clusters
    summary[f"Prevalence_clusters_{threshold}"] = n_clusters_present
    summary[f"Prevalence_clusters_{threshold}_pct"] = n_clusters_present_pct
    columns2check.append(f"Prevalence_clusters_{threshold}_pct")

    pd.DataFrame(allele_div).to_csv(f"{tag}_{threshold}_alleles.tsv", index=False, sep="\t")
    pd.DataFrame(loci_prev).to_csv(f"{tag}_{threshold}_prevalence.tsv", index=False, sep="\t")

    return summary, columns2check

def pass_loci(summary, columns2check, thr):
    """ this function adds a new column informing the user of the threshold at which each locus passes the minimum cluster prevalence required """

    def check_row(row):
        PASS = [col for col in columns2check if row[col] >= float(thr)]
        final_names = []
        for v in PASS:
            v = v.split("Prevalence_clusters_")[1]
            v = v.split("_pct")[0]
            final_names.append(v)
        return ",".join(final_names) if PASS else "-"

    summary["Selected"] = summary.apply(check_row, axis=1)

    return summary 

def get_filtered_binary(binary, presence_percentage, thr, column, direction):
    """ filter the binary matrix according to the cluster prevalenceof each locus """
    
    presence_percentage = presence_percentage.sort_values(by=column, ascending=False)
    if direction == "higher":
        loci_list = presence_percentage.loc[presence_percentage[column] >= float(thr), presence_percentage.columns[0]].tolist()
    elif direction == "lower":
        loci_list = presence_percentage.loc[presence_percentage[column] < float(thr), presence_percentage.columns[0]].tolist()
    binary = binary[loci_list]
    
    threshold_name = column.split("Prevalence_clusters_")[1]
    threshold_name = threshold_name.split("_pct")[0]
    print("\tNumber of loci " + str(threshold_name) + ": " + str(len(loci_list)))

    return binary, threshold_name

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="PanTree - Evaluate loci distribution across a species tree")
    parser.add_argument("-a", "--allele", dest="allele", required=True, help="[MANDATORY] Path to the input allele matrix (TSV). Missing data must be represented by 0.")
    parser.add_argument("-p", "--partitions", dest="partitions", required=False, help="[MANDATORY] Path to the input partitions table (TSV). We recomend using the partitions table obtained with ReporTree.")
    parser.add_argument("-tree", "--tree", dest="tree", help="[MANDATORY] Path to the phylogenetic tree (NWK). It can be a single-linkage tree obtained with ReporTree.")
    parser.add_argument("-s", "--column2select", dest="column_select", required=True, help="[MANDATORY] Comma-separated list of the thresholds (column names) of the partitions table that will be used to \
                        flag loci dispersed throughout the tree (e.g. --column2select single-1x1.0,single-10x1.0,single-100x1.0).")
    parser.add_argument("-cp", "--cluster-prevalence", dest="cluster_prevalence", required=True, help="[MANDATORY] Percentage of clusters/groups were a given locus is present, so it \
                        can be considered as elegible for the wgMLST.")
    parser.add_argument("-o", "--output", dest="output", required=True, help="[MANDATORY] Tag for output files.")
    parser.add_argument("--order", dest="order", choices=["yes", "no"], default="no", help="Order loci by the percentage of presence - descending [default: no].")
    parser.add_argument("--threshold", dest="thr", default="", help="Only analyze loci with a prevalence lower than this threshold (value between 0 and 1).")

    # Parse arguments
    args = parser.parse_args()

    # starting logs	----------
    
    log_name = f"{args.output}.log"
    log = open(log_name, "w+")

    print_log("\n******************** running PanTree.py ********************\n", log)
    print_log(f"Version {version} (Last updated {last_updated})", log)
    print_log("Command: " + " ".join(sys.argv), log)
	
    start = datetime.datetime.now()
    print_log(f"Start: {start}\n", log)

    print_log("Generating binary matrix...", log)
    binary_matrix, allelic_matrix = generate_binary(args.allele, f"{args.output}_binary.tsv")

    print_log("Ordering loci...", log)
    binary_matrix, allelic_matrix, presence_percentage = order_loci(binary_matrix, allelic_matrix, args.output, args.order, args.thr)

    print_log("Ordering samples...", log)
    binary_matrix, allelic_matrix = order_samples(binary_matrix, allelic_matrix, args.tree, args.output)

    print_log("Generating the prevalence heatmap...", log)
    prevalence_heatmap(binary_matrix, args.output, args.order)

    print_log("Determining loci prevalence and diversity per threshold...", log)
    partitions = pd.read_table(args.partitions)
    thresholds4selection = str(args.column_select).split(",")
    columns2check = []
    for threshold in thresholds4selection:
        print_log(f"\tProcessing threshold: {threshold}", log)
        presence_percentage, columns2check = get_loci_prevalence_threshold(threshold, allelic_matrix, partitions, presence_percentage, args.output, columns2check)
    
    print_log("Finding loci passing the cluster prevalence criteria...", log)
    presence_percentage.to_csv(f"{args.output}_summary.tsv", index = False, header=True, sep ="\t")
    
    presence_percentage = pass_loci(presence_percentage, columns2check, args.cluster_prevalence)
    
    presence_percentage.to_csv(f"{args.output}_summary.tsv", index = False, header=True, sep ="\t")

    print_log("Generating prevalence heatmap for each loci list...", log)
    for column in columns2check:
        binary, threshold_name = get_filtered_binary(binary_matrix, presence_percentage, args.cluster_prevalence, column, "higher")
        prevalence_heatmap(binary, args.output + "_" + str(threshold_name), args.order)
    
    print_log("Generating prevalence heatmap for each list of non-selected loci...", log)
    for column in columns2check:
        binary, threshold_name = get_filtered_binary(binary_matrix, presence_percentage, args.cluster_prevalence, column, "lower")
        prevalence_heatmap(binary, args.output + "_removed_" + str(threshold_name), args.order)


    # done	----------
	
    end = datetime.datetime.now()
	
    elapsed = end - start
    print_log("\n------------------------------------------------------------\n", log)
    print_log("PanTree.py is done! If you found any issue please contact us!!\n", log)
    print_log(f"End: {end}", log)
    print_log(f"Total elapsed time: {elapsed}", log)
    log.close()

if __name__ == "__main__":
    main()
