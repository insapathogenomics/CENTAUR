#!/usr/bin/env	python3

"""
DownTree

By Veronica Mixao
@INSA
"""

import sys
import argparse
import textwrap
import datetime as datetime
import pandas
import numpy

version = "1.0.0"
last_updated = "2025-05-19"

# functions	----------

def print_log(message, log):
	""" print messages in the terminal and in the log file """

	print(message)
	print(message, file = log)

def check_arguments(args):
    """ check if all arguments provided are valid """

    if args.clusters == "":
        sys.exit("DownTree needs a file with clustering information. Please use the \"--clusters\" argument!")
    if args.metadata == "":
        sys.exit("DownTree needs a file with assembly information. Please use the \"--metadata\" argument!")
    if args.column_name == "":
        sys.exit("DownTree needs to know which column must be used to infer the assembly quality. Please use the \"--column_name\" argument!")
    if args.quality_threshold == "":
        sys.exit("DownTree needs to know which assembly quality threshold to use. Please use the \"--quality-threshold\" argument!")        

def get_threshold(mx, n, s, grouping_column, contigs_order, contigs_threshold, log):
    """ get the threshold providing the most similar number of partitions 
    to the desired final number of assemblies """
	
    thresholds = []
    selected_threshold = ""
    abs_difference = 1000000000
    check = True

    if contigs_order == "low":
        flt_mx = mx[mx[mx.columns[-1]] <= contigs_threshold]
    elif contigs_order == "high":
        flt_mx = mx[mx[mx.columns[-1]] >= contigs_threshold]
    mx = mx.reset_index(drop=True)
    flt_mx = flt_mx.reset_index(drop=True)

    if grouping_column != "":
        for col in mx.columns[1:-1]:
            if check:
                thresholds.append(col)
            if col == grouping_column:
                check = False
                mx[col].replace("", numpy.nan, inplace=True)
                mx.dropna(subset=[col], inplace=True)
                flt_mx[col].replace("", numpy.nan, inplace=True)
                flt_mx.dropna(subset=[col], inplace=True)
                partitions_expected = mx[col].values.tolist()
                partitions_real = flt_mx[col].values.tolist()
                expected_total_representatives = 0
                real_total_representatives = 0
                for partition_expected in set(partitions_expected):
                    partition_expected_n = partitions_expected.count(partition_expected)
                    if partition_expected_n <= int(s):
                        expected_total_representatives += partition_expected_n
                    else:
                        expected_total_representatives += int(s)
                for partition_real in set(partitions_real):
                    partition_real_n = partitions_real.count(partition_real)
                    if partition_real_n <= int(s):
                        real_total_representatives += partition_real_n
                    else:
                        real_total_representatives += int(s)
                expected_samples = expected_total_representatives
                expected_clusters = len(set(partitions_expected))
                real_samples = real_total_representatives
                real_clusters = len(set(partitions_real))
                selected_threshold = col
    
    else:
        for col in mx.columns[1:-1]:
            if check:
                mx[col].replace("", numpy.nan, inplace=True)
                mx.dropna(subset=[col], inplace=True)
                flt_mx[col].replace("", numpy.nan, inplace=True)
                flt_mx.dropna(subset=[col], inplace=True)
                partitions_expected = mx[col].values.tolist()
                partitions_real = flt_mx[col].values.tolist()
                expected_total_representatives = 0
                real_total_representatives = 0
                for partition_expected in set(partitions_expected):
                    partition_expected_n = partitions_expected.count(partition_expected)
                    if partition_expected_n <= int(s):
                        expected_total_representatives += partition_expected_n
                    else:
                        expected_total_representatives += int(s)
                for partition_real in set(partitions_real):
                    partition_real_n = partitions_real.count(partition_real)
                    if partition_real_n <= int(s):
                        real_total_representatives += partition_real_n
                    else:
                        real_total_representatives += int(s)
                abs_difference_thr = abs(n - real_total_representatives)
                if abs_difference_thr < abs_difference and real_total_representatives > int(n):
                    abs_difference = abs_difference_thr
                    selected_threshold = col
                    thresholds.append(col)
                    expected_samples = expected_total_representatives
                    expected_clusters = len(set(partitions_expected))
                    real_samples = real_total_representatives
                    real_clusters = len(set(partitions_real))
                else:
                    check = False

    print_log("\nThreshold selected for analysis: " + str(selected_threshold), log)
    print_log("\t*Expected " + str(expected_samples) + " samples representing " + str(expected_clusters) + " genomic groups (partitions).", log)
    print_log("Applying quality filtering...", log)
    print_log("\t*Estimates to select " + str(real_samples) + " samples representing " + str(real_clusters) + " (" + str(real_clusters/expected_clusters*100) + "%) genomic groups (partitions).", log)

    return thresholds, selected_threshold, expected_samples, expected_clusters

def get_samples(mx,thresholds,selected_threshold,s,contigs_col,contigs_order,contigs_threshold,log):
    """ select the assemblies """

    samples_choice = {}
    cluster_threshold = {}
    samples_selected = []
    cluster_len = {}
    mx[selected_threshold].replace("", numpy.nan, inplace=True)
    mx.dropna(subset=[selected_threshold], inplace=True)
    partitions = mx[selected_threshold].values.tolist()
    for partition in set(partitions):
        samples_choice[partition] = []
        flt_mx = mx[mx[selected_threshold] == partition]
        original_samples = flt_mx[flt_mx.columns[0]].values.tolist()
        cluster_len[partition] = len(original_samples)
        if contigs_order == "low":
            flt_mx = flt_mx.sort_values(contigs_col, ascending=True)
            flt_mx = flt_mx[flt_mx[contigs_col] <= contigs_threshold]
        elif contigs_order == "high":
            flt_mx = flt_mx.sort_values(contigs_col, ascending=False)
            flt_mx = flt_mx[flt_mx[contigs_col] >= contigs_threshold]
        samples = flt_mx[flt_mx.columns[0]].values.tolist()
        contigs = flt_mx[contigs_col].values.tolist()
        if len(samples) == 0:
            samples_choice[partition].append("NO SAMPLE PASSED THE QC CRITERIA")
        elif len(samples) <= int(s): # cluster with up to maximum samples per group
            for i in range(0,len(samples)):
                samples_choice[partition].append(samples[i])
                samples_selected.append(samples[i])
                cluster_threshold[partition] = selected_threshold
        else: # many samples and we need subsampling
            representatives = []
            if len(thresholds) >= 2:
                for w in range(2,len(thresholds)+1):
                    test_col = -w
                    sub_clusters = []
                    for val in flt_mx[thresholds[test_col]].values.tolist():
                        if val not in sub_clusters:
                            sub_clusters.append(val)
                    sub_cluster_len = {}
                    sub_cluster_comp = {}
                    for sub_cluster in sub_clusters:
                        sub_mx = flt_mx[flt_mx[thresholds[test_col]] == sub_cluster]
                        sub_samples = sub_mx[sub_mx.columns[0]].values.tolist()
                        sub_cluster_comp[sub_cluster] = sub_samples
                        sub_cluster_len[sub_cluster] = len(sub_samples)
                    sub_cluster_len = dict(sorted(sub_cluster_len.items(), key=lambda item: item[1], reverse=True))
                    if len(sub_clusters) > 1:
                        for sub_cluster in sub_cluster_len.keys():
                            sub_mx = flt_mx[flt_mx[thresholds[test_col]] == sub_cluster]
                            sub_samples = sub_mx[sub_mx.columns[0]].values.tolist()
                            if len(set(sub_samples) & set(representatives)) == 0 and len(representatives) < int(s): # we did not select a representative from this cluster yet
                                samples_choice[partition].append(sub_samples[0])
                                samples_selected.append(sub_samples[0])
                                cluster_threshold[partition] = thresholds[test_col]
                                representatives.append(sub_samples[0])
                    if w == len(thresholds) and len(representatives) < int(s):
                        j = 1
                        sub_clusters = set(sub_clusters)
                        exclude_subclusters = set()
                        length_representatives = len(representatives)
                        while length_representatives < int(s):
                            for sub_cluster in sub_cluster_len.keys():
                                sub_mx = flt_mx[flt_mx[thresholds[test_col]] == sub_cluster]
                                sub_samples = sub_mx[sub_mx.columns[0]].values.tolist()
                                if len(sub_samples) >= int(j+1):
                                    if sub_samples[j] not in representatives and length_representatives < int(s):
                                        representatives.append(sub_samples[j])
                                        samples_selected.append(sub_samples[j])
                                        cluster_threshold[partition] = thresholds[test_col]
                                        length_representatives = len(representatives)
                                    else:
                                        exclude_subclusters.add(sub_cluster)
                                else:
                                    exclude_subclusters.add(sub_cluster)
                            if len(exclude_subclusters) == len(sub_clusters):
                                length_representatives = int(s)
                            j += 1
            else:
                for i in range(0,int(s)):
                    representatives.append(samples[i])
                    samples_selected.append(samples[i])
                    cluster_threshold[partition] = selected_threshold
            samples_choice[partition] = representatives

    return samples_choice, cluster_threshold, samples_selected, cluster_len

def get_closest_represented_cluster(samples_choice, mx, selected_threshold, samples_selected):
    """ Determine, for each cluster or singleton not represented, what is the closest 
    cluster or singleton represented in the subdataset """

    represented_mx = mx[mx[mx.columns[0]].isin(samples_selected)]
    thr_failed_represented = {}
    i = 0
    for col in mx.columns:
        if col == selected_threshold:
               col_index = i
        i += 1
    for cluster in samples_choice.keys():
        if samples_choice[cluster][0] == "NO SAMPLE PASSED THE QC CRITERIA":
            flt_mx = mx[mx[selected_threshold] == cluster]
            for col_value in mx.columns[col_index+1:]:
                if cluster not in thr_failed_represented.keys():
                    cluster_name = flt_mx[col_value].values.tolist()[0]
                    clusters_represented = represented_mx[col_value].values.tolist()
                    if cluster_name in clusters_represented:
                        thr_failed_represented[cluster] = col_value
                else:
                    break
    return  thr_failed_represented
     
def get_results(selected_threshold, samples_choice, tag, cluster_threshold, thr_failed_represented, cluster_len, log):
    """ Determine in the end how many samples were selected and how many clusters are represented """

    observed_samples = []
    observed_clusters = []
    with open(tag + "_representatives_per_cluster.tsv", "w+") as outfile:
        print("selected_grouping_variable\tcluster_name\tcluster_size\tn_selected_samples\tselected_samples\tgrouping_variable_w_resolution\tthreshold_w_representation", file = outfile)
        for cluster in samples_choice.keys():
            if samples_choice[cluster][0] != "NO SAMPLE PASSED THE QC CRITERIA":
                observed_clusters.append(cluster)
                print(str(selected_threshold) + "\t" + str(cluster) + "\t" + str(cluster_len[cluster]) + "\t" + str(len(samples_choice[cluster])) + "\t" + ",".join(samples_choice[cluster]) + "\t" + cluster_threshold[cluster] + "\t" + str("-"), file = outfile)
                for sample in samples_choice[cluster]:
                    observed_samples.append(sample)
            else:
                 print(str(selected_threshold) + "\t" + str(cluster) + "\t" + str(cluster_len[cluster]) + "\t" + str(0) + "\t" + str("NO SAMPLE PASSED THE QC CRITERIA") + "\t" + str("-") + "\t" + str(thr_failed_represented[cluster]), file = outfile)
    
    print_log("\t*Selected " + str(len(observed_samples)) + " samples representing " + str(len(observed_clusters)) + " genomic groups (partitions).", log)
    
    return observed_samples, observed_clusters
        
# running the pipeline	----------

def main():
    
	# argument options	----------
    
    parser = argparse.ArgumentParser(prog="DownTree.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                                  DownTree.py                                  #
									#                                                                             #
									############################################################################### 
									
									          DownTree: Bioinformatics tool for Genome Diversity Extraction          
									
									This tool performs dataset downsampling based on its genetic diversity and on 
									the quality of the respective genome assemblies. Representative samples will be
									selected according to the dataset clustering patterns, as far as their 
									assemblies pass the quality control parameter defined by the user.
                                                                                        
									Note: This tool assumes that the columns with the clustering information are
									ordered from the highest resolution (i.e. more clusters) to the lowest.
									          
									-------------------------------------------------------------------------------"""))
	
    ## general parameters
	
    versioning = parser.add_argument_group("Version")
    versioning.add_argument("-v", "--version", dest="version", action="store_true", help="Print version and exit")
	
    group0 = parser.add_argument_group("Main options", "DownTree.py main options")
    group0.add_argument("-c", "--clusters", dest="clusters", default="", type=str, help="[MANDATORY] File with the clustering data in *.tsv format. The first column must correspond to the sample name, and the other(s) to \
                        the clustering classification. To take the most advantage of this tool, we recommend that you use the *_partitions.tsv file of ReporTree. Still, this tool will also work with other \
                        classifications, such as, for example, MLST, Serotype, etc.")
    group0.add_argument("-m", "--metadata", dest="metadata", default="", type=str, help="[MANDATORY] File with the assembly quality information in *.tsv format.")
    group0.add_argument("-g", "--grouping_column", dest="grouping_column", default="", type=str, help="[OPTIONAL] Column of the '--clusters' file that must be used to select samples based on genetic diversity. \
                        If this column is not specified, this tool will search for the column that provides the number of isolates that is closer to the number specified in '--number_assemblies'.")
    group0.add_argument("-n", "--number_assemblies", dest="number_assemblies", default=5000, type=int, help="Approximate number of representative assemblies to be selected from the whole dataset. If \
                        '--grouping_column' is specified, this argument will be ignored. (default: 5000)")
    group0.add_argument("-s", "--samples_per_group", dest="samples_per_group", default=3, type=int, help="Number of representative assemblies to be selected per cluster. (default: 3)")
    group0.add_argument("-col", "--column_name", dest="column_name", default="", type=str, help="[MANDATORY] Name of the column of the '--metadata' file that must be used to infer the assembly quality (data in this column \
                        must be numeric).")
    group0.add_argument("-qt", "--quality-threshold", dest="quality_threshold", default="", help="[MANDATORY] Maximum/minimum (depending on the '--order' argument) number in the column specified with \
                        '--column_name' that corresponds to the threshold between an elegible and a non-elegible sample for selection.")
    group0.add_argument("-order", "--order", dest="order", default="low", type=str, help="This argument can take one of two options: 'low' - indicating that lower values represent higher quality (e.g. contigs); or \
                        'high' - indicating that higher values represent higher quality (e.g. N50). (default: low)")
    group0.add_argument("-t", "--tag", dest="tag", default="DownTree", type=str, help="Tag for output file names. (default: DownTree)")
    group0.add_argument("--evaluation", dest="evaluation", action="store_true", help="Performs just the evaluation of the dataset and the clustering data. If you set this argument, no dataset will be generated.")
    args = parser.parse_args()

    # check if version	----------
	
    if args.version:
        print("version:", version, "\nlast_updated:", last_updated)
        sys.exit()   

    # check arguments	----------
	
    check_arguments(args)

    # start log ----------

    log_name = args.tag + ".log"
    log = open(log_name, "w+")

    print_log("\n******************** running DownTree.py ********************\n", log)
    print_log("version " + str(version) + " last updated on " + str(last_updated) + "\n", log)
    print_log(" ".join(sys.argv), log)
	
    start = datetime.datetime.now()
    print_log("start: " + str(start), log)

    # run pipeline	----------
    
    clustering_original = pandas.read_table(args.clusters, dtype = str)
    if args.grouping_column != "" and args.grouping_column not in clustering_original.columns:
        sys.exit("The indicated '--grouping_column' cannot be found in the clusters table!!")
    metadata = pandas.read_table(args.metadata)
    if args.column_name not in metadata.columns:
        sys.exit("The indicated '--column_name' cannot be found in the metadata table!!")
	
    clustering = clustering_original.set_index(clustering_original[clustering_original.columns[0]], drop = True)
    metadata = metadata.set_index(metadata[metadata.columns[0]], drop = True)
    df = pandas.concat([clustering, metadata[args.column_name]], axis=1, join="inner")
    
    # check if evaluation	----------
    
    thresholds, selected_threshold, expected_samples, expected_clusters = get_threshold(df, int(args.number_assemblies), int(args.samples_per_group), args.grouping_column, args.order, int(args.quality_threshold), log)
    if not args.evaluation:
        df2print = pandas.concat([clustering[clustering.columns[0]],clustering[selected_threshold], metadata[metadata.columns [1:]]], axis=1, join="inner")
        print_log("Selecting representative samples for each cluster...", log)
        samples_choice, cluster_threshold, samples_selected, cluster_len = get_samples(df,thresholds,selected_threshold,int(args.samples_per_group),args.column_name,args.order,int(args.quality_threshold),log)
        thr_failed_represented = get_closest_represented_cluster(samples_choice, clustering_original, selected_threshold, samples_selected)
        print_log("Working on reports...", log)
        observed_samples, observed_clusters = get_results(selected_threshold, samples_choice, args.tag, cluster_threshold, thr_failed_represented, cluster_len, log)
        selected_for_schema = []
        original_samples = df[df.columns[0]].values.tolist()
        for s in original_samples:
            if s not in observed_samples:
                selected_for_schema.append("no")
            else:
                selected_for_schema.append("yes")
        df2print.insert(2, "selected_for_subdataset", selected_for_schema)
        df2print.to_csv(args.tag + "_metadata_w_selection.tsv", index = False, header=True, sep ="\t")
        
        represented_in_schema = []
        for cluster_to_report in df2print[selected_threshold].values.tolist():
            if cluster_to_report not in thr_failed_represented.keys():
                represented_in_schema.append("yes")
            else:
                represented_in_schema.append(str(thr_failed_represented[cluster_to_report]))
        df2print.insert(3, "represented_in_subdataset", represented_in_schema)
        df2print.to_csv(args.tag + "_metadata_w_selection.tsv", index = False, header=True, sep ="\t")
         
    print_log("\nDone!", log)
    end = datetime.datetime.now()
    elapsed = end - start
    print_log("\nend: " + str(end) + "\n", log)
    print_log("Time elapsed: " + str(elapsed), log)
    print_log("\n*************************************************\n", log)
    log.close()

if __name__ == "__main__":
    main()