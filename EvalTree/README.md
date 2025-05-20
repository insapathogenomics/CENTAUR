# EvalTree
Toolbox for comparative clustering evaluation of Whole-Genome Sequencing (WGS) pipelines for routine bacterial surveillance.
Evaltree is a tool developed within the Centaur framework that aims to promote integrated epidemiological surveillance in a One-Health programme.
Thus, contributing to a more efficient response to outbreak scenarios and preventing the spread of infectious diseases, any bacterial pathogen X.

It can be useful in investigating at least four scenarios in genomic surveillance:

- Compare the performance of WGS pipelines between laboratories.

- Compare the resolution between well-established schemes and new ones.

- Promoting the long-term sustainability of pipelines (updates, versions, parameters).

- Comparison between cg/wgMLST schemes and traditional typing methods (serotypes, sequence-type, clonal complex).

EvalTree aims to support laboratories by:

- Automating routine tasks, reducing the risk of error and saving time.

- Simplifying the evaluation of congruence between two pipelines, eliminating the need to run scripts individually, and promoting reproducibility and confidence in the results.

- Generating reports with interactive plots, allowing users explore the data. This makes interpretation less time-consuming and demanding, while also facilitating data sharing.

- Contributing to the standardization of results, promoting communication between national reference laboratories and other EU members, and strengthening collaboration and result comparability.

- Being easily implemented by laboratories with lower expertise and to have a simple installation process.


## Implementation

EvalTree is a command-line tool implemented in Python 3.8, able to create user-friendly reports by orchestrating the scripts developed by Mixão et. al., 2025. 
It is highly recommended that inputs be originated using ReporTree software available in [insapathogenomics/ReporTree](https://github.com/insapathogenomics/ReporTree) to ensure optimal functionality. 

## Inputs

EvalTree is a tool designed primarily to make comparisons between two inputs.
The inputs can be folders (a ReporTree folder is recommended) or files (partition files or other types of files with classifications (e.g., sequence-type, serotypes, etc.) are recommended).
When ReporTree folders are provided, the tool can process 5 types of files (described below), according to the instructions provided on the command line. 
Optionally, only one ReporTree folder can be entered, but only the functionalities of visualising clustering plots resulting from partitions_summary or sample_of_interest and the pipeline characterisation section are available.

- Description files:
    - partitions.tsv: The main input. It contains all genetic partitions (clusters) across each allele distance threshold.
      
    - clusterComposition.tsv: File that categorizes samples as part of a cluster or as singletons for each threshold, and reports the size of each cluster and the singleton.
      
    - partitions_summary.tsv: Describe per threshold and cluster associated metadata (e.g., country, source). File also report changes in the clusters across the thresholds, including: nomenclature change (indicates if a cluster has changed compared to the previous partition (kept, split, merged, etc.), n_increases (indicating the number of new samples added to the cluster), and size of the cluster. Nomenclature change is better explained in [insapathogenomics/ReporTree](https://github.com/insapathogenomics/ReporTree)
      
    - sample_of_interest_partitions_summary.tsv:  Similar to partitions_summary.tsv, but exclusively for the samples of interest.
      
    - stableRegions.tsv: Identify stability regions in each pipeline, where the cluster composition is stable across five consecutive thresholds and with an Adjusted Wallace coefficient above 0.99.

## Main output files

An HTML-based report is generated, containing interactive Plotly visualisations built from the accessory files produced during script orchestration.
The document adopts a hierarchical accordion structure, organized into two main sections (characterization and congruence) and two optional sections (clustering and outbreaks).

- Pipeline characterization:
    - Displays the executed command line and the name of the pipeline.
    - Shows the number of samples analysed and the thresholds used.
    - Includes a scarlet plot illustrating the number of partitions per threshold.

- Congruence between two pipelines:
    - Compares the clustering results of pipelines across all resolution levels, visualized using a heatmap. 
    - A bar plot of stability regions is also included.
    - A tendency plot highlights the best correspondence point that produces the highest agreement thresholds between the pipelines. 

- Clustering plots:
    - For each pipeline, clustering plots are displayed for the selected threshold(s) (--plots-threshold).
    - Results can be split by one or more categories (--column-plots).
    - The pie chart plots include one or multiple segmentes, where each segment represents the frequency of samples within a given category.

- Outbreak plots:
    - Evaluates cluster composition for a given threshold (--threshold_outbreak) and pipeline, illustrated in a heatmap.

Other important file:
- All_correspondence.tsv: It allows identification of corresponding thresholds between the methods used in each pipeline.
  
## Description of the main modules of EvalTree

  ## Characterization
  - Requires partitions.tsv file generated by ReporTree software. 
  - The number of existing partitions is calculated for each allelic distance threshold, reflecting how the samples are grouped at a given threshold.
    
  ## Clustering visualisation
  - Clustering is presented when one or more ReporTree folders are provided as input, each containing one of the files partitions_summary or sample_of_interest.
  - The user must provide the following arguments for this section to be included in the report:
    -  plots_threshold (-pt; allelic distance threshold)
    -  columns_plots (-cp; defines which categories will be analysed)
  - By default, the program assumes that the file is partitions_summary and sets n_cluster (number of plots) to 3.
  - If the file is sample_of_interest, the argument plots_summary (-ps) must be provided with the option sample_of_interest. In this case, only the plots where the cluster size increases will be shown. The n_cluster argument is not applicable in this case.

  ## Congruence
   - Requires partitions.tsv file generated by ReporTree software.
  - The congruence of two pipelines is evaluated from the orchestration of 2 scripts,  comparing_partitions_v2.py available in [insapathogenomics/ComparingPartitions](https://github.com/insapathogenomics/ComparingPartitions) and get_best_correspondence.py available in [insapathogenomics/WGS_cluster_congruence](https://github.com/insapathogenomics/WGS_cluster_congruence), developed by Mixão et. al.,2025. 
  - The comparing_partitions_v2.py script was developed based on the previous comparing_partitions.py methodology [jacarrico/ComparingPartitions](https://github.com/jacarrico/ComparingPartitions). This new version introduces two options: single methods and between methods.
    
      - The between_methods option:
          - The aim is to compare 2 methods (e.g., minimum spanning tree and hierarchical cluster single linkage) and check how well the clustering results generated by different analysis methods coincide/are congruent with each other at a given threshold. This analysis is important for assessing the consistency and reliability of typing methods.
          - The comparison of methods is determined by calculating the congruence score.  
          - The congruence score ranges from 0 (no congruence) to 3 (total congruence), because it results from the sum of AWCa→b with AWCb→a and Adjust Rand. The file final_score.tsv return this information.
          - The AWC is the probability that two samples that group together by one method at a given threshold will also be together in the other method. This is done for all thresholds. 
          - The AR is a measure of agreement between typing methods and varies between 0 and 1.
            
     - The single method option:
          - It's a range of thresholds that provides consistent clustering results.
          - It is calculated using the neighbourhood Adjust Wallace Coefficient (nAWC) with the script comparing partitions v2.py, which returns this information with the file metrics.tsv.
          - The script's stability option evaluates the number and composition of clusters obtained by a given method (MST, for example).
          - Subsequent comparisons are made between consecutive thresholds in a sequence, where each threshold is compared with its predecessor. The nAWC is calculated for each pair of thresholds (‘n + 1’ → ‘n’).
          - A region is stable when a sequence of at least 5 thresholds has comparison pairs with an AWC above 0.99. nAWC values range up to 1.
          - Useful for assessing the population structure of the species.
     - For both options can be applied a threshold argument (--threshold) on the partition file, which filters the matrice for these range thresholds.

- Get_best_correspondence.py available in [insapathogenomics/WGS_cluster_congruence](https://github.com/insapathogenomics/WGS_cluster_congruence) was used to make comparisons between two different pipelines in order to assess what is the threshold with the most similar clustering results. Only comparisons that produced CS ≥ 2.85 were considered possible ‘corresponding points’.
    - The score used in this script can be changed with the score argument (--score) until a range between 0-3.
        
  ## Outbreak
  - This module is only available when the user introduces the argument threshold_outbrek (-to) and both inputs are ReporTree folders containing the cluster composition file.
  - Evaltree calls another script named stats_outbreak_analysis.py available in [insapathogenomics/WGS_cluster_congruence](https://github.com/insapathogenomics/WGS_cluster_congruence) to make the outbreaks cluster composition analysis.
  - This script determines the number of clusters identified in a pipeline at a given threshold that could be detected with the same composition by another pipeline at a similar or even higher threshold.
  - If the user wants to re-analyse the threshold outbreak with other values, they can use the repeat threshold outbreaks option (-rto) to repeat just this analysis, returning a new HTML report.


### Installation with conda

```bash

git clone https://github.com/jg-pereira/CENTAUR/EvalTree.git
cd CENTAUR/EvalTree
mkdir scripts
cd scripts
```
Warning: Make sure you are inside the EvalTree/scripts directory before running the commands below.

```bash
git clone https://github.com/insapathogenomics/WGS_cluster_congruence
git clone https://github.com/insapathogenomics/ComparingPartitions
conda env create --name evaltree --file=evaltree_env.yml
```

Run pytest to check that your installation was successful in the EvalTree directory:

```bash
pytest
```

Run EvalTree:
```bash
python evaltree.py -h
```

### Dependencies:
## Usage

  ```bash
  -h, --help            Show this help message and exit
  -i1 INPUT1, --input1 INPUT1
                        [MANDATORY] Specifies the first input type (folder or file), requiring the full path. The folder must contain the partition matrix file with clustering data, and is highly recommended
                        to be a Reportree output folder containing all relevant analysis files. Alternatively, the file can be a traditional sequence-type matrix or a partition matrix. Using either of these
                        input types enables the analysis.
  -i2 INPUT2, --input2 INPUT2
                        [OPTIONAL] Specifies the second input type (folder or file), requiring the full path. The folder must contain the partition matrix file with clustering data, and is highly recommended
                        to be a Reportree output folder containing all relevant analysis files. Alternatively, the file can be a traditional sequence-type matrix or a partition matrix. Using either of these
                        input types enables the analysis.
  -o OUTPUT, --output OUTPUT
                        [MANDATORY] Specifies the output directory for storing all analysis results.
  -list {partitions_summary,sample_of_interest}, --list {partitions_summary,sample_of_interest}
                        [OPTIONAL] Specify the names of the columns present in the *_partition_summary.tsv or *_SAMPLE_OF_INTERES_partition_summary.tsv file.

Congruence

  -s SCORE, --score SCORE
                        [OPTIONAL] Define a minimum score to consider two partitions (one from each pipeline) as corresponding. The score accepts values between 0 and 3. Partition - It refer to the number of
                        identical clusters that exist at the same threshold.
  -t THRESHOLD, --threshold THRESHOLD
                        [OPTIONAL] Defines an integer range to select or filter threshold columns from the partition matrix file. A filtered partition matrix, containing only the selected columns, will be
                        created and used for subsequent analysis. Ranges are specified using a hyphen to separate the minimum and maximum values (e.g., 10-20). If this option is not set, the script will
                        perform clustering for all possible thresholds in the range 0 to the maximum threshold.
 
 Visualising clustering

  -ps {partitions_summary,sample_of_interest}, --plots_summary {partitions_summary,sample_of_interest}
                        [OPTIONAL] Specifies the type of cluster characterization file (*_partition_summary.tsv or *_SAMPLES_OF_INTEREST_partitions_summary.tsv), both of which are expected to be located
                        within a Reportree results folder. Using the partition_summary option, the largest clusters present in the file will be characterized. Alternatively, the samples_of_interest option
                        will characterize all clusters, including those resulting from the addition of new samples (kept increase, new, new (increase), new (merge_increase), new (split_increase), new
                        (split_merge_increase)).
  -n N_CLUSTER, --n_cluster N_CLUSTER
                        [OPTIONAL] Specifies the number of top clusters to be displayed from the *_partition_summary.tsv file, which must be located within a Reportree results folder. This argument is not
                        applicable when using the samples_of_interest option.
  -cp COLUMNS_PLOTS, --columns_plots COLUMNS_PLOTS
                        [OPTIONAL] Name(s) of the column(s) to process the characterization of the clustering data in the selected file (specified by the plots_summary argument). For multiple column names,
                        indicate them separated by commas without spaces (e.g., column1,column2).
  -pt PLOTS_THRESHOLD, --plots_threshold PLOTS_THRESHOLD
                        [OPTIONAL] Identifies the integer threshold(s) to be applied to the file specified by the plots_summary argument. For multiple thresholds, indicate them separated by commas without
                        spaces (e.g., X,Y,Z). This generates a pie chart showing the clustering data for the specified threshold(s), according to the columns_plot argument.
  -pcn PLOTS_CATEGORY_NUMBER, --plots_category_number PLOTS_CATEGORY_NUMBER
                        [OPTIONAL] Determines the number of plot categories in the *_partition_summary.tsv or *_sample_of_interest_partition_summary.tsv file that are intended to be collapsed into the Other
                        category for visualization in the cluster plots. When there are more than 5 slices (default), they will be combined into one category named Other
  -pcp PLOTS_CATEGORY_PERCENTAGE, --plots_category_percentage PLOTS_CATEGORY_PERCENTAGE
                        [OPTIONAL] Determines the percentage of plot categories in the *_partition_summary.tsv or *_sample_of_interest_partition_summary.tsv file that are intended to be collapse into the
                        Othercategory for visualization in the cluster plots. Slices plots with a lower percentage than the entered plots_category_percentage will be combined into one category named Others
Outbreak 

-to THRESHOLD_OUTBREAK, --threshold_outbreak THRESHOLD_OUTBREAK
                        [OPTIONAL] Determine the number of clusters identified in one pipeline at a given threshold that will exist with the same composition in another pipeline at the same or a higher
                        threshold. Full attention, this argument has its own structure: two threshold (strings-methods) and the type of comparison is either equal (defined by , ) or lower_equal (defined by
                        <= ) Threshold1: Threshold at which the genetic clusters must be identified for the pipeline of interest. Threshold2: Threshold at which the genetic clusters must be searched in the
                        other pipelines. Comparison (equal or lower equal): - equal: Used to assess whether a cluster is detected at a given threshold by another pipeline. Use a comma , to separate
                        threshold1,threshold2. Example of expression: MST-7x1.0,MST-7x1.0. - lower_equal: Used to assess whether a cluster is detected up to a given threshold in another pipeline. Use <=
                        between threshold1<=threshold2. Example of expression: MST-7x1.0<=MST-9x1.0. For multiple pair of threshold values, use ; as a separator. Example of expression:
                        "MST-7x1.0,MST-7x1.0;<=MST-7x1.0,MST-10x1.0" represents two pair of threshold values.

  -rto, --repeat_threshold_outbreak
                        [OPTIONAL] This argument can only be used after of a previous analysis of threshold_outbreak.

```

### Simple EvalTree command line with  two input ReporTree folders 
```bash
evaltree.py -i1 input1 -i2 input2 -o output -ps partitions_summary -pt MST-7x1.0 -cp name_column -to "MST-7x1.0,MST-7x1.0;<=MST-7x1.0,MST-9x1.0"
```
### Simple EvalTree command line with  two files 
```bash
evaltree.py -i1 X_partitions.tsv -i2 Y_partitions.tsv -o output
```
  
## Citation

If you run EvalTree, please cite the publication:

Mixão, V., Pinto, M., Brendebach, H., Sobral, D., Santos, J. D., Radomski, N., Uldall, A. S. M., Bomba, A., Pietsch, M., Bucciacchio, A., de Ruvo, A., Castelli, P., Iwan, E., Simon, S., Coipan, C. E., Linde, J., Petrovska, L., Kaas, R. S., Joensen, K. G., Nielsen, S. H., Kiil, K., Lagesen, K., Di Pasquale, A., Gomes, J. P., Deneke, C., Tausch, S. H., & Borges, V. (2025). Multi-country and intersectoral assessment of cluster congruence between pipelines for genomics surveillance of foodborne pathogens. Nature Communications, 16, Article 3961. https://doi.org/10.1038/s41467-025-59246-8

EvalTree relies on the work of other developers. So you must cite:

1) ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true 


## Funding
This work was supported by funding from the xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
