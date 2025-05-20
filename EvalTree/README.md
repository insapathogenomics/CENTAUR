# EvalTree
_EvalTree.py_ is a tool developed in the frame of the CENTAUR project (supported by the European ISIDORe initiative) that, among other features, **compares cluster composition at all possible resolution levels between two typing solutions, and their performance in identifying similar outbreak signals**.

## Possible applications
- Comparison of the WGS-based pipelines in place at two different laboratories, facilitating inter-laboratory communication and cooperation and promoting large-scale External Quality Assessments
- Comparison of two different versions of the same WGS-based pipeline, providing confidence in the implementation of pipeline updates, promoting the long.term sustainability of pipelines
- Comparison of different WGS-based pipeline components - useful during cg/wgMLST schema creation
- Comparison of a WGS-based pipeline and traditional typing classification, providing information about backwards comparability and, thus, supporting the technological transition to WGS

## Implementation

_EvalTree.py_ is a command-line tool, implemented in Python 3.8, able to provide user-friendly reports of cluster congruence assessments between two typing pipelines. This tool orchestrates the methodology developed by [Mixão et al. 2025](https://www.nature.com/articles/s41467-025-59246-8).

## Input
_EvalTree.py_ is designed primarily to make comparisons between the typing classification obtained by two different methods (e.g. cgMLST, sequence type, serotype, clonal complex, etc.). For each pipeline/typing system under comparison, one of the following input types can be provided:
- Folder with the [ReporTree](https://github.com/insapathogenomics/ReporTree) run
- TSV file with the clustering information at one (e.g. sequence type or serotype) or multiple levels (e.g. cgMLST clustering)
  
_EvalTree.py_ is able to accept typing information from different origins as input, as far as a TSV file in which rows correspond to samples and column(s) correspond to the typing classification method/level is provided for each method under comparison. Still, if you used [ReporTree](https://github.com/insapathogenomics/ReporTree) to obtain your cgMLST clusters, you can provide your output ReporTree folder as input to _EvalTree.py_ and this tool will use ReporTree output files for additional analyses:
- _partitions.tsv_: The main input. It contains clustering information at all possible distance thresholds, and will be used for the cluster congruence assessment with the other method.
- _clusterComposition.tsv_: File that summarizes the clusters and singletons identified at each threshold level. This file can be used for an in-depth analysis at potential outbreak-level (as specified by the user).
- _partitions_summary.tsv_: Characterization of genetic clusters with the available metadata. This file can be used in order to also provide a graphical visualization of the composition of the clusters, according to user's specifications.
- _SAMPLES_OF_INTEREST_partitions_summary.tsv_ - similar to _partitions_summary.tsv_ but exclusively for the samples of interest. If the user specifies that the graphical visualization should only be performed for clusters with "samples of interst", this file is used by _EvalTree.py_ instead of the _partitions_summary.tsv_.
- _stableRegions.tsv_: File that reports the stability regions observed in a pipeline. This file is used by _EvalTree.py_for a graphical representation of these regions in the pipeline characterization.

## Output
The main output file of _EvalTree.py_ is an **HTML report** containing interactive Plotly visualisations of all the results generated from an analysis. This report is divided into two main sections:    
  
**1. Pipeline characterization** (this section is provided for each typing method under comparison)  
- Displays the name of the pipeline
- Shows the number of samples and thresholds used
- Includes a scarlet plot illustrating the number of partitions (clusters) per threshold/classification system

If a ReporTree folder is provided as input, this section may also show:
- A graphical representation of the stability regions identified in each pipeline
- Pie charts with the characterization of the genetic clusters observed at a (or multiple) user-defined threshold(s) according to the metadata variable indicated by the user

**2. Congruence analysis**  
- Heatmap indication the congruence score obtained from the pairwise comparisons at all possible threshold levels
- A scatterplot indicating the inter-pipeline corresponding points and their respective trend line 

If a ReporTree folder is provided as input, this section may also show:
- Graphical representation of the ability of the two pipelines to detect the same outbreak signals

**Source files**   
Besides the HTML report, this tool provides a wide variety of source files that are generated throughout the analysis. Examples of these files are:
- _All_correspondence.tsv_: It reports the inter-pipeline corresponding thresholds between two pipelines
  
## Description of the main modules of EvalTree

  ## Characterization
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


## Installation with conda

```bash
git clone https://github.com/jg-pereira/CENTAUR.git
cd CENTAUR/EvalTree/
mkdir scripts
cd scripts
git clone https://github.com/insapathogenomics/WGS_cluster_congruence
git clone https://github.com/insapathogenomics/ComparingPartitions
```

To create the conda environment:
```bash
cd ..
conda env create --name EvalTree --file=EvalTree_env.yml
```

Run pytest to check that your installation was successful in the EvalTree directory:
```bash
pytest
```

Run EvalTree:
```bash
python EvalTree.py -h
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
                        threshold1,threshold2. Example of expression: "MST-7x1.0,MST-7x1.0". - lower_equal: Used to assess whether a cluster is detected up to a given threshold in another pipeline. Use <=
                        between threshold1<=threshold2. Example of expression: MST-7x1.0<=MST-9x1.0. For multiple pair of threshold values, use ; as a separator. Example of expression:
                        "MST-7x1.0,MST-7x1.0;<=MST-7x1.0,MST-10x1.0" represents two pair of threshold values.

  -rto, --repeat_threshold_outbreak
                        [OPTIONAL] This argument can only be used after of a previous analysis of threshold_outbreak.

```

### A simple EvalTree command line example using two input ReporTree folders 
```bash
python EvalTree.py -i1 input1 -i2 input2 -o output -ps partitions_summary -pt MST-7x1.0 -cp name_column -to "MST-7x1.0,MST-7x1.0;<=MST-7x1.0,MST-9x1.0"
```
### A simple EvalTree command line example using two files 
```bash
python EvalTree.py -i1 X_partitions.tsv -i2 Y_partitions.tsv -o output
```
  
## Citation

If you run EvalTree, please cite the publication:

[Mixão V et al. (2025). Multi-country and intersectoral assessment of cluster congruence between pipelines for genomics surveillance of foodborne pathogens. Nature Communications, 16, Article 3961. https://doi.org/10.1038/s41467-025-59246-8](https://doi.org/10.1038/s41467-025-59246-8)

EvalTree relies on the work of other developers. So you must also cite:
- ReporTree: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1
- ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true 


## Funding
This work was supported by the ISIDORe project (funding from the European Union’s Horizon Europe Research & Innovation Programme, Grant Agreement no. 101046133) and by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001.
