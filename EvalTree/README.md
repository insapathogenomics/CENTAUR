# EvalTree
_EvalTree.py_ is a tool developed in the frame of the CENTAUR project (supported by the European ISIDORe initiative) that, among other features, **compares cluster composition at all possible resolution levels between two typing solutions, and their performance in identifying similar outbreak signals**.

## Possible applications
- Comparison of the WGS-based pipelines in place at two different laboratories, facilitating inter-laboratory communication and cooperation and promoting large-scale External Quality Assessments
- Comparison of two different versions of the same WGS-based pipeline, providing confidence in the implementation of pipeline updates, promoting the long-term sustainability of pipelines
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
  
# Script orchestration
  - _EvalTree.py_ automatically orchestrates the execution of scripts, handling inputs and outputs to perform the following analysis:
  
  ### Congruence analysis
  - It uses **comparing_partitions_v2.py** script, available in [insapathogenomics/ComparingPartitions](https://github.com/insapathogenomics/ComparingPartitions)
  - The script was developed based on the previous comparing_partitions.py methodology [jacarrico/ComparingPartitions](https://github.com/jacarrico/ComparingPartitions)
  - The analysis method is **Between_methods**:
    
      -  For each pairwise pipeline comparison, a **congruence score (CS)** is computed, ranging from 0 (no congruence) to 3 (total congruence), to assess the consistency between pipelines.
      -  This score is calculated as follows:
          - Adjusted Wallace coefficient (AWC): Probability that two samples that cluster together using one method (at a given threshold level) also cluster together with another one (at a given threshold level). This calculation is performed for all thresholds in both directions (method A → method B and method B → method A)
          - Adjusted Rand (AR) coefficient: Measures the overall agreement between the typing methods
          - The final CS is defined as:
              CS = AWC A→B + AWC B→A + AR

  ### Identification of pipeline stability regions     
  - It also uses the **comparing_partitions_v2.py** script
  - Stability regions are defined as threshold ranges where clustering results are consistent 
  -  The analysis method is **Stability**:
      - neighborhood Adjusted Wallace coefficient (nAWC): Subsequent comparisons are made between consecutive thresholds (n + 1 → n) in a sequence, where each threshold is compared with its previous one  

 ### Identification of corresponding points
  - It uses **get_best_correspondence.py** script available in [insapathogenomics/WGS_cluster_congruence](https://github.com/insapathogenomics/WGS_cluster_congruence)
  - This analysis identifies, for each pipeline comparison, the threshold that provides the most similar clustering results in the other pipeline referred to as the best corresponding point. This is determined based on the highest CS scores
  -  Only comparisons yielding CS ≥ 2.85, which ensures a score ≥0.95 for each CS metric component, were considered as possible corresponding points 
      
 ### Outbreak
  - If outbreak analysis was specified, _EvalTree.py_ calls stats_outbreak_analysis.py script available in [insapathogenomics/WGS_cluster_congruence](https://github.com/insapathogenomics/WGS_cluster_congruence)
  - It determines the number of clusters identified in a pipeline at a given threshold that could be detected with the same composition by another pipeline at a similar or even higher threshold

# Installation and dependencies

### Dependencies:
- pandas
- plotly
- pytest
- kaleido
- scikit-learn
- scipy
- statsmodels

### Installation with conda

#### 1. Manual installation from GitHub repository

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

#### 2. Install via Conda

```bash
conda install insapathogenomics::evaltree
conda install -c joana.gomes evaltree
```

### Installation with PyPi
Run:
 ```bash
pip install evaltree
```
## Usage

  ```bash
  -h, --help            show this help message and exit
  -v, --version         [OPTIONAL] Specify the version number of EvalTree.
  -i1 INPUT1, --input1 INPUT1
                        [MANDATORY] Specify the first input type (folder or file) and requires the full path. If a
                        folder is provided, it must contain the partition table with clustering data at all possible
                        threshold levels. Alternatively, individual clustering table files from traditional methods or
                        WGS-based methods can be supplied.
  -i2 INPUT2, --input2 INPUT2
                        [OPTIONAL] Specify the second input type (folder or file) and requires the full path. If a
                        folder is provided, it must contain the partition table with clustering data at all possible
                        threshold levels. Alternatively, individual clustering table files from traditional methods or
                        WGS-based methods can be supplied.
  -o OUTPUT, --output OUTPUT
                        [OPTIONAL] Specify the output directory for storing all analyses results. If no folder is
                        provided, the program will automatically create one based on the prefix of the files.
  -s SCORE, --score SCORE
                        [OPTIONAL] Define a minimum score to consider two partitions (one from each pipeline) as
                        corresponding. The score accepts values between 0 (low congruence) and 3 (high congruence).
  -t THRESHOLD, --threshold THRESHOLD
                        [OPTIONAL] Define an integer threshold range to select or filter threshold columns from the
                        partition table file. A filtered partition table, containing only the selected columns, will
                        be created and used for subsequent analysis. Threshold ranges are specified using a hyphen to
                        separate the minimum and maximum values (e.g., 10-20). If this option is not set, the script
                        will perform clustering for all possible thresholds in the range 0 to the maximum threshold.
  -ps {partitions_summary,sample_of_interest}, --plots_summary {partitions_summary,sample_of_interest}
                        [OPTIONAL] Specify the type of cluster characterization file (partitions_summary.tsv or
                        SAMPLES_OF_INTEREST_partitions_summary.tsv. Both files are expected to be located within a
                        ReporTree folder. Using the partition_summary option, the largest clusters present in the file
                        will be characterized. Alternatively, the samples_of_interest option will characterize all
                        clusters, including those resulting from the addition of new samples (kept increase, new, new
                        (increase), new (merge_increase), new (split_increase), new (split_merge_increase)).
  -n N_CLUSTER, --n_cluster N_CLUSTER
                        [OPTIONAL] Specify the number of top clusters to be displayed from the partitions_summary.tsv
                        file, which must be located within a ReporTree results folder. This argument is not applicable
                        when using the sample_of_interest option.
  -cp COLUMNS_PLOTS, --columns_plots COLUMNS_PLOTS
                        [OPTIONAL] Name(s) of the column(s) to process the characterization of the clustering data in
                        the selected file (specified by the plots_summary argument). For multiple column names
                        (categories), indicate them separated by commas without spaces (e.g., column1,column2).
  -pt PLOTS_THRESHOLD, --plots_threshold PLOTS_THRESHOLD
                        [OPTIONAL] Identify the integer threshold(s) to be applied to the file specified by the
                        plots_summary argument. For multiple thresholds, indicate them separated by commas without
                        spaces (e.g., X,Y,Z). This generates a pie chart showing the clustering data for the specified
                        threshold(s), according to the columns_plot argument.
  -pcn PLOTS_CATEGORY_NUMBER, --plots_category_number PLOTS_CATEGORY_NUMBER
                        [OPTIONAL] Determine the number of plot categories in the partitions_summary.tsv or
                        SAMPLES_OF_INTEREST_partitions_summary.tsv file that are intended to be collapsed into the
                        Other category for visualization in the cluster plots. When there are more than 5
                        slices/segments (default), they will be combined into one category named "Other".
  -pcp PLOTS_CATEGORY_PERCENTAGE, --plots_category_percentage PLOTS_CATEGORY_PERCENTAGE
                        [OPTIONAL] Determine the percentage of plot categories in the partitions_summary.tsv or
                        SAMPLES_OF_INTEREST_partitions_summary.tsv file that are intended to be collapse into the
                        "Other" category for visualization in the cluster plots. Slices/segments plots with a lower
                        percentage than the entered plots_category_percentage will be combined into one category named
                        "Other".
  -to THRESHOLD_OUTBREAK, --threshold_outbreak THRESHOLD_OUTBREAK
                        [OPTIONAL] Determine the number of clusters identified in one pipeline (rows) at a given
                        threshold that were also detected with the exact same cluster composition at a (or up to a)
                        given threshold in the other pipeline (columns).
                        This argument has a specific structure: it requires two thresholds and the type of comparison (fixed or dynamic threshold). 
                        Threshold1: Threshold at which the genetic clusters are identified in the pipeline of interest.
                        Threshold2: Threshold at which the genetic clusters are searched in the other pipeline.
                        Comparison (fixed or dynamic):
                        - Fixed: Direct cluster comparison at specific threshold between two pipelines. Use a comma , to separate threshold1 and threshold2. Example: 7,7.
                        - Dynamic: Range of threshold to find genetic clusters with similar composition in the other pipeline. Use <= to indicate a dynamic threshold range. Example of expression: 7,<=9.
                        For multiple threshold pairs, separate them with a semicolon. Example: "7,7;<=7,10" represents two threshold pairs.
  -list {partitions_summary,sample_of_interest}, --list {partitions_summary,sample_of_interest}
                        [OPTIONAL] Specify the names of the columns present in the partitions_summary.tsv or
                        SAMPLES_OF_INTEREST_partitions_summary.tsv file.
  -rto, --repeat_threshold_outbreak
                        [OPTIONAL] This argument can only be used after of a previous analysis of threshold_outbreak
                        and generate a second HTML report.
  -n_stab N_STABILITY, --n_stability N_STABILITY
                        [OPTIONAL] Range of threshold which the cluster composition can be similar.
  -thr_stab THR_STABILITY, --thr_stability THR_STABILITY
                        [OPTIONAL] The neighborhood Adjusted Wallace Coefficient (nAWC) threshold used to determine if
                        a clustering threshold is considered consistent or stable.
   -ttc TRADITIONAL_TYPING_CATEGORY, --traditional_typing_category TRADITIONAL_TYPING_CATEGORY
                        [OPTIONAL] Select the values of the category of traditional typing (e.g., serotype, sequence
                        type), with high number of samples.
```
### Examples of command-line usage with:
  #### Manual Installation using:
  - Two input files
```bash
python EvalTree.py -i1 X_partitions.tsv -i2 Y_partitions.tsv 
```
 - Two ReporTree folders to generate clustering visualizations
```bash
python EvalTree.py -i1 input1 -i2 input2 -o output -ps partitions_summary -pt 7 -cp name_column -to "7,7;<=7,9"
```
  #### Installation via Conda OR PyPi using:
  - Two input files
```bash
evaltree i1 X_partitions.tsv -i2 Y_partitions.tsv -o output
```

## Citation

If you run EvalTree, please cite the publication:

[Mixão V et al. (2025). Multi-country and intersectoral assessment of cluster congruence between pipelines for genomics surveillance of foodborne pathogens. Nature Communications, 16, Article 3961. https://doi.org/10.1038/s41467-025-59246-8](https://doi.org/10.1038/s41467-025-59246-8)

EvalTree relies on the work of other developers. So you must also cite:
- ReporTree: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1
- ComparingPartitions: https://journals.asm.org/doi/10.1128/jcm.02536-05?permanently=true 


## Funding
This work was supported by the ISIDORe project (funding from the European Union’s Horizon Europe Research & Innovation Programme, Grant Agreement no. 101046133) and by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001.
