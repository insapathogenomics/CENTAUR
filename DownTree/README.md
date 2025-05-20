# DownTree

_DownTree.py_ was designed to select a subset of assemblies from a large and diverse dataset that is still representative of the global dataset diversity. For this, _DownTree.py_ takes into consideration information about previously determined cgMLST clusters (or any other grouping variable of interest) and the assemblies' metrics.

## Possible applications
- Selection a subset of high-quality genomes representative of the species diversity for a cgMLST schema creation
- Selection of the high-quality genomes that represent a given populational group (e.g. Sequence Type, Serotype, etc.)
- Selection of the high-quality genome that represents a given outbreak cluster for posterior in-depth SNP-based analysis

## Input
- TSV file with clustering information (e.g. cgMLST clusters at all threshold levels, traditional typing classification, _k_-mer based clustering, etc.) from a large dataset of isolates
- TSV file with assembly metrics

_Note: if a previous cgMLST schema or a multi-sequence alignment are available, we strongly recomend that you use the partitions.tsv file of [ReporTree](https://github.com/insapathogenomics/ReporTree), which provides clustering information at all possible threshold levels for a dataset._  

## Output
- TSV file with information for each grouping variable (cluster) of the selected representative isolate(s)
- TSV file with assembly metrics including additional columns informing if a sample was selected as representative and if it belongs to a cluster/group that is represented in the final dataset.

## Installation

```
git clone https://github.com/vmixao/CENTAUR.git
cd CENTAUR/DownTree/
conda env create --name downtree --file=downtree.yml
python DownTree.py -h
```

## Usage
```
options:
  -h, --help            show this help message and exit

Version:
  -v, --version         Print version and exit

Main options:
  DownTree.py main options

  -c, --clusters CLUSTERS
                        [MANDATORY] File with the clustering data in *.tsv format. The first column must correspond to the sample name, and the other(s) to the clustering classification.
                        To take the most advantage of this tool, we recommend that you use the *_partitions.tsv file of ReporTree. Still, this tool will also work with other
                        classifications, such as, for example, MLST, Serotype, etc.
  -m, --metadata METADATA
                        [MANDATORY] File with the assembly quality information in *.tsv format.
  -g, --grouping_column GROUPING_COLUMN
                        [OPTIONAL] Column of the '--clusters' file that must be used to select samples based on genetic diversity. If this column is not specified, this tool will search
                        for the column that provides the number of isolates that is closer to the number specified in '--number_assemblies'.
  -n, --number_assemblies NUMBER_ASSEMBLIES
                        Approximate number of representative assemblies to be selected from the whole dataset. If '--grouping_column' is specified, this argument will be ignored.
                        (default: 5000)
  -s, --samples_per_group SAMPLES_PER_GROUP
                        Number of representative assemblies to be selected per cluster. (default: 3)
  -col, --column_name COLUMN_NAME
                        [MANDATORY] Name of the column of the '--metadata' file that must be used to infer the assembly quality (data in this column must be numeric).
  -qt, --quality-threshold QUALITY_THRESHOLD
                        [MANDATORY] Maximum/minimum (depending on the '--order' argument) number in the column specified with '--column_name' that corresponds to the threshold between an
                        elegible and a non-elegible sample for selection.
  -order, --order ORDER
                        This argument can take one of two options: 'low' - indicating that lower values represent higher quality (e.g. contigs); or 'high' - indicating that higher values
                        represent higher quality (e.g. N50). (default: low)
  -t, --tag TAG         Tag for output file names. (default: DownTree)
  --evaluation          Performs just the evaluation of the dataset and the clustering data. If you set this argument, no dataset will be generated.

```

## _DownTree.py_ running options
_DownTree.py_ can be in one of two ways:
1. The user determines a '--grouping_column' that must be used to assess genetic diversity. _DownTree.py_ then selects the _s_ genome assemblies representative of each of the groups of the column.
2. The user determined the target '--number_assemblies' and _DownTree.py' identifies the grouping column providing the closest number of representatives taking into consideration the defined _s_ number of genome assemblies per group and the assemblies passing the '--quality-threshold'.

_Note: As only assemblies passing the '--quality-threshold' can be selected, it is possible that some groups do not have any representative sequence._

#### What is the '--evaluation' argument?
This argument analyzes the dataset and gives you an overview of i) how many assemblies will be selected; ii) what is the percentage of samples represented in the subdataset; and iii) what is the percentage of clusters of the grouping column determined by the user (option 1) or selected py the script (option 2) that is represented in the subdataset.
This allows the user of having a previous idea of the final results that will be obtained without the need of waiting for the script to run, and adjust the parameters accordingly.

#### How does _DownTree.py_ select representative isolates of groups/clusters with more than _s_ samples?
If the user provides as input a partitions.tsv table obtained with ReporTree, _DownTree.py_ will analyze all columns providing higher resolution than the '--grouping_column' and prioritize the selection of the more divergent samples.

![image](https://github.com/user-attachments/assets/af0351ef-ea19-458b-87e2-8be234461a1f)


## Example of command line
### Option 1
```
 python DownTree.py -c ReporTree_partitions.tsv -m assembly_metrics.tsv -g "MST-4x1.0" -s 3 -col contigs -qt 100 -order low -t DownTree
```

### Option 2
```
python DownTree.py -c ReporTree_partitions.tsv -m assembly_metrics.tsv -n 5000 -s 3 -col contigs -qt 100 -order low -t DownTree
```

## Funding
This work was supported by the ISIDORe project (funding from the European Union’s Horizon Europe Research & Innovation Programme, Grant Agreement no. 101046133) and by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001.
