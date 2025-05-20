# PanTree

_PanTree.py_ was designed to evaluate loci prevalence and allele diversity, and map their presence/absence in the “species” tree topology, providing evidence for an informed decision on the inclusion of accessory loci during the creation of a novel wgMLST schema.

## Possible applications
- Evaluation of loci prevalence and distribution during the creation of a novel cg/wgMLST schema

![image](https://github.com/user-attachments/assets/5522281a-6adf-4154-9236-e4e18f6f5ec6)


## Input
- Allelic matrix (TSV) of the cg/wgMLST schema under development
- Species tree (NWK)
- Clustering information based on a previously used typing method

_Note: if a previous cgMLST schema or a multi-sequence alignment are available, we strongly recomend that you use the tree and the partitions.tsv file of [ReporTree](https://github.com/insapathogenomics/ReporTree)_  

## Output
- TSV file with information about the prevalence of each cg/wgMLST loci in the dataset, their allelic diversity and distribution across clusters identified at level X

## Installation

```
git clone https://github.com/vmixao/CENTAUR.git
cd CENTAUR/PanTree/
conda env create --name reportree --file=pantree.yml
python PanTree.py -h
```

## Usage
```
options:
  -h, --help            show this help message and exit
  -a ALLELE, --allele ALLELE
                        [MANDATORY] Path to the input allele matrix (TSV). Missing data must be represented by 0.
  -p PARTITIONS, --partitions PARTITIONS
                        [MANDATORY] Path to the input partitions table (TSV). We recomend using the partitions table obtained with ReporTree.
  -tree TREE, --tree TREE
                        [MANDATORY] Path to the phylogenetic tree (NWK). It can be a single-linkage tree obtained with ReporTree.
  -s COLUMN_SELECT, --column2select COLUMN_SELECT
                        [MANDATORY] Comma-separated list of the thresholds (column names) of the partitions table that will be used to flag loci dispersed throughout the tree (e.g. --column2select single-1x1.0,single-10x1.0,single-100x1.0).
  -cp CLUSTER_PREVALENCE, --cluster-prevalence CLUSTER_PREVALENCE
                        [MANDATORY] Percentage of clusters/groups were a given locus is present, so it can be considered as elegible for the wgMLST.
  -o OUTPUT, --output OUTPUT
                        [MANDATORY] Tag for output files.
  --order {yes,no}      Order loci by the percentage of presence - descending [default: no].
  --threshold THR       Only analyze loci with a prevalence lower than this threshold (value between 0 and 1).
```

## Example of command line
```
python PanTree.py -a alleles.tsv -p ReporTree_partitions.tsv -tree Reportree_tree.nwk -s "single-1x1.0,single-10x1.0,single-100x1.0" -cp 0.1 -o PanTree --order yes --threshold 0.98
```

## Funding
This work was supported by the ISIDORe project (funding from the European Union’s Horizon Europe Research & Innovation Programme, Grant Agreement no. 101046133) and by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001.
