# PanTree

_PanTree.py_ was designed to evaluate loci prevalence and allele diversity, and map their presence/absence in the “species” tree topology, providing evidence for an informed decision on the inclusion of accessory loci during the creation of a novel wgMLST schema.

## Possible applications
- Evaluation of loci prevalence and distribution during the creation of a novel cg/wgMLST schema

![image](https://github.com/user-attachments/assets/8815037e-9f3a-4de8-a98e-73e83e9c872a)

## Input
- Allelic matrix (TSV) of the cg/wgMLST schema under development
- Species tree (NWK)
- Clustering information based on a previously used typing method

_Note: if a previous cgMLST schema or a multi-sequence alignment are available, we strongly recomend that you use the tree and the partitions.tsv file of [ReporTree](https://github.com/insapathogenomics/ReporTree)_  

## Output
- TSV file with information about the prevalence of each cg/wgMLST loci in the dataset, their allelic diversity and distribution across clusters identified at level X

## Usage


