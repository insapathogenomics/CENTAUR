# ACDTree

This repository contains a novel set of tools (currently under development in the frame of CENTAUR project, supported by the European ISIDORe initiative) that address key bioinformatics challenges faced throughout the life-cycle of a cg/wgMLST pipeline. Each of these tools will have soon their own GitHub page. Stay tuned!

## Motivation
Routine genomics-based surveillance of bacterial pathogens is increasingly relying on cg/wgMLST approaches. Currently, there is a need for bioinformatics solutions that facilitate the development, deployment and application of new cg/wgMLST workflows, while ensuring their long-term sustainability and inter-laboratory comparability in the Era of global surveillance. In this repository, we present a **novel set of tools that address key bioinformatics challenges faced throughout the life-cycle of a cg/wgMLST pipeline**, illustrating their added-values across the different stages of a **brand-new implementation of routine genomics surveillance for a bacterial pathogen X**.

<p align="center">
  <img src="https://github.com/user-attachments/assets/693c9fa2-8ed0-4bf8-b273-d65cefcadfd3" alt="pipeline_life_cycle" />
</p>

This toolkit is composed of independent tools, whose usage can be potentiated if combined:

## _DownTree_
_DownTree.py_ was designed to select a subset of assemblies from a large and diverse dataset that is still representative of the global dataset diversity. For this, _DownTree.py_ takes into consideration information about previously determined cgMLST clusters (or any other grouping variable of interest) and the assemblies' metrics.

#### Possible applications
- Selection a subset of high-quality genomes representative of the species diversity for a cgMLST schema creation
- Selection of the high-quality genomes that represent a given populational group (e.g. Sequence Type, Serotype, etc.)
- Selection of the high-quality genome that represents a given outbreak cluster for posterior in-depth SNP-based analysis

#### Input
- TSV file with clustering information (e.g. cgMLST clusters at all threshold levels, traditional typing classification, _k_-mer based clustering, etc.) from a large dataset of isolates
- TSV file with assembly metrics

_Note: if a previous cgMLST schema or a multi-sequence alignment are available, we strongly recomend that you use the partitions.tsv file of [ReporTree](https://github.com/insapathogenomics/ReporTree), which provides clustering information at all possible threshold levels for a dataset._  

#### Output
- TSV file with information for each grouping variable (cluster) of the selected representative isolate(s)
- TSV file with assembly metrics including additional columns informing if a sample was selected as representative and if it belongs to a cluster/group that is represented in the final dataset.

You can find more details about the usage of this tool in [DownTree's page](https://github.com/insapathogenomics/CENTAUR/blob/main/DownTree/README.md).

## _PanTree_
_PanTree.py_ was designed to evaluate loci prevalence and allele diversity, and map their presence/absence in the “species” tree topology, providing evidence for an informed decision on the inclusion of accessory loci during the creation of a novel wgMLST schema.

#### Possible applications
- Evaluation of loci prevalence and distribution during the creation of a novel cg/wgMLST schema

#### Input
- Allelic matrix (TSV) of the cg/wgMLST schema under development
- Species tree (NWK)
- Clustering information based on a previously used typing method

_Note: if a previous cgMLST schema or a multi-sequence alignment are available, we strongly recomend that you use the tree and the partitions.tsv file of [ReporTree](https://github.com/insapathogenomics/ReporTree), which provides clustering information at all possible threshold levels for a dataset._  _  

#### Output
- TSV file with information about the prevalence of each cg/wgMLST loci in the dataset, their allelic diversity and distribution across clusters identified at level X
- Graphical representation of loci presence/absence throughout a species tree 

You can find more details about the usage of this tool in [PanTree's page](https://github.com/insapathogenomics/CENTAUR/blob/main/PanTree/README.md).

## _EvalTree_
_EvalTree.py_ was designed to automate a [previously developed methodology](https://github.com/insapathogenomics/WGS_cluster_congruence) for evaluation of congruence between two clustering analysis ([Mixão et al. 2025](https://www.nature.com/articles/s41467-025-59246-8)).

#### Possible applications
- Comparison of the WGS-based pipelines in place at two different laboratories, facilitating inter-laboratory communication and cooperation and promoting large-scale External Quality Assessments
- Comparison of two different versions of the same WGS-based pipeline, providing confidence in the implementation of pipeline updates
- Comparison of different WGS-based pipeline components - useful during schema creation
- Comparison of a WGS-based pipeline and traditional typing classification, providing information about backwards comparability and, thus, supporting the technological transition to WGS  

#### Input
For each pipeline/typing system under comparison, one of the following:
- Folder with the [ReporTree](https://github.com/insapathogenomics/ReporTree) run
- TSV file with the clustering information at one or multiple levels

#### Output
- HTML report of the congruence between cluster composition at all possible resolution levels between two typing solutions, and their performance in identifying similar outbreak signals

You can find more details about the usage of this tool in [EvalTree's page](https://github.com/insapathogenomics/CENTAUR/blob/main/EvalTree/README.md).

## _ReporTree_

_ReporTree_ is a flexible pipeline that facilitates the [detection of genetic clusters and their linkage to epidemiological data (Mixão et al. 2023)](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1).

#### Possible applications
- Monitoring genetic clusters and perform their characterization in routine surveillance
- Determine genetic clusters at all possible threshold levels using a cg/wgMLST or a SNP-based approach
- Perform a dynamic zoom-in analysis of an outbreak cluster, extending the static cgMLST with accessory loci that are present in this cluster

**NOTE: [ReporTree](https://github.com/insapathogenomics/ReporTree) outputs can be used as input for any of the above-mentioned tools (_DownTree.py_, _PanTree.py_ and _EvalTree.py_)!**
  
You can find more details about the usage of this tool in [ReporTree's page](https://github.com/insapathogenomics/ReporTree).

## Funding
This work was supported by the ISIDORe project (funding from the European Union’s Horizon Europe Research & Innovation Programme, Grant Agreement no. 101046133) and by national funds through FCT - Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001.
