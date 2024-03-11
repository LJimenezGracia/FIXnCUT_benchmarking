# FIXnCUT_benchmarking

This is the GitHub repository for the **FixNCut: Single-cell genomics through reversible tissue fixation and dissociation** manuscript.


## Abstract

The use of single-cell technologies for clinical applications requires disconnecting sampling from downstream processing steps. Early sample preservation can further increase robustness and reproducibility by avoiding artifacts introduced during specimen handling. We present FixNCut, a methodology for the reversible fixation of tissue followed by dissociation that overcomes current limitations. We applied FixNCut to human and mouse tissues to demonstrate the preservation of RNA integrity, sequencing library complexity, and cellular composition, while diminishing stress-related artifacts. Besides single-cell RNA sequencing, FixNCut is compatible with multiple single-cell and spatial technologies, making it a versatile tool for robust and flexible study designs.


## Code implementation

The repository is organized into the following folder tree, which contains all the necessary files and scripts to perform the detailed tasks and reproduce all our results.

* **01_cellranger_mapping** --> It includes an overview of the project data information, including which samples and 10X libraries were generated. Also, it contains the scripts needed to create a folder directory to perform the sequencing read mapping to the reference genome. Finally, it includes R markdown notebooks to perform a general quality control on the raw sequencing reads from each library, considering different organism and tissue types independently.

* **02_QC** --> All scripts to predict doublets and downsampling reads for each 10X generated library. Also, R markdown notebooks to perform the quality control and the first pre-processing, including data normalization, scaling, dimensionality reduction and integration (when needed), considering different organism and tissue types independently.

* **03_clustering_annotation** --> All R markdown notebooks to decide the optimal clustering resultion, to find differential expressed markers for each clustering, and to assign a biological-relevant identity to each cluster, considering different organism and tissue types independently.

* **04_GEX_analysis** --> Here, you can find all the code used to performed further downstream analysis on the processed data, considering different organism and tissue types independently. In this sense, it includes all the scripts and R markdown notebooks to perform a general quality control on cleaned data, pseudobulk gene correlation for library and cell-type, cell composition analysis, and differential expression analysis (DEA) followed by gene set enrichment analysis (GSEA) & Gene Ontology (GO) using multiple databases. Finally, we also performed an analysis of well-defined gene signature across protocol conditions.

### Package versions

The (most important) packages and versions needed to reproduce the full analysis are listed below:

* CellRanger (v 6.1.1) was used to mapped single-cell RNA-seq reads (10X Genomics) to the reference genome.

*--- in R (v 4.1.0) ---*
* Seurat (v 4.0.0)
* SeuratObject (v 4.0.1)
* SingleCellExperiment (v 1.12.0)
* SummarizedExperiment (v 1.20.0)
* Harmony (v 1.0)
* tidyverse (v 1.3.1)
* dplyr (v 1.0.5)
* ggplot2 (v 3.3.3)
* ggpubr (v 0.4.0)
* dittoSeq (v 1.7.0)
* DT (v 0.17)
* edgeR (v 3.32.1)
* MatrixGenerics (v 1.2.1)
* stringr (v 1.4.0)
* limma (v 3.46.0)
* RColorBrewer (v 1.1.2)

*--- in Python (v 3.8.3) ---*
* IPython (v 7.18.1)
* anndata (v 0.7.6)
* pandas (v 1.1.2)
* matplotlib (v 3.3.1)
* sccoda (v 0.1.2.post1)
* scrublet (v 0.2.1)


## Data accessibility

* The complete raw data (FASTQ files) generated in this study have been submitted to the NCBI Gene Expression Omnibus (GEO) under accession number [GSE229944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229944).
* The count matrices and metadata are deposited at [Zenodo](https://zenodo.org/doi/10.5281/zenodo.7837623).

All data will be available upon publication.


## Code accessibility

You can easily download a copy of all the files contained in this repository by:

* Cloning the git repository using the following command in the terminal:

`git clone https://github.com/LJimenezGracia/FIXnCUT_benchmarking.git`

* Downloading a .ZIP archive [HERE](https://github.com/LJimenezGracia/FIXnCUT_benchmarking/archive/refs/heads/main.zip).
