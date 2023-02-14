# FIXnCUT_benchmarking

This is the GitHub repository for the **FixNCut: Digesting fixed tissues for single-cell genomics** manuscript.

*This manuscript is current under review process.*


## Abstract

The expansion of single-cell applications in the era of personalized medicine and the increasing complexity and decentralization of current studies require disconnecting the time and site of sampling from downstream processing steps. In addition, cellular transcriptomic profiles are dynamic and can change under external stressors, altering the natural state of the cells. In light of this, early sample preservation improves robustness and reproducibility by avoiding artifacts introduced during sample handling. Here, we present FixNCut, an approach for reversibly fixing tissue followed by dissociation to overcome current limitations in the generation of single-cell data. We have applied FixNCut prior to single-cell RNA sequencing for mouse lung and colon tissue, sample types with inherent sampling challenges. Further, we processed human colon biopsies from pathological context to provide evidence of the clinical utility of this methodology. We showed that reversible fixation followed by dissociation enables the removal of time and location constraints, while preserving the RNA integrity, library complexity and cellular composition as well as diminishing stress-related artifacts associated with sample processing. Moreover, fixed cells can be stained with antibodies and lipid modified oligos (LMOs), as proved by cytofluorimetric analysis, thus allowing fluorescence-activated cell sorting (FACS) or hashing prior to single-cell analysis. In principle, FixNCut is compatible with most standard single-cell technologies, making it a versatile protocol to enable robust and flexible study designs.


## Code implementation

The repository is organized into the following folder tree, which contain the necessary files to perform all detailed tasks and obtained all our results.

* **01_cellranger_mapping** -->
* **02_QC** -->
* **03_clustering_annotation** -->
* **04_GEX_analysis** -->
* **bin** -->


## Data accessibility

* The complete raw data (FASTQ files) generated in this study have been submitted to the NCBI Gene Expression Omnibus (GEO) under accession number XXXXX.
* The count matrices and metadata are deposited at [Zenodo](). 


## Code accessibility

You can easily download a copy of all the files contained in this repository by:

* Cloning the git repository using the following command in the terminal:

`git clone https://github.com/LJimenezGracia/FIXnCUT_benchmarking.git`

* Downloading a .ZIP archive [HERE](https://github.com/LJimenezGracia/FIXnCUT_benchmarking/archive/refs/heads/main.zip).
