# ---
# author: "Laura Jiménez Gracia"
# date: 2021-04-09
# ---
# This R script allows to find specific biomarkers for each of the defined clusters.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(magrittr)

## Paths
path_r_objects <- here::here("03_clustering_annotation/results/R_objects")
path_r_tables <- here::here("03_clustering_annotation/results/tables")

## Load Seurat object
seurat_obj <- readRDS(paste0(path_r_objects, "/omniscope39_m_LUNG_healthy_stress_clustering_resolutions.rds"))
seurat_obj


# Define Cluster biomarkers
## Ensure we are using RNA assay
DefaultAssay(seurat_obj) <- "RNA"

## Select cluster resolution
Idents(seurat_obj) <- "RNA_snn_res.0.3"

## Find diferentially expressed features
seurat_obj_markers <- seurat_obj %>% 
  FindAllMarkers(
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox")

## Filtering out markers with adjusted p-value > 0.05
seurat_obj_markers %<>%
  dplyr::filter(p_val_adj < 0.05)


# Save cell-type biomarkers
## Export in .rds object
saveRDS(seurat_obj_markers, 
        file = paste0(path_r_objects, "/biomarkers_omniscope39_m_LUNG_healthy_stress_resolution0.3.rds"))

## Export in .xlsx format
### Prepare data: sort and filter
biomarkers_df <- seurat_obj_markers %>%
  dplyr::arrange(cluster, desc(abs(avg_log2FC))) %>%
  dplyr::filter(avg_log2FC > 0.5)

### Convert biomarkers dataframe into a biomarker list by cluster
biomarkers_list <- purrr::map(levels(biomarkers_df$cluster),
                              ~ biomarkers_df[biomarkers_df$cluster == .x, ])
names(biomarkers_list) <- levels(biomarkers_df$cluster)

openxlsx::write.xlsx(biomarkers_list,
                     file = paste0(path_r_tables, "/biomarkers_omniscope39_m_LUNG_healthy_stress_resolution0.3.xlsx"))
