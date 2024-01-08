# ---
# author: "Laura Jiménez Gracia"
# date: 2021-12-04
# ---
# This R script allows to perform a Differential Expression Analysis (DEA)
# between the different groups associated to a particular sample_protocol, and
# to find transcriptional signatures associated to that sample_protocol.


# Pre-processing
## Load packages
library(tidyverse)
library(Seurat)
library(ggrepel)

## Parameters
# Paths
path_r_objects_in <- here::here("03_clustering_annotation/results/R_objects")
path_r_tables <- here::here("04_GEX_analysis/4_DE_GO_GSEA/results/tables")
path_r_objects_out <- here::here("04_GEX_analysis/4_DE_GO_GSEA/results/R_objects")

# Functions
source(here::here("bin/utils.R"))

# Define comparison parameters
comparison_name_list <- c("Fixed_vs_Fresh")
group_1_list <- c("Fixed")
group_2_list <- c("Fresh")

## Load data
seurat_obj <- readRDS(paste0(path_r_objects_in, "/FIXnCUT_MouseLung_fixed_perfused_clustering_annotation_cleaned.rds"))


for (index in 1:length(comparison_name_list)) {
  comparison_name <- comparison_name_list[[index]]
  group_1 <- group_1_list[[index]]
  group_2 <- group_2_list[[index]]
  
  # Define comparison paths
  path_object <- paste0(path_r_objects_out, "/FIXnCUT_MouseLung_fixed_perfused_DEgenes_", comparison_name, ".rds")
  path_file <- paste0(path_r_tables, "/FIXnCUT_MouseLung_fixed_perfused_DEgenes_", comparison_name, ".csv")
  path_file_split <- paste0(path_r_tables, "/FIXnCUT_MouseLung_fixed_perfused_DEgenes_", comparison_name, "_UPDOWN.xlsx")
  
  # DE analysis considering all cell-types together
  ## Subset data
  cells_to_keep <- (seurat_obj$sample_protocol == group_1 | seurat_obj$sample_protocol == group_2)
  seurat_obj_sub <- seurat_obj[, cells_to_keep]
  Idents(seurat_obj_sub) <- "sample_protocol"
  
  ## Find DE genes
  DE_genes <- find_DEgenes(
    seurat = seurat_obj_sub,
    ident_1 = group_1,
    ident_2 = group_2,
    test_de = "MAST",
    #latent_vars = "replicate",
    threshold_pvaladj = 0.05
  )
  
  # Save results
  ## DE genes as RDS object
  saveRDS(DE_genes, file = path_object)
  
  ## List of DE genes as csv file
  write.csv(DE_genes, file = path_file, row.names=FALSE)
  
  
  ## SPLIT up- & down-regulated
  DE_genes_up <- DE_genes %>% 
    dplyr::filter(is_significant & avg_log2FC > 0) %>% 
    dplyr::select(gene, avg_log2FC, p_val_adj)
  
  DE_genes_down <- DE_genes %>% 
    dplyr::filter(is_significant & avg_log2FC < 0)%>% 
    dplyr::select(gene, avg_log2FC, p_val_adj)
  
  name_group_1 <- paste0(group_1, " UP")
  name_group_2 <- paste0(group_2, " DOWN")
  DE_genes_split <- list(DE_genes_up, DE_genes_down)
  names(DE_genes_split) <- c(name_group_1, name_group_2)
  
  # Save results
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_split, file = path_file_split)
  
}
