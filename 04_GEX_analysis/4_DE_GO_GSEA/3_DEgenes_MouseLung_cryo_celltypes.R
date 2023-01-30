# ---
# author: "Laura Jiménez Gracia"
# date: 2021-12-04
# ---
# This R script allows to perform a Differential Expression Analysis (DEA)
# between the different groups associated to a particular condition, and
# to find transcriptional signatures associated to that condition.


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
comparison_name_list <- c("Fixed_vs_Cryo", "FixedCryo_vs_Cryo", "FixedCryo_vs_Fixed")
group_1_list <- c("Fixed", "Fixed+Cryopreserved", "Fixed+Cryopreserved")
group_2_list <- c("Cryopreserved", "Cryopreserved", "Fixed")

## Load data
seurat_obj <- readRDS(paste0(path_r_objects_in, "/FIXnCUT_MouseLung_cryopreserved_clustering_annotation_cleaned-NoNeutrophils.rds"))

for (index in 1:length(comparison_name_list)) {
  comparison_name <- comparison_name_list[[index]]
  group_1 <- group_1_list[[index]]
  group_2 <- group_2_list[[index]]
  
  # Define comparison paths
  path_object <- paste0(path_r_objects_out, "/FIXnCUT_MouseLung_cryo_DEgenes_cells_", comparison_name, ".rds")
  path_file <- paste0(path_r_tables, "/FIXnCUT_MouseLung_cryo_DEgenes_cells_", comparison_name, ".csv")
  path_file_split <- paste0(path_r_tables, "/FIXnCUT_MouseLung_cryo_DEgenes_cells_", comparison_name, "_UPDOWN.xlsx")
  
  # DE analysis considering cell-types separately
  ## Subset data
  cells_to_keep <- (seurat_obj$sample_protocol == group_1 | seurat_obj$sample_protocol == group_2)
  seurat_obj_subset <- seurat_obj[, cells_to_keep]
  
  DE_genes_bycell <- purrr::map(levels(seurat_obj_subset$celltypes), function(celltype) {
    seurat_obj_subset <- seurat_obj_subset[, seurat_obj_subset$celltypes == celltype]    
    Idents(seurat_obj_subset) <- "sample_protocol"
    
    ## Find DE genes
    DE_genes <- find_DEgenes_celltypes(
      seurat = seurat_obj_subset,
      ident_1 = group_1,
      ident_2 = group_2,
      test_de = "MAST",
      threshold_pvaladj = 0.05
    )
    
    })
  names(DE_genes_bycell) <- levels(seurat_obj_subset$celltypes)
  
  # Save results
  ## DE genes as RDS object
  saveRDS(DE_genes_bycell, file = path_object)
  
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_bycell, file = path_file)


  
  ## SPLIT up- & down-regulated
  DE_genes_split <- purrr::map(names(DE_genes_bycell), function(celltype) {
    # DEA
    DE_genes_up <- DE_genes_bycell[[celltype]] %>% 
      dplyr::filter(is_significant & avg_log2FC > 0.25) %>% 
      dplyr::select(gene, avg_log2FC, p_val_adj)
    
    DE_genes_down <- DE_genes_bycell[[celltype]] %>% 
      dplyr::filter(is_significant & avg_log2FC < -0.25)%>% 
      dplyr::select(gene, avg_log2FC, p_val_adj)
    
    DE_genes_list <- list(DE_genes_up, DE_genes_down)
    names(DE_genes_list) <- c("UP", "DOWN")
    return(DE_genes_list)
    }
    )
  names(DE_genes_split) <- names(DE_genes_bycell)
  DE_genes_split_list <- do.call(c, DE_genes_split)
  DE_genes_split_list <- DE_genes_split_list %>%
    discard(is.null)

  # Save results
  ## Export list of DE genes
  openxlsx::write.xlsx(DE_genes_split_list, file = path_file_split)
}


