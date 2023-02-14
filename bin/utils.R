# ---
# author: "Laura Jiménez Gracia"
# date: 2021-03-20
# ---
# This R script is used to define functions used in multiple scripts.

################################################################################
############################ QC Cellranger mapping #############################
################################################################################

clean_cellranger_metrics_df <- function(metrics_df) {
  colnames(metrics_df) <- str_remove(colnames(metrics_df), ":")
  colnames(metrics_df) <- str_remove(colnames(metrics_df), "\\(")
  colnames(metrics_df) <- str_remove(colnames(metrics_df), "\\)")
  colnames(metrics_df) <- str_replace_all(colnames(metrics_df), "-", "_")
  colnames(metrics_df) <- str_replace_all(colnames(metrics_df), ",", "")
  metrics_df <- as.data.frame(metrics_df)
  for (col in colnames(metrics_df)) {
    metrics_df[, col] <- str_remove(metrics_df[, col], regex("_\\(\\d+.\\d+\\%\\)"))
    if (any(str_detect(na.omit(metrics_df[, col]), "%"))) {
      metrics_df[, col] <- as.double(str_remove(metrics_df[, col], "%"))
    }
  }
  metrics_df
}


table_cellranger_metrics_gex <- function(metrics_gex_df, cellranger_version) {
  
  metrics_gex_df <- metrics_gex_df[order(metrics_gex_df$library_name), ]
  table_metrics_gex <- metrics_gex_df[, c(#"subproject_folder",
                                          "library_name",
                                          "chemistry",
                                          "lib__Number_of_reads", 
                                          "lib__Estimated_number_of_cells", 
                                          "lib__Fraction_reads_in_cells",
                                          "sample__Median_reads_per_cell",
                                          "sample__Confidently_mapped_to_exonic_regions",
                                          "sample__Median_UMI_counts_per_cell",
                                          "sample__Median_genes_per_cell")]
  
  
  table_metrics_gex %>%
    gt::gt() %>%
    fmt_percent(columns = c("lib__Fraction_reads_in_cells", "sample__Confidently_mapped_to_exonic_regions"), 
                scale_values = FALSE, decimals = 1) %>%
    fmt_number(columns = "lib__Number_of_reads", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
    tab_header(
      title = md("**GEX QC metrics**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      #subproject_folder = md("**Subproject**"),
      library_name = md("**Library**"),
      lib__Number_of_reads = md("**Number of Reads**"),
      lib__Estimated_number_of_cells = md("**Estimated Number of Recovered Cells**"),
      lib__Fraction_reads_in_cells = md("**Fraction of Reads in Cells**"),
      sample__Median_reads_per_cell = md("**Median Reads per Cell**"),
      sample__Confidently_mapped_to_exonic_regions = md("**Fraction of Reads Mapped to Exonic Reads**"),
      sample__Median_UMI_counts_per_cell = md("**Median UMI counts per Cell**"),
      sample__Median_genes_per_cell = md("**Median Genes per Cell**"),
      chemistry = md("**Chemistry**")
    )  
}


table_cellranger_metrics_vdj.t <- function(metrics_vdj.t_df, cellranger_version)  {
  table_metrics_vdj <- metrics_vdj.t_df[, c(#"subproject_folder", 
                                            "library_name",
                                            "lib__Number_of_reads", 
                                            "lib__Estimated_number_of_cells", 
                                            "lib__Fraction_reads_in_cells",
                                            "lib__Mean_reads_per_cell",
                                            "lib__Mean_used_reads_per_cell",
                                            "lib__Reads_mapped_to_any_VDJ_gene",
                                            "sample__Number_of_cells_with_productive_V_J_spanning_pair")]
  
  table_metrics_vdj %>%
    gt::gt() %>%
    fmt_percent(columns = c("lib__Fraction_reads_in_cells", "lib__Reads_mapped_to_any_VDJ_gene"), 
                scale_values = FALSE, decimals = 1) %>%
    fmt_number(columns = "lib__Number_of_reads", scale_by = 1 / 1E6, pattern = "{x}M") %>% 
    tab_header(
      title = md("**V(D)J-T QC metrics**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      #subproject_folder = md("**Subproject**"),
      library_name = md("**Library**"),
      lib__Number_of_reads = md("**Number of Reads**"),
      lib__Estimated_number_of_cells = md("**Estimated Number of Recovered Cells**"),
      lib__Fraction_reads_in_cells = md("**Fraction of Reads in Cells**"),
      lib__Mean_reads_per_cell = md("**Mean Reads per Cell**"),
      lib__Mean_used_reads_per_cell = md("**Mean (used) Reads per Cell**"),
      lib__Reads_mapped_to_any_VDJ_gene = md("**Fraction of Reads Mapped to any VDJ gene**"),
      sample__Number_of_cells_with_productive_V_J_spanning_pair = md("**Cells With Productive V-J Spanning Pair**")
    )
}


table_cellranger_annotation_vdj.t <- function(metrics_vdj.t_df, cellranger_version) {
  
  metrics_vdj.t_df %>% 
    select("subproject_folder", 
           "library_name",
           "lib__Estimated_number_of_cells",
           "sample__Number_of_cells_with_productive_V_J_spanning_pair",
           "sample__Cells_with_productive_V_J_spanning_TRA_TRB_pair", 
           "sample__Cells_with_productive_TRA_contig",
           "sample__Cells_with_productive_TRB_contig", 
           "sample__Paired_clonotype_diversity") %>% 
    gt::gt() %>%
    fmt_percent(columns = c("sample__Cells_with_productive_V_J_spanning_TRA_TRB_pair",
                            "sample__Cells_with_productive_TRA_contig", 
                            "sample__Cells_with_productive_TRB_contig"), 
                scale_values = FALSE, decimals = 1) %>%
    tab_header(
      title = md("**V(D)J-T annotation**"),
      subtitle = (cellranger_version)
    ) %>%
    cols_label(
      subproject_folder = md("**Subproject**"),
      library_name = md("**Library**"),
      lib__Estimated_number_of_cells = md("**Estimated Number of Recovered Cells**"),
      sample__Number_of_cells_with_productive_V_J_spanning_pair = md("Cells"),
      sample__Cells_with_productive_V_J_spanning_TRA_TRB_pair = md("Fraction"),
      sample__Cells_with_productive_TRA_contig = md("TRA"),
      sample__Cells_with_productive_TRB_contig = md("TRB"),
      sample__Paired_clonotype_diversity = md("**Paired Clonotype Diversity**")
    ) %>% 
    tab_spanner(
      label = md("**Productive contig**"),
      columns = vars(
        sample__Cells_with_productive_TRA_contig,
        sample__Cells_with_productive_TRB_contig)
    )   %>% 
    tab_spanner(
      label = md("**Productive V-J Spanning Pair**"),
      columns = vars(
        sample__Cells_with_productive_V_J_spanning_TRA_TRB_pair,
        sample__Number_of_cells_with_productive_V_J_spanning_pair)
    )
}



# Exploration of mapped/unmapped library' reads
plot_library_pct_reads <- function(df) {
  # Filter df
  cellranger_metrics_gex_df_sub <- df[c("library_name",
                                        "lib_unmapped", 
                                        "lib__Confidently_mapped_to_intergenic_regions",
                                        "lib__Confidently_mapped_to_intronic_regions",
                                        "lib__Confidently_mapped_to_exonic_regions")]
  # Rename df columns
  colnames(cellranger_metrics_gex_df_sub) <- c('Libraries','Unmapped','Intergenic', "Intronic", "Exonic")
  
  # Prepare df
  cellranger_metrics_gex_df_sub <- cellranger_metrics_gex_df_sub %>% 
    gather(key = "mapping_type", value = "percentage", 'Unmapped':'Exonic')
  
  # Plot
  gg_qc_map <- ggplot(cellranger_metrics_gex_df_sub, 
                      aes(x = Libraries, y = percentage, fill = mapping_type)) +
    geom_bar(position="stack", stat="identity") +
    coord_flip() + 
    theme_classic() +
    scale_fill_manual(values = c("#8DB580", "#FFEF9F", "#855A5C", "#8A8E91")) +
    labs(x = "",
         y = "Percentage of sequenced reads") +
    theme(axis.title = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          #axis.text.x = element_text(hjust = 1, angle = 90),
          strip.placement = "outside",
          strip.background = element_rect(colour = NA),
          legend.position = "top",
          legend.text = element_text(size = 12),
          legend.title=element_blank()) +
    guides(fill = guide_legend(reverse = TRUE))
  gg_qc_map
}



################################################################################
################ QC Hashtag Demultiplexing & Doublet prediction ################
################################################################################

table_genodemux_summary <- function(metrics_gex_df_hashed) {
  
  metrics_gex_df_hashed <- metrics_gex_df_hashed[, c("subproject_folder", "library_name", "demux", "demux_Freq")] %>%
    unique() %>% 
    spread(key = demux, value = demux_Freq) %>% 
    mutate(totalcells = singlet+doublet+unassigned)
  
  demux_vars <- c("subproject_folder",
                  "library_name",
                  "totalcells",
                  "singlet",
                  "doublet",
                  "unassigned")
  
  demux_global_cells <- metrics_gex_df_hashed[, demux_vars]
  
  demux_global_df <- demux_global_cells %>%
    mutate(Singlet = singlet) %>% 
    mutate(Singlet_pct = singlet/totalcells) %>%
    mutate(Multiplet = doublet) %>% 
    mutate(Multiplet_pct = doublet/totalcells) %>%    
    mutate(Unassigned = unassigned) %>% 
    mutate(Unassigned_pct = unassigned/totalcells) %>%     
    mutate(Multiplet_rate = totalcells*0.08/10000)
  demux_global_df <- demux_global_df[c("subproject_folder", "library_name", "totalcells", 
                                       "Singlet", "Singlet_pct",
                                       "Multiplet", "Multiplet_pct", "Multiplet_rate",
                                       "Unassigned", "Unassigned_pct")]
  rownames(demux_global_df) <- 1:nrow(demux_global_df)
  
  demux_global_df %>% 
    gt() %>%
    fmt_percent(
      columns = vars(Multiplet_pct, Multiplet_rate, Singlet_pct, Unassigned_pct),
      decimals = 1
    ) %>% 
    tab_header(
      title = md("**Genotype Demultiplexing QC**"),
      subtitle = md("CellSNP-lite + VIREO")
    ) %>%
    cols_label(
      subproject_folder = md("**Subproject**"),
      library_name = md("**Library**"),
      totalcells = md("**# Total cells**"),
      Multiplet = md("**# cells**"),
      Multiplet_pct = md("**Fraction**"),
      Multiplet_rate = md("*Expected*"),
      Singlet = md("**# cells**"),
      Singlet_pct = md("**Fraction**"),
      Unassigned = md("**# cells**"),
      Unassigned_pct = md("**Fraction**"),
    ) %>% 
    tab_spanner(
      label = md("**Singlets**"),
      columns = vars(
        Singlet,
        Singlet_pct)) %>% 
    tab_spanner(
      label = md("**Multiplets**"),
      columns = vars(
        Multiplet,
        Multiplet_pct,
        Multiplet_rate))  %>% 
    tab_spanner(
      label = md("**Unassigned**"),
      columns = vars(
        Unassigned,
        Unassigned_pct)
    )
}


table_demultiplexing_summary <- function(demux_global_cells) {
  demux_global_prop <- as.data.frame.matrix(prop.table(demux_global_cells, margin = 1))
  
  demux_global_df <- cbind(as.data.frame.matrix(demux_global_cells), 
                           Doublet_pct=demux_global_prop$Doublet,
                           Negative_pct=demux_global_prop$Negative,
                           Singlet_pct=demux_global_prop$Singlet)
  demux_global_df <- cbind("gem_id" = rownames(demux_global_df), demux_global_df)
  demux_global_df <- demux_global_df %>%
    group_by("gem_id") %>% 
    mutate(total_cells = Doublet+Singlet+Negative) %>% 
    mutate(doublet_rate = total_cells*0.08/10000)
  demux_global_df <- demux_global_df[c("gem_id", "total_cells", 
                                       "Singlet", "Singlet_pct",
                                       "Doublet", "Doublet_pct", "doublet_rate",
                                       "Negative", "Negative_pct")]
  rownames(demux_global_df) <- 1:nrow(demux_global_df)
  
  demux_global_df %>% 
    gt() %>%
    fmt_percent(
      columns = vars(Doublet_pct, doublet_rate, Singlet_pct, Negative_pct),
      decimals = 1
    ) %>% 
    tab_header(
      title = md("**Hashing Demultiplexing QC**"),
    ) %>%
    cols_label(
      gem_id = md("**GEM ID**"),
      total_cells = md("**# total cells**"),
      Doublet = md("**# cells**"),
      Doublet_pct = md("**Fraction**"),
      doublet_rate = md("*Expected*"),
      Singlet = md("**# cells**"),
      Singlet_pct = md("**Fraction**"),
      Negative = md("**# cells**"),
      Negative_pct = md("**Fraction**"),
    ) %>% 
    tab_spanner(
      label = md("**Singlets**"),
      columns = vars(
        Singlet,
        Singlet_pct)) %>% 
    tab_spanner(
      label = md("**Doublets**"),
      columns = vars(
        Doublet,
        Doublet_pct,
        doublet_rate))  %>% 
    tab_spanner(
      label = md("**Negative**"),
      columns = vars(
        Negative,
        Negative_pct)
    )
}


gg_nothashed_scrublet_pred <- function(seurat_metadata, gem_id, scrublet_threshold) {
  seurat_metadata %>%
    ggplot(aes(x = gem_id,
               y = scrublet_doublet_scores,
               fill = scrublet_predicted_doublet,
               colour = scrublet_predicted_doublet)) +
    geom_violin(alpha = 0.1) +
    geom_jitter(height = 0, width = 0.1) + 
    theme_bw() +
    ylim(0, 1) +    
    geom_hline(yintercept = scrublet_threshold, linetype = "dashed", colour = "red") +
    geom_text(aes(0, scrublet_threshold, label = scrublet_threshold,
                  vjust = -1, hjust = -0.4)) +
    labs(title = gem_id,
         x = "",
         y = "Doublet Score",
         fill = "Predicted Doublet",
         colour = "Predicted Doublet") +
    theme(axis.title = element_text(size = 14),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "bottom"
    )
}

gg_nothashed_umis <- function(seurat_obj) {
  ###
  # Plots the RNA UMIs distribution of singlets, doublets and negative 10X barcodes
  #
  # Parameters
  # seurat_obj
  #
  # Return
  # It returns a ggplot violin plot object with UMI distribution for HTO-global classification
  ###
  Idents(seurat_obj) <- "scrublet_predicted_doublet"
  VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.1, log = TRUE) +
    labs(x = "Predicted Doublet", y = "Number of UMIs") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12),
          plot.title = element_blank(),
          legend.position = "none"
    )
}


  
################################################################################
############################## QC Gene Expression ############################## 
################################################################################

gg_gex_horizontal_boxplot <- function(df,
                               categorical_var,
                               continuous_var,
                               fill,
                               title,
                               ylab,
                               decreasing = FALSE) {
  unordered_lev <- unique(df[[categorical_var]])
  means_cont <- purrr::map_dbl(unordered_lev, function(x) {
    mean(df[[continuous_var]][df[[categorical_var]] == x])
  })
  names(means_cont) <- unordered_lev
  ordered_lev <- names(means_cont)[order(means_cont, decreasing = decreasing)]
  df[[categorical_var]] <- factor(df[[categorical_var]], levels = ordered_lev)
  df %>%
    ggplot(aes_string(categorical_var, continuous_var, fill = fill)) +
    geom_boxplot() +
    labs(title = title, x = "", y = ylab) +
    theme_bw() +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = "bottom") +
    coord_flip()
}


table_qc_gex <- function(qc_metadata_df, subtitle) {
  qc_metadata_df %>%
    group_by(library_name) %>%
    dplyr::summarise(num_cells = n(),
                     median_counts_cell = round(median(nCount_RNA)),
                     median_genes_cell = round(median(nFeature_RNA)),
                     average_mitochondrial_fraction = mean(pct_mt) / 100
    ) %>%
    gt() %>%
    tab_header(
      title = md("**GEX QC**"),
      subtitle = subtitle
    ) %>%
    cols_label(
      library_name = md("**Library**"),
      num_cells = md("**Number of Cells**"),
      median_counts_cell = md("**Median UMI counts per Cell**"),
      median_genes_cell = md("**Median Genes per Cell**"),
      average_mitochondrial_fraction = md("**Mean Mitochondrial Fraction**")
    ) %>%
    fmt_percent(
      columns = vars(average_mitochondrial_fraction),
      decimals = 1
    )
}


################################################################################
############################# QC Doublet exclusion #############################
########################### by Ramon Massoni-Badosa ############################
################################################################################

plot_histogram_doublets <- function(df, x, x_lab, bins) {
  df %>%
    ggplot(aes_string(x)) +
    geom_histogram(bins = bins) +
    xlab(x_lab) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}

plot_density_doublets <- function(df, x, x_lab, color, color_lab) {
  df %>%
    ggplot(aes_string(x = x, color = color)) +
    geom_density() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = x_lab, color = color_lab) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11)
    )
}


plot_boxplot_doublets <- function(df, x, y, fill, y_lab) {
  df %>%
    ggplot(aes_string(x, y, fill = fill)) +
    geom_boxplot() +
    labs(x = "", y = y_lab) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 11)
    )
}


plot_scatter_doublets <- function(df, x, y, x_lab, y_lab) {
  df %>%
    ggplot(aes_string(x, y)) +
    geom_point(size = 0.1, alpha = 0.25) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = x_lab, y = y_lab) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
}


feature_plot_doublets <- function(seurat_obj, feature) {
  p <- FeaturePlot(
    seurat_obj,
    features = feature,
    cols = viridis::inferno(10),
    pt.size = 0.3
  ) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.text = element_text(size = 11)
    )
  p
}


################################################################################
################################ QC Integration ################################ 
################################################################################


gg_lisi_scores_by_confounder <- function(lisi_scores_df, confounder, title){
  lisi_scores_df %>% 
    ggplot(aes_string(x = "correction", y = confounder, fill = "correction")) +
    geom_violin() +
    stat_summary(fun = median, geom = "crossbar", size = 0.5, color= "black") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = title,
         x = "", 
         y = "LISI score") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 13, hjust = 0.5),
      strip.text = element_text(color = "black", size = 12),
      axis.title.x = element_text(color = "black", size = 13),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(size = 11)
    )    
}


################################################################################
################################# GEX analysis #################################
################################################################################

find_DEgenes <- function(seurat, ident_1, ident_2, test_de, threshold_pvaladj, latent_vars=NULL) {
  # Find DE genes
  de_genes <- Seurat::FindMarkers(
    seurat,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_de, # default wilcox
    logfc.threshold = 0, # default 0.25
    min.pct = 0.1, # default 0.1,
    latent.vars = latent_vars #default null
  )
  
  # Tag and sort significant DE genes
  de_genes <- de_genes %>% 
    rownames_to_column(var = "gene") %>%
    mutate(is_significant = ifelse(p_val_adj < threshold_pvaladj, TRUE, FALSE)) %>% 
    dplyr::arrange(desc(is_significant), desc(abs(avg_log2FC)))
  
  # Return DE genes
  de_genes
}

gg_volcano_DEgenes <- function(DEgenes_df, title, ident_1, ident_2, threshold_log2FC, threshold_adjpval) {
  # Pre-process data
  DEgenes_df <- DEgenes_df %>% 
    mutate(is_sig_direction = case_when(
      is_significant & avg_log2FC > threshold_log2FC ~ "UP",
      is_significant & avg_log2FC < -threshold_log2FC ~ "DOWN",
      TRUE ~ "FALSE"))
  
  # Subset sDE gene names
  subset_sig <- DEgenes_df %>% 
    dplyr::filter(is_significant & abs(avg_log2FC) > threshold_log2FC)
  top_up <- as.character(subset_sig$gene[subset_sig$avg_log2FC > 0][1:15])
  top_down <- subset_sig %>% 
    dplyr::arrange(avg_log2FC)
  top_down <- as.character(top_down$gene[1:15])
  subset_sig <- subset_sig %>% 
    dplyr::filter(gene %in% c(top_up, top_down))

  # Volcano plot
  volcano_plot <- DEgenes_df %>%
    ggplot(aes(avg_log2FC, -1 * log10(p_val_adj), color = is_sig_direction)) +
    geom_point(size = 3, alpha=0.7) +
    scale_color_manual(
      values = c("UP" = "red", "FALSE" = "grey", "DOWN" = "blue"),
      labels = c("UP" = "Up-regulated genes",
                 "FALSE" ="non-sDE genes",
                 "DOWN" = "Down-regulated genes"
                 )) +    
    geom_hline(yintercept = threshold_adjpval, color = "black", linetype = "dashed") +
    geom_vline(xintercept = threshold_log2FC, color = "black", linetype = "dashed") +
    geom_vline(xintercept = -threshold_log2FC, color = "black", linetype = "dashed") +
    geom_text_repel(data = subset_sig, aes(label = gene), color = "black", max.overlaps = 30, size = 6) +
    labs(title = title,
         x = "Log2(FC)",
         y = "-Log10(Adj p-value)",
         color = "") +
    theme_classic() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          legend.text = element_text(size = 17),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  volcano_plot
}


find_DEgenes_celltypes <- function(seurat, ident_1, ident_2, test_de, threshold_pvaladj, min_cells=3, latent_vars=NULL) {
  # Find DE genes
  de_genes <- Seurat::FindMarkers(
    seurat,
    ident.1 = ident_1,
    ident.2 = ident_2,
    test.use = test_de, # default wilcox
    logfc.threshold = 0, # default 0.25
    min.pct = 0.1, # default 0.1
    min.cells.group = min_cells, # default 3
    latent.vars = latent_vars #default NULL
  )
  
  # Tag and sort significant DE genes
  de_genes <- de_genes %>% 
    rownames_to_column(var = "gene") %>%
    mutate(is_significant = ifelse(p_val_adj < threshold_pvaladj, TRUE, FALSE)) %>% 
    dplyr::arrange(desc(is_significant), desc(abs(avg_log2FC)))
  
  # Return DE genes
  de_genes
}



run_enrichR_forGO <- function(database, DE_genes, adj_pval, number_enriched_genes, num_GO_size_min, num_GO_size_max, OR) {
  
  # Run enrichR
  enrichR_results <- enrichr(DE_genes, databases = database)[[1]]
  
  # Term name
  enrichR_results$term_name <- lapply(enrichR_results[["Term"]], function(x) {
    x <- str_split(x, pattern = "\\(", n = Inf)[[1]][1]
    x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
    x
  }
  )
  
  # Term GO code
  enrichR_results$code_GO <- lapply(enrichR_results[["Term"]], function(x) {
    str_extract_all(x, "(?<=\\().*(?=\\))")[[1]]}
  )  
  
  # Number of DE genes enriched for a specific GO
  enrichR_results$number_genes <- lapply(enrichR_results[["Overlap"]], function(x) {
    str_split(x, pattern = "/", n = Inf)[[1]][1]})
  
  # Number of genes within a specific GO
  enrichR_results$GO_size <- lapply(enrichR_results[["Overlap"]], function(x) {
    str_split(x, pattern = "/", n = Inf)[[1]][2]})
  
  # Percentage overlap
  enrichR_results$pct_overlap <- sapply(enrichR_results$Overlap, function(x) 
    eval(parse(text=x)))
  
  # Convert dataframe results
  enrichR_results <- as.tibble(enrichR_results) %>% 
    transform(term_name = as.character(term_name),
              code_GO = as.character(code_GO),
              number_genes = as.numeric(number_genes),
              GO_size = as.numeric(GO_size))
  
  # Tag GO as reliable TRUE/FALSE
  enrichR_results <- enrichR_results %>%
    mutate(is_reliable = ifelse(
      (Adjusted.P.value < adj_pval & 
         number_genes >= number_enriched_genes & 
         GO_size >= num_GO_size_min &
         GO_size <= num_GO_size_max & 
         Odds.Ratio  >= OR), TRUE, FALSE)) %>%
    dplyr::arrange(desc(is_reliable), Adjusted.P.value)
  
  enrichR_results
}

gg_enrichR_GO <- function(comparison, database, enrichR_results_filtered, n_top) {
  
  enrichR_results_filtered <- enrichR_results_filtered[1:n_top,][, c("Odds.Ratio", "term_name", "Adjusted.P.value", "number_genes")]
  
  gg_plot <- enrichR_results_filtered %>%
    ggplot(., aes(
      x = Odds.Ratio,
      y = fct_reorder(term_name, -Adjusted.P.value),
      size = number_genes,
      color = Adjusted.P.value)) +
    geom_point() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "red") +
    scale_color_gradient2(
      name = "Adjusted p-value",
      low = "#bf212f",
      mid = "#006f3c",
      high = "lightgrey") +
    theme_minimal() +
    labs(title = comparison,
         x = "Odds Ratio",
         y = paste0(database, "Term"),
         size = "# Genes") +
    theme(
      axis.text.y = element_text(size=12)
    )
  
  gg_plot
}



# This function converts a vector of gene identifiers to other identifiers 
convert_geneid <- function(gene_vec,
                           gene_from  = "SYMBOL", # "SYMBOL", "ENSEMBL", "ENTREZID", "UNIGENE"
                           gene_to  = "ENTREZID", # "SYMBOL", "ENSEMBL", "ENTREZID", "UNIGENE"
                           annotation = "org.Hs.eg.db" # org.Hs.eg.db or org.Mm.eg.db
                           ) {
  
  # Run separate function for human or mouse
  if (annotation == "org.Hs.eg.db") {
    # Convert symbols to ENTREZID
    de_entrezid <- AnnotationDbi::mapIds(
      x = org.Hs.eg.db::org.Hs.eg.db,
      keys = unlist(gene_vec),
      #   multiVals = "first",
      column = gene_to,
      keytype = gene_from
      )
    
  } else if (annotation == "org.Mm.eg.db") {
    # Convert symbols to ENTREZID
    de_entrezid <- AnnotationDbi::mapIds(
      x = org.Mm.eg.db::org.Mm.eg.db,
      keys = unlist(gene_vec),
      #   multiVals = "first",
      column = gene_to,
      keytype = gene_from
      )
    
  }
  
  de_entrezid <- de_entrezid[!is.na(de_entrezid)]
  return(de_entrezid)
  
}

#' This function carries out a GO enrichment analysis between a set of
#' differentially expressed genes within the specified gene universe.

GO_enrichment <- function(DE_genes,
                               gene_universe,
                               gene_from  = "SYMBOL", # SYMBOL, ENSEMBL, ENTREZID, or UNIGENE
                               gene_to = "ENTREZID", # SYMBOL, ENSEMBL, ENTREZID, or UNIGENE
                               pvalue_cutoff = 0.05,
                               test_direction = "over",
                               annotation = "org.Hs.eg.db", # org.Hs.eg.db or org.Mm.eg.db
                               ontology = "BP" # BP, CC or MF
                          ) {

  # Convert Target genes from symbol/ensembl to entrezID
  target_enrich <- convert_geneid(gene_vec = unique(DE_genes),
                                  gene_from = gene_from,
                                  gene_to = gene_to,
                                  annotation = annotation)
  
  # Convert Universe genes from symbol/ensembl to entrezID
  univers_g <- convert_geneid(gene_vec = unique(gene_universe),
                              gene_from = gene_from,
                              gene_to = gene_to,
                              annotation = annotation)
  
  # Create a GOHyperGParams instance
  params <- methods::new(
    "GOHyperGParams",
    geneIds = target_enrich,
    universeGeneIds = univers_g,
    annotation = annotation,
    ontology = ontology,
    pvalueCutoff = pvalue_cutoff,
    testDirection = test_direction)
  
  # Carry out hyper geometric test
  hgOver <- GOstats::hyperGTest(params)
  hgOver_df <- summary(hgOver)
  go_results <- hgOver_df %>% 
    arrange(desc(OddsRatio))
  return(go_results)
  
}

# 
# # Computing parameters
# genes <- rownames(seurat[["RNA"]]@data)
# avg_expr <- rowMeans(as.matrix(seurat[["RNA"]]@data))
# log2_fc <- apply(as.matrix(seurat[["RNA"]]@data), 1, function(x) {
#   mean_KO <- mean(x[seurat$hash.ID == ident_1]) + 1
#   mean_wt <- mean(x[seurat$hash.ID == ident_2]) + 1
#   log2(mean_KO / mean_wt)
# })
# 
# # Set up dataframe
# meanexp_df <- data.frame(
#   gene = genes, 
#   average_expression = avg_expr, 
#   log2_fc = log2_fc)
# deg_general <- merge(
#   deg_general,
#   meanexp_df, 
#   by.x = "gene", by.y = "gene") %>%
#   
#   
#   
#   # Functions
#   

# 
# 
# 
# DEanalysis_MAplot_general <- function(DEgenes_df, ident, ident_1, ident_2, threshold_log2FC, threshold_avexpr) {
#   # Subset sDE gene names
#   deg_subset_sig <- dplyr::filter(DEgenes_df, is_significant & 
#                                     log(average_expression + 1) > threshold_avexpr & 
#                                     abs(log2_fc) > threshold_log2FC)
#   
#   deg_subset_sig_up <- deg_subset_sig %>% 
#     dplyr::arrange(desc(log2_fc)) %>% 
#     filter(log2_fc > 0, preserve = TRUE)
#   deg_subset_sig_up_names <- as.character(
#     slice_max(deg_subset_sig_up, order_by = log2_fc, n = 10)$gene)
#   
#   deg_subset_sig_down <- deg_subset_sig %>% 
#     dplyr::arrange(log2_fc) %>% 
#     filter(log2_fc < 0, preserve = TRUE)
#   deg_subset_sig_down_names <- as.character(
#     slice_min(deg_subset_sig_down, order_by = log2_fc, n = 10)$gene)
#   
#   deg_subset_sig <- deg_subset_sig %>% 
#     dplyr::filter(gene %in% c(deg_subset_sig_up_names, deg_subset_sig_down_names))
#   
#   # MA plot
#   maplot_general <- DEgenes_df %>% 
#     ggplot(aes(log2(average_expression + 1), log2_fc, color = is_significant)) +
#     geom_point(size = 1) +
#     geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
#     geom_smooth(method = "loess", color = "darkblue") +
#     labs(title = paste0(ident_1, " vs ", ident_2),
#          x = "Log2(Average expression)",
#          y = "Log2FC",
#          color = "") +
#     theme_classic() +
#     scale_color_manual(
#       values = c("darkgrey", "lightgreen"), 
#       labels = c("non-sDE genes",
#                  "sDE genes")) +
#     theme(axis.title = element_text(size = 11),
#           legend.text = element_text(size = 11),
#           plot.title = element_text(hjust = 0.5, face = "bold")) +
#     geom_text_repel(data = deg_subset_sig, aes(label = gene), color = "black", max.overlaps = 20)  
#   
#   # Save ggplot
#   ggsave(filename = paste0("results/DEanalysis/cd45pos_maplot_DEAgeneral_", ident, ".png"),
#          plot = maplot_general)
#   
#   print(maplot_general)
# }


################################################################################
############################### General Functions ############################## 
################################################################################


cpdb_convert_seurat_to_h5ad <- function(seurat_obj, intermediate_h5_file_path){
  # Remove unwanted data assays and slots
  seurat_obj <- DietSeurat(
    object = seurat_obj,
    counts = FALSE,
    data = TRUE,
    scale.data = FALSE,
    assays = "RNA"
  )
  
  # Convert-seurat-to-anndata-object
  SeuratDisk::SaveH5Seurat(seurat_obj, 
                           filename = intermediate_h5_file_path)
  
  # Then, convert the h5 file into an AnnData object
  SeuratDisk::Convert(intermediate_h5_file_path, 
                      dest = "h5ad", 
                      overwrite = TRUE)
  
  # Remove Seurat h5 intermediate object
  base::system2('rm', args = intermediate_h5_file_path)
}



convert_genes_human_to_mouse <- function(hsapiens_gene_list){
  ###
  # Converts human genes to mouse genes
  #
  # Parameters
  # hsapiens_gene_list: a list of human genes to be converted
  #
  # Return
  # It returns a list of mouse genes
  ###
  require(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol",
                   values = hsapiens_gene_list,
                   mart = human,
                   attributesL = c("mgi_symbol"),
                   martL = mouse,
                   uniqueRows = T)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  row.names(genesV2) <- genesV2$HGNC.symbol
  mmusculus_gene_list <- genesV2[x, 2 ]
  
  return(mmusculus_gene_list)
}


convert_genes_mouse_to_human <- function(mmusculus_gene_list){
  ###
  # Converts mouse genes to human genes
  #
  # Parameters
  # mmusculus_gene_list: a list of mouse genes to be converted
  #
  # Return
  # It returns a list of human genes
  ###
  require(biomaRt)
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol",
                   values = mmusculus_gene_list,
                   mart = mouse,
                   attributesL = c("hgnc_symbol"),
                   martL = human,
                   uniqueRows = T)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  row.names(genesV2) <- genesV2$HGNC.symbol
  hsapiens_gene_list <- genesV2[x, 2 ]
  
  return(hsapiens_gene_list)
}



# 
#' This function takes in a Seurat 3.0 object and returns a named list with 2
#' objects formated to be loaded in the ShinyApp:
#' https://singlecellgenomics-cnag-crg.shinyapps.io/Annotation/
#' 
#' 1st the metadata + the coordinates of the 2D embeding, 
#' the later with names coord_x and coord_y, and second the expression matrix selected.
#'
#' @param se_obj Object of class Seurat from which we want to extract the information.
#' @param assay Object of class Character indicating from which assay to extract expression data.
#' @param slot Object of class Character indicating from which slot to extract expression data, by default data.
#' @param reduction Object of class Character indicating from which dimensionality reduction we want to extract the coordinates.
#' @return This function returns a named list, the first position contains the joint metadata + 2D embeding and the second contains the expression data matrix. 2D dimensions are labelles as coord_x and coord_y.
#' @export
#' @examples
#'
#'
prepare_se_obj <- function(se_obj,
                           assay,
                           slot = "data",
                           reduction = "umap") {
  suppressMessages(library(dplyr))
  suppressMessages(library(Seurat))
  suppressMessages(library(tibble))
  suppressMessages(library(SummarizedExperiment))
  
  # Input check
  if (! is(se_obj, "Seurat")) stop("Object passed is not a Seurat object;")
  if (! assay %in% Seurat::Assays(se_obj)) stop("assay not in the Seurat object's available assays.")
  if (! slot %in% c("counts", "data", "scale.data")) warning("slot not in the Seurat object's assays.")
  if (! reduction %in% names(se_obj@reductions)) stop("reduction not in the Seurat object's available reductions.")
  
  # Extract 2D cordinates
  embed_df <- se_obj@reductions[[reduction]]@cell.embeddings %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode")
  
  # Change 2D coordinate names to coord_x and coord_y
  colnames(embed_df) <- c("barcode", "coord_x", "coord_y")
  
  # Join metadata with coordinates
  metadata_df <- se_obj@meta.data %>%
    data.frame() %>%
    tibble::rownames_to_column("barcode") %>%
    dplyr::left_join(embed_df, by = "barcode")
  
  # Extract expression data
  assay_data <- Seurat::GetAssayData(se_obj,
                                     slot = slot,
                                     assay = assay)
  
  return(list(metadata = metadata_df,
              expr_data = assay_data))
}



################################################################################
############################### General Variables ############################## 
################################################################################

