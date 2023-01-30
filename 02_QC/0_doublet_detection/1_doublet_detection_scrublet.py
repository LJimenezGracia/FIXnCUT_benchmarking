#! /usr/bin/env python

"""
Author: "Laura Jiménez Gracia"
Date: 2021-01-19

Doublet Detection with Scrublet

Description:
This script reads the filtered_feature_bc_matrix from a specific GEM channel (library) 
and runs scrublet. Moreover, it saves a .csv file with the doublet score and doublet 
prediction for each cell, the histograms of their distribution and the UMAPs.

Paramater: 
    Expected doublet rate -- From the 10X Genomics protocols, it is known that 
        there is a linear relationship between the % of multiplets and 
        the # of recovered cells per 10x Chromium Chip Channel (GEM well),
        with a slope of 8% multiplets for 10,000 cells.
"""

# Import packages
import os
import pandas as pd
import numpy as np
import scipy.io
import scrublet as scr
import matplotlib.pyplot as plt

# Initialize variables
subproject = "THREETR_27"
subproject_path = "/mnt/cnagclus/projects/FIXnCUT/01_cellranger_mapping/subprojects/{}".format(subproject)
metadata_path = "/mnt/cnagclus/projects/FIXnCUT/01_cellranger_mapping/data/FIXnCUT_metadata.csv"
scrublet_results_path = "/mnt/cnagclus/projects/FIXnCUT/02_QC/0_doublet_detection/results"


def main():
    """
    """
    
    # Define list of libraries to merge 'metrics'
    metadata_df = pd.read_csv(metadata_path, sep=",", header=0)
    mask = (metadata_df["subproject"] == subproject)
    libraries = metadata_df.loc[mask, "gem_id"]
    libraries_list = list(set(libraries))

    # Create results folder
    if not os.path.isdir(scrublet_results_path):
        os.mkdir(scrublet_results_path)

    for gem_id in libraries_list:
        print("\n################## ", gem_id, " ##################")

        # Load data for GEX + TCR
        counts_path = "{}/jobs/{}/{}/outs/per_sample_outs/{}/count/sample_feature_bc_matrix/matrix.mtx.gz".format(
            subproject_path, gem_id, gem_id, gem_id)
        barcodes_path = "{}/jobs/{}/{}/outs/per_sample_outs/{}/count/sample_feature_bc_matrix/barcodes.tsv.gz".format(
            subproject_path, gem_id, gem_id, gem_id)               
        counts_matrix = scipy.io.mmread(counts_path).T.tocsc()
        barcodes_df = pd.read_csv(barcodes_path, header=None)
        n_cells = counts_matrix.shape[0]
	    #n_genes = counts_matrix.shape[1]
        print("Number of cells: {}".format(n_cells))

        # Parameters
        expected_doublet_rate = round((n_cells * 0.08 / 10000), 3)
        print("Expected multiplet rate: {}%\n".format(round(expected_doublet_rate*100, 1)))

        ## Initialize scrublet object
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)

        # Run doublet simulation (with default values)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=3,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30)
        doublet_scores = np.round(doublet_scores, decimals = 3)

        ## Patch to avoid failing scrublet whem no automatic threshold is found
        ## setting threshold to 1, and no doublets will be detected
        if scrub.predicted_doublets_ is None:
            print("No doublets predicted using default method.")
            scrub.call_doublets(threshold=1)

        # Create data frame
        scrublet_doubl_dict = {"barcodes": barcodes_df[0].values,
                            "scrublet_doublet_scores": doublet_scores,
                            "scrublet_predicted_doublet": predicted_doublets}
        scrublet_doubl_df = pd.DataFrame(scrublet_doubl_dict)

        # Get 2D embedding to visualize the results
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

        # Plots
        hist = scrub.plot_histogram()

        umap = scrub.plot_embedding('UMAP', order_points=True)

        # Save all data
        if not os.path.exists(scrublet_results_path):
            os.mkdir(scrublet_results_path)

        scrublet_doubl_df.to_csv("{}/{}__{}_doublet_prediction_scrublet.csv".format(
            scrublet_results_path, subproject, gem_id), index=False)

        hist[0].savefig("{}/{}__{}_doublet_prediction_scrublet_hist.png".format(
            scrublet_results_path, subproject, gem_id), dpi = 100)       
        
        #umap[0].savefig("{}/{}__{}_doublet_prediction_scrublet_umap.png".format(
        #    scrublet_results_path, subproject, gem_id), dpi = 100)
            
        plt.close('all')

if __name__ == "__main__":
    main()
