#! /usr/bin/env python3

"""
Author: "Laura Jiménez Gracia"
Date: 2021-11-05

This script merges performance metrics of cellranger multi [>= version 6] for all libraries (GEM ids) 
in this subproject into a single file cellranger_mapping_metrics.csv located in /results folder.
"""

# Load packages
import os
import pandas as pd
import numpy as np
import config as cfg


def main():

    # Get paths
    project_path = cfg.project_path
    subprojects_path = cfg.subprojects_path
    metadata_path = cfg.metadata_path

    # Define list of libraries to merge 'metrics'
    metadata_df = pd.read_csv(metadata_path, sep=",", header=0)

    # Define dataframes
    summary_gex_lib_dfs = pd.DataFrame()
    summary_gex_sample_dfs = pd.DataFrame()
    summary_gex = pd.DataFrame()
    summary_vdjt_lib_dfs = pd.DataFrame()
    summary_vdjt_sample_dfs = pd.DataFrame()
    summary_vdjt = pd.DataFrame()


    subproject_list = np.unique(metadata_df["subproject"])
    for subproject in subproject_list:
        gem_id_list = list(metadata_df.loc[(metadata_df["subproject"] == subproject), "gem_id"].unique())
        for gem_id in gem_id_list:
            ## FOR GEX + VDJ libraries
            metrics_path=f"{subprojects_path}/{subproject}/jobs/{gem_id}/{gem_id}/outs/per_sample_outs/{gem_id}/metrics_summary.csv"
            if os.path.exists(metrics_path):
                metrics_df = pd.read_csv(metrics_path)
                metrics_df.columns = metrics_df.columns.str.replace(' ','_')
                for i in metrics_df.columns:
                    metrics_df[i] = metrics_df[i].str.replace(' ','_')
            
                # Parsing GEX information
                # Subsetting GEX LIBRARY
                gex_lib_df = metrics_df[
                    (metrics_df['Library_or_Sample'] == 'Library') &
                    (metrics_df['Library_Type'] == 'Gene_Expression')][['Metric_Name','Metric_Value']]
                gex_lib_df = gex_lib_df.set_index('Metric_Name')
                gex_lib_df.index.names = [None]
                gex_lib_df.columns = [gem_id]
                gex_lib_df_t = gex_lib_df.transpose()
                gex_lib_df_t.reset_index(inplace=True)
                gex_lib_df_t = gex_lib_df_t.rename(columns = {'index':'gem_id'})
                summary_gex_lib_dfs = summary_gex_lib_dfs.append(gex_lib_df_t)

                # Subsetting GEX SAMPLE
                gex_sample_df = metrics_df[
                    (metrics_df['Library_or_Sample'] == 'Sample') &
                    (metrics_df['Library_Type'] == 'Gene_Expression')][['Metric_Name','Metric_Value']]
                gex_sample_df = gex_sample_df.set_index('Metric_Name')
                gex_sample_df.index.names = [None]
                gex_sample_df.columns = [gem_id]
                gex_sample_df_t = gex_sample_df.transpose()
                gex_sample_df_t.reset_index(inplace=True)
                gex_sample_df_t = gex_sample_df_t.rename(columns = {'index':'gem_id'})
                summary_gex_sample_dfs = summary_gex_sample_dfs.append(gex_sample_df_t)


                # Parsing VDJ-T information
                chemistry = metadata_df.loc[(metadata_df["gem_id"] == gem_id), "chemistry"].unique()
                if (chemistry == "SC5P-R2"):
                    # Subsetting VDJ-T LIBRARY
                    vdjt_lib_df = metrics_df[
                        (metrics_df['Library_or_Sample'] == 'Library') &
                        (metrics_df['Library_Type'] == 'VDJ_T')][['Metric_Name','Metric_Value']]
                    vdjt_lib_df = vdjt_lib_df.set_index('Metric_Name')
                    vdjt_lib_df.index.names = [None]
                    vdjt_lib_df.columns = [gem_id]
                    vdjt_lib_df_t = vdjt_lib_df.transpose()
                    vdjt_lib_df_t.reset_index(inplace=True)
                    vdjt_lib_df_t = vdjt_lib_df_t.rename(columns = {'index':'gem_id'})
                    summary_vdjt_lib_dfs = summary_vdjt_lib_dfs.append(vdjt_lib_df_t)

                    # Subsetting VDJ-T SAMPLE
                    vdjt_sample_df = metrics_df[
                        (metrics_df['Library_or_Sample'] == 'Sample') &
                        (metrics_df['Library_Type'] == 'VDJ_T')][['Metric_Name','Metric_Value']]
                    vdjt_sample_df = vdjt_sample_df.set_index('Metric_Name')
                    vdjt_sample_df.index.names = [None]
                    vdjt_sample_df.columns = [gem_id]
                    vdjt_sample_df_t = vdjt_sample_df.transpose()
                    vdjt_sample_df_t.reset_index(inplace=True)
                    vdjt_sample_df_t = vdjt_sample_df_t.rename(columns = {'index':'gem_id'})
                    summary_vdjt_sample_dfs = summary_vdjt_sample_dfs.append(vdjt_sample_df_t)
      

    # Formatting dataframes
    if not summary_gex_lib_dfs.empty:
        summary_gex_lib_dfs.drop_duplicates(subset = "gem_id", inplace=True)
        summary_gex_lib_dfs = summary_gex_lib_dfs.T.drop_duplicates().T
        summary_gex_lib_dfs = summary_gex_lib_dfs.add_prefix('lib__')
        summary_gex_lib_dfs = summary_gex_lib_dfs.rename(columns={"lib__gem_id": "gem_id"})
        summary_gex_sample_dfs.drop_duplicates(subset = "gem_id", inplace=True)
        summary_gex_sample_dfs = summary_gex_sample_dfs.T.drop_duplicates().T
        summary_gex_sample_dfs = summary_gex_sample_dfs.add_prefix('sample__')
        summary_gex_sample_dfs = summary_gex_sample_dfs.rename(columns={"sample__gem_id": "gem_id"})

        summary_gex = summary_gex_lib_dfs.merge(summary_gex_sample_dfs, how = "left", on="gem_id")

    if not summary_vdjt_lib_dfs.empty:
        summary_vdjt_lib_dfs.drop_duplicates(subset = "gem_id", inplace=True)
        summary_vdjt_lib_dfs = summary_vdjt_lib_dfs.T.drop_duplicates().T
        summary_vdjt_lib_dfs = summary_vdjt_lib_dfs.add_prefix('lib__')
        summary_vdjt_lib_dfs = summary_vdjt_lib_dfs.rename(columns={"lib__gem_id": "gem_id"})
        summary_vdjt_sample_dfs.drop_duplicates(subset = "gem_id", inplace=True)
        summary_vdjt_sample_dfs = summary_vdjt_sample_dfs.T.drop_duplicates().T
        summary_vdjt_sample_dfs = summary_vdjt_sample_dfs.add_prefix('sample__')
        summary_vdjt_sample_dfs = summary_vdjt_sample_dfs.rename(columns={"sample__gem_id": "gem_id"})

        summary_vdjt = summary_vdjt_lib_dfs.merge(summary_vdjt_sample_dfs, how = "left", on="gem_id")

    # Save combined metrics
    # Create results directory
    results_directory = f"{project_path}/results"
    if not os.path.exists(results_directory):
        os.mkdir(results_directory)

    if not summary_gex.empty:
        summary_gex.to_csv(f"{results_directory}/cellranger_mapping_metrics_GEX.csv", header = True, index = None)
    
    if not summary_vdjt.empty:
        summary_vdjt.to_csv(f"{results_directory}/cellranger_mapping_metrics_VDJT.csv", header = True, index = None)



if __name__ == '__main__':
    exit(main())
