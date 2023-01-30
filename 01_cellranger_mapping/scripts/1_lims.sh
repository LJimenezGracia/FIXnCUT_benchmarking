#!/bin/bash

# Load required modules
module load PYTHON/2.7.5
module load lims/1.2

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project

subprojects="SCGTEST_32,SCGTEST_33,SCGTEST_34,OMNISCOPE_39,THREETR_25,OMNISCOPE_56,OMNISCOPE_58,SALAS_47,THREETR_26,THREETR_27" # separate by comma
/scratch/project/production/DAT/apps/LIMSQ/limsq -sp $subprojects | sed 's/;/\t/g' > ../data/project_FIXnCUT_info.tsv
