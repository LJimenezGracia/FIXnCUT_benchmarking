#!/usr/bin/env bash

#SBATCH --job-name="SCGTEST33_LUNG"
#SBATCH --workdir=/scratch/devel/ljimenez/projects/FIXnCUT/02_QC/1_downsampling_reads

#SBATCH --error=logs/slurm_%x_%J.err
#SBATCH --output=logs/slurm_%x_%J.out

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00

#SBATCH --mail-type=END       
#SBATCH --mail-user=laura.jimenez@cnag.crg.eu

echo [`date "+%Y-%m-%d %T"`] started job on $HOSTNAME

ulimit -n 16000
export HDF5_USE_FILE_LOCKING="FALSE"
export TENX_IGNORE_DEPRECATED_OS=1 

/scratch/groups/singlecell/software/cellranger/6.1.1/cellranger aggr --id=SCGTEST33_LUNG --csv=/scratch/devel/ljimenez/projects/FIXnCUT/02_QC/1_downsampling_reads/data/SCGTEST33_LUNG_aggr.csv --normalize=mapped --nosecondary

echo [`date "+%Y-%m-%d %T"`] job finished
