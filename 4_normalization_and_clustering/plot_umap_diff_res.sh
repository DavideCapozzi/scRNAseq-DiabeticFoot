#!/bin/bash
#SBATCH -A sens2025518
#SBATCH --job-name=plot_umap_diff_res
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=10-00:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=plot_umap_diff_res_%j.log

module load bioinfo-tools
module load Seurat/5.2.1  

Rscript /home/fdann/Desktop/proj/sc-res/2_new_analysis/4_normalization_and_clustering/3_plot_umap_cluster_different_res.R
