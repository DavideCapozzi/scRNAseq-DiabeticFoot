#!/bin/bash

# --- Opzioni SLURM ---
#SBATCH --job-name=violin_plots          # Nome del job
#SBATCH --output=logs/violin_%j.log      # File di log (l'output del terminale finirà qui)    # File di errori
#SBATCH --nodes=1                        # Numero di nodi
#SBATCH --error=logs/violin_%j.err
#SBATCH --ntasks=1                       # Numero di task
#SBATCH --cpus-per-task=8                # Numero di CPU (Seurat beneficia di un po' di multi-threading)
#SBATCH --mem=32G                        # Memoria RAM (gli oggetti Seurat sono pesanti, 32GB dovrebbero bastare)
#SBATCH --time=02:00:00                  # Tempo massimo (2 ore sono più che sufficienti per dei plot)

# Oppure se il cluster usa i moduli:
module load R/4.2.1 

# --- Esecuzione ---
echo "Inizio Job alle: $(date)"

# Vai nella directory del progetto
cd /lustre/home/fdannunzio/diabetic_foot

# Crea una cartella per i log se non esiste
mkdir -p logs

# Lancia lo script R (assumendo che tu l'abbia salvato come script_violin.R)
Rscript 5_make_violin_plot_selected_proteins.R

echo "Fine Job alle: $(date)"