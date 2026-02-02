#!/bin/bash
#SBATCH -A sens2025518
#SBATCH --job-name=cellranger
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=336:00:00
#SBATCH --mail-user=federica.dannunzio@uniroma1.it
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=cellranger_%j.log

# ============================================
CELLRANGER_DIR="/home/fdann/Desktop/wharf/fdann/fdann-sens2025518/pre-processing/cellranger-9.0.1/bin"
export PATH="$CELLRANGER_DIR:$PATH"

# ==============================================
# 1. REFERENCE GENOME SETUP
# ==============================================
REF_DIR="/home/fdann/Desktop/wharf/fdann/fdann-sens2025518/pre-processing"  # Reference will be stored here
REF_PATH="$REF_DIR/refdata-gex-GRCh38-2024-A"

# ==============================================
# 2. DYNAMIC SAMPLE SHEET GENERATION
# ==============================================
SAMPLES_CSV="$(mktemp)" || { echo "Failed to create temp file"; exit 1; }
echo "sample_id,fastq_dir" > "$SAMPLES_CSV"

process_batch() {
    local dir="$1"
    [ -d "$dir" ] || { echo "Directory $dir not found"; return 1; }
    
    for r1 in "$dir"/*_R1_*.fastq.gz; do
        [ -e "$r1" ] || { echo "No R1 files in $dir"; break; }
        # Extract the sample ID, assuming it's the first part before '_S' in the filename
        sample=$(basename "$r1" | sed 's/_S[0-9]*.*//')
        echo "$sample,$dir" >> "$SAMPLES_CSV" || return 1
    done
}

# Process two FASTQ directories
process_batch "/home/fdann/Desktop/proj/fastq_ngisthlm01223" || exit 1
process_batch "/home/fdann/Desktop/proj/fastq_ngisthlm01379" || exit 1

echo "=== Sample Sheet ==="
cat "$SAMPLES_CSV"
echo "============================="

# ==============================================
# 3. SAMPLE PROCESSING 
# ==============================================
OUTPUT_DIR="/home/fdann/Desktop/proj/results"
mkdir -p "$OUTPUT_DIR" || { echo "Failed to create OUTPUT_DIR"; exit 1; }

echo "=== Generated sample sheet ==="
cat "$SAMPLES_CSV"
echo "============================="

while IFS=, read -r sample_id fastq_dir; do
    [ "$sample_id" == "sample_id" ] && continue

    OUTPUT_LOCATION_1="$OUTPUT_DIR/$sample_id"
    OUTPUT_LOCATION_2="$OUTPUT_DIR/$sample_id"
    # Check if already processed in EITHER location
    if [ -d "$OUTPUT_LOCATION_1" ] || [ -d "$OUTPUT_LOCATION_2" ]; then
        echo "=== Skipping $sample_id (already processed in one of the output directories) ==="
        continue
    fi

    echo "=== Processing $sample_id ==="
    
    cellranger count \
        --id "$sample_id" \
        --transcriptome "$REF_PATH" \
        --fastqs "$fastq_dir" \
        --sample "$sample_id" \
        --localcores "$SLURM_CPUS_PER_TASK" \
        --create-bam false 2>&1 | tee "$OUTPUT_DIR/${sample_id}.log"
    
    if [ -d "$sample_id" ]; then
        mv "$sample_id" "$OUTPUT_DIR/" && \
        echo "Success: $sample_id → $OUTPUT_DIR/"
    else
        echo "ERROR: $sample_id failed (see $OUTPUT_DIR/${sample_id}.log)"
    fi
done < <(tail -n +2 "$SAMPLES_CSV")

rm "$SAMPLES_CSV"
echo "=== Analysis complete ==="
echo "Results in: $OUTPUT_DIR"
echo "Log files:"
ls -lh "$OUTPUT_DIR"/*.log
