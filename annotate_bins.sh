#!/usr/bin/env bash
#SBATCH --job-name=bin_annotation
#SBATCH --array=0-99%10
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/slurm/bin_annotation_%A_%a.log
#SBATCH --error=logs/slurm/bin_annotation_%A_%a.err
#SBATCH --account=ofinkel

# Script to run Prodigal and MicrobeAnnotator on COMEBin and SemiBin bins
# One array task per bin directory (sample + binner combination)

set -euo pipefail

# ==================== USER SETTINGS ====================

# Base directory containing binning results
BINNING_BASE="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/binning"

# MicrobeAnnotator database location
MICROBEANNOTATOR_DB="/sci/backup/aerez/aerez/moshea/MicrobeAnnotator_DB/"

# Treatment to process (e.g., carR)
TREATMENT="${TREATMENT:-carR}"

# Number of threads for MicrobeAnnotator
THREADS=8

# =======================================================

echo "[$(date)] ====== Starting Bin Annotation for Array Task ${SLURM_ARRAY_TASK_ID} ======"

# Build array of bin directories to process
# Format: sample_name|binner|bin_directory
BIN_DIRS=()

# Find all sample directories for this treatment
for sample_dir in "${BINNING_BASE}/${TREATMENT}"/*; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi

    sample_name=$(basename "$sample_dir")

    # Check for COMEBin bins
    comebin_dir="${sample_dir}/comebin/comebin_res/comebin_res_bins"
    if [ ! -d "$comebin_dir" ]; then
        comebin_dir="${sample_dir}/comebin/comebin_res_bins"
    fi

    if [ -d "$comebin_dir" ] && [ "$(ls -A ${comebin_dir}/*.fa 2>/dev/null)" ]; then
        BIN_DIRS+=("${sample_name}|comebin|${comebin_dir}")
    fi

    # Check for SemiBin bins
    semibin_dir="${sample_dir}/semibin/output_bins"
    if [ -d "$semibin_dir" ] && [ "$(ls -A ${semibin_dir}/*.fa 2>/dev/null)" ]; then
        BIN_DIRS+=("${sample_name}|semibin|${semibin_dir}")
    fi
done

# Check if array task ID is valid
if [ ${SLURM_ARRAY_TASK_ID} -ge ${#BIN_DIRS[@]} ]; then
    echo "[$(date)] Array task ${SLURM_ARRAY_TASK_ID} has no bin directory to process (only ${#BIN_DIRS[@]} directories found)"
    exit 0
fi

# Get bin directory for this task
IFS='|' read -r SAMPLE_NAME BINNER BIN_DIR <<< "${BIN_DIRS[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date)] Sample: $SAMPLE_NAME"
echo "[$(date)] Binner: $BINNER"
echo "[$(date)] Bin directory: $BIN_DIR"

# Set up output directory
OUTPUT_DIR="${BINNING_BASE}/${TREATMENT}/${SAMPLE_NAME}/annotations/${BINNER}"
PROTEIN_DIR="${OUTPUT_DIR}/proteins"
ANNOTATION_DIR="${OUTPUT_DIR}/microbeannotator"

mkdir -p "$PROTEIN_DIR"
mkdir -p "$ANNOTATION_DIR"

# ==================== STEP 1: Run Prodigal ====================

echo "[$(date)] ====== Running Prodigal on bins ======"

bin_count=0
protein_files=()

for bin_file in "${BIN_DIR}"/*.fa; do
    if [ ! -f "$bin_file" ]; then
        continue
    fi

    base=$(basename "$bin_file")
    stem="${base%.*}"

    out_genes="${PROTEIN_DIR}/${stem}-genes.fna"
    out_prot="${PROTEIN_DIR}/${stem}-proteins.faa"

    echo "[$(date)]   Processing bin: $base"

    if prodigal -i "$bin_file" -d "$out_genes" -a "$out_prot" 2>&1; then
        echo "[$(date)]   ✓ Prodigal completed for $base"
        protein_files+=("$out_prot")
        ((bin_count++))
    else
        echo "[$(date)]   ✗ Prodigal failed for $base"
    fi
done

echo "[$(date)] Prodigal completed for $bin_count bins"

if [ ${#protein_files[@]} -eq 0 ]; then
    echo "[$(date)] ERROR: No protein files generated"
    exit 1
fi

# ==================== STEP 2: Run MicrobeAnnotator ====================

echo "[$(date)] ====== Running MicrobeAnnotator ======"

# Change to protein directory for MicrobeAnnotator
cd "$PROTEIN_DIR"

# Run MicrobeAnnotator on all protein files
echo "[$(date)] Running MicrobeAnnotator on ${#protein_files[@]} protein files"

if microbeannotator \
    -i *.faa \
    -m diamond \
    -d "$MICROBEANNOTATOR_DB" \
    -p "$THREADS" \
    -t 5 \
    -o "$ANNOTATION_DIR" \
    --cluster both 2>&1 | tee "${ANNOTATION_DIR}/microbeannotator.log"; then

    echo "[$(date)] ✓ MicrobeAnnotator completed successfully"
else
    echo "[$(date)] ✗ MicrobeAnnotator failed"
    exit 1
fi

# ==================== Summary ====================

echo "[$(date)] ====== Annotation Summary ======"
echo "[$(date)] Sample: $SAMPLE_NAME"
echo "[$(date)] Binner: $BINNER"
echo "[$(date)] Bins processed: $bin_count"
echo "[$(date)] Protein files: ${PROTEIN_DIR}"
echo "[$(date)] Annotations: ${ANNOTATION_DIR}"
echo "[$(date)] ====== Annotation completed successfully ======"
