#!/usr/bin/env bash
#SBATCH --job-name=bin_annotation
#SBATCH --array=0-99%10
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --account=ofinkel

# Script to run Prodigal and MicrobeAnnotator on COMEBin and SemiBin bins
# One array task per bin directory (sample + binner combination)

set -euo pipefail

# Set default values for variables that config utilities might check
export SAMPLE_INFO_FILE="${SAMPLE_INFO_FILE:-}"
export TREATMENTS_FILE="${TREATMENTS_FILE:-}"
export TOTAL_SAMPLES="${TOTAL_SAMPLES:-0}"

# Source configuration and utilities
if [ -n "${PIPELINE_SCRIPT_DIR:-}" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# ==================== USER SETTINGS ====================

# Base directory containing binning results
BINNING_BASE="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/binning"

# MicrobeAnnotator database location
MICROBEANNOTATOR_DB="/sci/backup/aerez/aerez/moshea/MicrobeAnnotator_DB/"

# Treatment to process (e.g., carR)
TREATMENT="${TREATMENT:-carR}"

# MicrobeAnnotator parallelization settings
PROCESSES=10              # Number of parallel processes (-p)
THREADS_PER_PROCESS=5     # Threads per process (-t)

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
echo "[$(date)] Bin directory: $BIN_DIR"
echo "[$(date)] Protein output directory: $PROTEIN_DIR"

# Count total bins first
total_bins=$(ls -1 "${BIN_DIR}"/*.fa 2>/dev/null | wc -l)
echo "[$(date)] Found $total_bins bin files to process"

if [ $total_bins -eq 0 ]; then
    echo "[$(date)] ERROR: No .fa files found in $BIN_DIR"
    exit 1
fi

# List first few bins
echo "[$(date)] First 5 bins:"
ls -1 "${BIN_DIR}"/*.fa | head -5

bin_count=0
failed_count=0
protein_files=()

# Temporarily disable exit on error for Prodigal loop
set +e

for bin_file in "${BIN_DIR}"/*.fa; do
    if [ ! -f "$bin_file" ]; then
        echo "[$(date)]   Skipping non-file: $bin_file"
        continue
    fi

    base=$(basename "$bin_file")
    stem="${base%.*}"

    out_genes="${PROTEIN_DIR}/${stem}-genes.fna"
    out_prot="${PROTEIN_DIR}/${stem}-proteins.faa"

    echo "[$(date)]   [$((bin_count + failed_count + 1))/$total_bins] Processing bin: $base"

    # Run prodigal and capture output
    if prodigal -i "$bin_file" -d "$out_genes" -a "$out_prot" > "${PROTEIN_DIR}/${stem}-prodigal.log" 2>&1; then
        if [ -f "$out_prot" ] && [ -s "$out_prot" ]; then
            echo "[$(date)]   ✓ Prodigal completed for $base ($(grep -c '^>' "$out_prot" 2>/dev/null || echo 0) proteins)"
            protein_files+=("$out_prot")
            ((bin_count++))
        else
            echo "[$(date)]   ✗ Prodigal ran but no proteins generated for $base"
            ((failed_count++))
        fi
    else
        echo "[$(date)]   ✗ Prodigal failed for $base (check ${PROTEIN_DIR}/${stem}-prodigal.log)"
        ((failed_count++))
    fi
done

# Re-enable exit on error
set -e

echo "[$(date)] Prodigal summary:"
echo "[$(date)]   Successful: $bin_count bins"
echo "[$(date)]   Failed: $failed_count bins"
echo "[$(date)]   Total: $total_bins bins"

if [ ${#protein_files[@]} -eq 0 ]; then
    echo "[$(date)] ERROR: No protein files generated"
    exit 1
fi

# ==================== STEP 2: Run MicrobeAnnotator ====================

echo "[$(date)] ====== Running MicrobeAnnotator ======"
echo "[$(date)] Step 1: Preparing to activate conda environment..."

# Initialize conda if not already done
export CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"
echo "[$(date)] CONDA_BASE: $CONDA_BASE"

if [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    echo "[$(date)] ERROR: conda.sh not found at ${CONDA_BASE}/etc/profile.d/conda.sh"
    exit 1
fi

echo "[$(date)] Step 2: Sourcing conda.sh..."
source "${CONDA_BASE}/etc/profile.d/conda.sh"

echo "[$(date)] Step 3: Activating microbeannotator environment..."
if ! conda activate microbeannotator; then
    echo "[$(date)] ERROR: Failed to activate microbeannotator conda environment"
    echo "[$(date)] Available environments:"
    conda env list
    exit 1
fi

echo "[$(date)] Step 4: Checking if microbeannotator command is available..."
if ! command -v microbeannotator &> /dev/null; then
    echo "[$(date)] ERROR: MicrobeAnnotator not available in conda environment"
    echo "[$(date)] Current conda environment: $CONDA_DEFAULT_ENV"
    echo "[$(date)] PATH: $PATH"
    conda deactivate
    exit 1
fi

echo "[$(date)] Step 5: MicrobeAnnotator found at: $(which microbeannotator)"
echo "[$(date)] Step 6: Verifying protein directory: $PROTEIN_DIR"
if [ ! -d "$PROTEIN_DIR" ]; then
    echo "[$(date)] ERROR: Protein directory not found: $PROTEIN_DIR"
    conda deactivate
    exit 1
fi

echo "[$(date)] Step 7: Changing to protein directory..."
cd "$PROTEIN_DIR" || {
    echo "[$(date)] ERROR: Failed to change to protein directory"
    conda deactivate
    exit 1
}

echo "[$(date)] Step 8: Listing protein files..."
ls -lh *.faa | head -5
echo "[$(date)] Total protein files: $(ls -1 *.faa 2>/dev/null | wc -l)"

# Create annotation output directory
echo "[$(date)] Step 9: Creating annotation output directory: $ANNOTATION_DIR"
mkdir -p "$ANNOTATION_DIR"

# Run MicrobeAnnotator on all protein files
echo "[$(date)] Step 10: Starting MicrobeAnnotator on ${#protein_files[@]} protein files"
echo "[$(date)] Using $PROCESSES processes with $THREADS_PER_PROCESS threads each"
echo "[$(date)] Database: $MICROBEANNOTATOR_DB"
echo "[$(date)] Output directory: $ANNOTATION_DIR"
echo "[$(date)] Current working directory: $(pwd)"
echo "[$(date)] Starting at: $(date)"

set +e  # Temporarily disable exit on error to capture output
microbeannotator \
    -i *.faa \
    -m diamond \
    -d "$MICROBEANNOTATOR_DB" \
    -p "$PROCESSES" \
    -t "$THREADS_PER_PROCESS" \
    -o "$ANNOTATION_DIR" \
    --cluster both 2>&1 | tee "${ANNOTATION_DIR}/microbeannotator.log"

EXIT_CODE=${PIPESTATUS[0]}
set -e  # Re-enable exit on error

echo "[$(date)] MicrobeAnnotator finished with exit code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "[$(date)] ✓ MicrobeAnnotator completed successfully"
    conda deactivate
else
    echo "[$(date)] ✗ MicrobeAnnotator failed with exit code: $EXIT_CODE"
    echo "[$(date)] Check log file: ${ANNOTATION_DIR}/microbeannotator.log"
    conda deactivate
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
