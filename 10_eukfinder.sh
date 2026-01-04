#!/bin/bash
#SBATCH --job-name=eukfinder
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# 10_eukfinder.sh - Run EukFinder on largest bins from binning stage

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get bin info from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Use bins list file from environment variable if set, otherwise use largest bins list
if [ -n "$BINS_LIST_FILE" ]; then
    # Using custom bins list (e.g., all bins)
    BINS_FILE="$BINS_LIST_FILE"
    echo "Using bins list: $BINS_FILE"
else
    # Default to largest bins list
    BINS_FILE="${OUTPUT_DIR}/largest_bins_list.txt"
    echo "Using largest bins list: $BINS_FILE"
fi

if [ ! -f "$BINS_FILE" ]; then
    echo "ERROR: Bins file not found: $BINS_FILE"
    echo "Please run identify_largest_bins.sh or identify_all_bins.sh first"
    exit 1
fi

# Get the bin for this task ID (skip comment lines)
BIN_INFO=$(grep -v "^#" "$BINS_FILE" | sed -n "$((TASK_ID + 1))p")

if [ -z "$BIN_INFO" ]; then
    echo "No bin found for array index $TASK_ID"
    exit 0
fi

# Parse bin information
IFS='|' read -r BIN_PATH TREATMENT SAMPLE BINNER BIN_NAME BIN_SIZE NUM_CONTIGS <<< "$BIN_INFO"

echo "====== Starting EukFinder for bin ======"
echo "Bin: $BIN_NAME"
echo "Sample: $SAMPLE"
echo "Treatment: $TREATMENT"
echo "Binner: $BINNER"
echo "Size: $BIN_SIZE bp"
echo "Contigs: $NUM_CONTIGS"
echo ""

# Initialize
init_conda

# Create output prefix for this specific bin
OUTPUT_PREFIX="${SAMPLE}_${BINNER}_$(basename "$BIN_NAME" .fa)"

# Setup output directories - one directory per bin
EUKFINDER_OUTPUT_DIR="${OUTPUT_DIR}/eukfinder/${TREATMENT}/${OUTPUT_PREFIX}"
mkdir -p "$EUKFINDER_OUTPUT_DIR"

# Create log directory
mkdir -p "${LOG_DIR}/eukfinder/${TREATMENT}"

# Setup temporary directory for this bin
TEMP_DIR=$(mktemp -d -p "${WORK_DIR}" eukfinder_${SAMPLE}_${BINNER}_XXXXXX)
log "Created temporary directory: $TEMP_DIR"

# Check if already processed
CHECKPOINT_FILE="${EUKFINDER_OUTPUT_DIR}/.${OUTPUT_PREFIX}.done"

if [ -f "$CHECKPOINT_FILE" ]; then
    log "EukFinder already completed for $OUTPUT_PREFIX"
    rm -rf "$TEMP_DIR"
    exit 0
fi

# Function to run EukFinder on a bin
run_eukfinder() {
    local bin_file="$1"
    local output_prefix="$2"
    local work_dir="$3"

    log "Running EukFinder on $bin_file"

    # Activate EukFinder conda environment
    activate_env eukfinder

    # Check if EukFinder is available
    if ! command -v eukfinder &> /dev/null; then
        log "ERROR: eukfinder command not available in conda environment"
        conda deactivate
        return 1
    fi

    # Verify database paths are configured
    if [ -z "$EUKFINDER_CENTRIFUGE_DB" ] || [ "$EUKFINDER_CENTRIFUGE_DB" == "/path/to/centrifuge/database" ]; then
        log "ERROR: EUKFINDER_CENTRIFUGE_DB not configured in 00_config_utilities.sh"
        conda deactivate
        return 1
    fi

    if [ -z "$EUKFINDER_PLAST_DB" ] || [ "$EUKFINDER_PLAST_DB" == "/path/to/plast/database" ]; then
        log "ERROR: EUKFINDER_PLAST_DB not configured in 00_config_utilities.sh"
        conda deactivate
        return 1
    fi

    if [ -z "$EUKFINDER_PLAST_ID_MAP" ] || [ "$EUKFINDER_PLAST_ID_MAP" == "/path/to/plast/id_map" ]; then
        log "ERROR: EUKFINDER_PLAST_ID_MAP not configured in 00_config_utilities.sh"
        conda deactivate
        return 1
    fi

    log "Using Centrifuge database: $EUKFINDER_CENTRIFUGE_DB"
    log "Using PLAST database: $EUKFINDER_PLAST_DB"
    log "Using PLAST ID map: $EUKFINDER_PLAST_ID_MAP"

    # Change to working directory (EukFinder creates outputs in current directory)
    cd "$work_dir" || {
        log "ERROR: Could not change to working directory: $work_dir"
        conda deactivate
        return 1
    }

    # EukFinder parameters
    local threads=${EUKFINDER_THREADS:-${SLURM_CPUS_PER_TASK:-48}}
    local chunks=${EUKFINDER_CHUNKS:-6}
    local taxonomy_update=${EUKFINDER_TAXONOMY_UPDATE:-False}
    local evalue=${EUKFINDER_EVALUE:-0.01}
    local pid=${EUKFINDER_PID:-60}
    local cov=${EUKFINDER_COV:-30}
    local mhlen=${EUKFINDER_MHLEN:-100}

    log "EukFinder parameters:"
    log "  Threads: $threads"
    log "  Chunks: $chunks"
    log "  Taxonomy update: $taxonomy_update"
    log "  E-value: $evalue"
    log "  PID threshold: $pid"
    log "  Coverage threshold: $cov"
    log "  Minimum hit length: $mhlen"

    # Run EukFinder
    log "Running: eukfinder long_seqs -l $bin_file -o $output_prefix --mhlen $mhlen --cdb $EUKFINDER_CENTRIFUGE_DB -n $threads -z $chunks -t $taxonomy_update -p $EUKFINDER_PLAST_DB -m $EUKFINDER_PLAST_ID_MAP -e $evalue --pid $pid --cov $cov"

    eukfinder long_seqs \
        -l "$bin_file" \
        -o "$output_prefix" \
        --mhlen "$mhlen" \
        --cdb "$EUKFINDER_CENTRIFUGE_DB" \
        -n "$threads" \
        -z "$chunks" \
        -t "$taxonomy_update" \
        -p "$EUKFINDER_PLAST_DB" \
        -m "$EUKFINDER_PLAST_ID_MAP" \
        -e "$evalue" \
        --pid "$pid" \
        --cov "$cov" \
        2>&1 | tee "${LOG_DIR}/eukfinder/${TREATMENT}/${output_prefix}_eukfinder.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    if [ $exit_code -eq 0 ]; then
        # Verify that results were actually created
        if [ -d "Eukfinder_results" ] || [ -d "Intermediate_data" ]; then
            log "EukFinder completed successfully"
            return 0
        else
            log "ERROR: EukFinder exited with code 0 but produced no results"
            log "This usually indicates a database configuration issue"
            log "Check the log file at: ${LOG_DIR}/eukfinder/${TREATMENT}/${output_prefix}_eukfinder.log"
            return 1
        fi
    else
        log "ERROR: EukFinder failed with exit code: $exit_code"
        return 1
    fi
}

# Copy results from working directory to final output
copy_eukfinder_results() {
    local work_dir="$1"
    local output_dir="$2"
    local output_prefix="$3"

    log "Copying EukFinder results to final output directory..."

    # Copy Eukfinder_results directory
    if [ -d "${work_dir}/Eukfinder_results" ]; then
        mkdir -p "${output_dir}/Eukfinder_results"
        cp -r "${work_dir}/Eukfinder_results/"* "${output_dir}/Eukfinder_results/" 2>/dev/null || true
        log "Copied Eukfinder_results"
    fi

    # Copy Intermediate_data directory
    if [ -d "${work_dir}/Intermediate_data" ]; then
        mkdir -p "${output_dir}/Intermediate_data"
        cp -r "${work_dir}/Intermediate_data/"* "${output_dir}/Intermediate_data/" 2>/dev/null || true
        log "Copied Intermediate_data"
    fi

    # Create a summary file
    create_eukfinder_summary "$output_dir" "$output_prefix"
}

# Create EukFinder summary
create_eukfinder_summary() {
    local output_dir="$1"
    local output_prefix="$2"
    local summary_file="${output_dir}/eukfinder_summary.txt"

    cat >> "$summary_file" << EOF

EukFinder Analysis: $output_prefix
=====================================

Date: $(date)
Bin: $BIN_NAME
Sample: $SAMPLE
Treatment: $TREATMENT
Binner: $BINNER
Original bin size: $BIN_SIZE bp
Number of contigs: $NUM_CONTIGS

Results:
EOF

    # Check if EukFinder created a summary_table.txt
    local results_dir="${output_dir}/Eukfinder_results"
    local summary_table="${results_dir}/summary_table.txt"

    if [ -f "$summary_table" ]; then
        # Parse the summary_table.txt file that EukFinder creates
        log "Parsing summary_table.txt"
        while IFS=$'\t' read -r group num_seq total_size; do
            # Skip the header line
            if [[ "$group" == "Group" ]]; then
                continue
            fi
            # Format and add to summary
            echo "  ${group}: ${num_seq} sequences, ${total_size} bp" >> "$summary_file"
        done < "$summary_table"
    elif [ -d "$results_dir" ]; then
        # Fallback: Count sequences in individual FASTA files
        log "Counting sequences in FASTA files"
        for category in Euk Unk EUnk Bact Arch Misc; do
            local cat_file="${results_dir}/${category}.fasta"
            if [ -f "$cat_file" ]; then
                local count=$(grep -c "^>" "$cat_file" 2>/dev/null || echo "0")
                local size=$(grep -v "^>" "$cat_file" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
                echo "  ${category}: $count sequences, $size bp" >> "$summary_file"
            else
                echo "  ${category}: No sequences" >> "$summary_file"
            fi
        done
    else
        echo "  No results directory found" >> "$summary_file"
    fi

    echo "" >> "$summary_file"
    log "Summary updated: $summary_file"
}

# Main execution
log "====== Starting EukFinder for $OUTPUT_PREFIX ======"

if run_eukfinder "$BIN_PATH" "$OUTPUT_PREFIX" "$TEMP_DIR"; then
    # Copy results to final output directory
    copy_eukfinder_results "$TEMP_DIR" "$EUKFINDER_OUTPUT_DIR" "$OUTPUT_PREFIX"

    # Create checkpoint
    touch "$CHECKPOINT_FILE"

    log "====== EukFinder completed for $OUTPUT_PREFIX ======"
else
    log "ERROR: EukFinder failed for $OUTPUT_PREFIX"
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Cleanup
rm -rf "$TEMP_DIR"
log "Cleaned up temporary directory"

exit 0
