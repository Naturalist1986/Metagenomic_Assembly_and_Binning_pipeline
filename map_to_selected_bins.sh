#!/bin/bash
#SBATCH --job-name=mag_map
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel

# map_to_selected_bins.sh - Map reads to selected/refined bins (MAGs)
# Maps to concatenated selected bins per treatment to measure MAG coverage

set -euo pipefail

# Source configuration and utilities
if [ -n "${PIPELINE_SCRIPT_DIR:-}" ]; then
    SCRIPT_DIR="$PIPELINE_SCRIPT_DIR"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Function for logging
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Get array task ID
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

log "Starting selected bins (MAG) mapping (Array index: $ARRAY_INDEX)"

# Calculate Java memory allocation
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    MEM_GB=$((SLURM_MEM_PER_NODE / 1024))
    JAVA_MEM=$((MEM_GB * 9 / 10))
else
    JAVA_MEM=96
fi
log "Using ${JAVA_MEM}G memory for Java processes"

# Define directories
SELECTED_BINS_DIR="${OUTPUT_DIR}/selected_bins"
MAPPING_DIR="${OUTPUT_DIR}/selected_bins_mapping"

# Get list of treatments from selected_bins directory
TREATMENTS=($(ls -1 "$SELECTED_BINS_DIR" 2>/dev/null | sort))

if [ ${#TREATMENTS[@]} -eq 0 ]; then
    log "ERROR: No treatments found in $SELECTED_BINS_DIR"
    exit 1
fi

# Get treatment for this array index
if [ $ARRAY_INDEX -ge ${#TREATMENTS[@]} ]; then
    log "ERROR: Array index $ARRAY_INDEX exceeds treatment count ${#TREATMENTS[@]}"
    exit 1
fi

TREATMENT="${TREATMENTS[$ARRAY_INDEX]}"

log "Processing treatment: $TREATMENT"

# Create output directory
TREATMENT_DIR="${MAPPING_DIR}/${TREATMENT}"
mkdir -p "$TREATMENT_DIR"

# Clean up any stale BBMap reference indexes from previous runs
if [ -d "${SCRIPT_DIR}/ref" ]; then
    log "Cleaning up stale BBMap reference directory (ignoring NFS errors)"
    rm -rf "${SCRIPT_DIR}/ref" 2>/dev/null || true
fi

# Check if already completed
SUMMARY_FILE="${TREATMENT_DIR}/mag_mapping_summary.txt"
if [ -f "$SUMMARY_FILE" ]; then
    log "MAG mapping already completed for $TREATMENT (found summary file)"
    exit 0
fi

# Find and concatenate selected bins
log "Collecting selected bins for treatment $TREATMENT..."

SELECTED_DIR="${SELECTED_BINS_DIR}/${TREATMENT}"

if [ ! -d "$SELECTED_DIR" ]; then
    log "ERROR: Selected bins directory not found: $SELECTED_DIR"
    exit 1
fi

CONCAT_REF="${TREATMENT_DIR}/selected_bins_concat.fasta"

# Concatenate all .fa files from the selected bins directory
# Handles both coassembly (*.fa directly) and individual (sample_subdir/*.fa) layouts
> "$CONCAT_REF"
bin_count=0

# First try direct .fa files (coassembly mode)
for bin_fa in "${SELECTED_DIR}"/*.fa; do
    if [ -f "$bin_fa" ]; then
        cat "$bin_fa" >> "$CONCAT_REF"
        bin_count=$((bin_count + 1))
    fi
done

# If none found, try sample subdirectories (individual assembly mode)
if [ $bin_count -eq 0 ]; then
    for sample_dir in "${SELECTED_DIR}"/*/; do
        if [ -d "$sample_dir" ]; then
            for bin_fa in "${sample_dir}"*.fa; do
                if [ -f "$bin_fa" ]; then
                    cat "$bin_fa" >> "$CONCAT_REF"
                    bin_count=$((bin_count + 1))
                fi
            done
        fi
    done
fi

if [ $bin_count -eq 0 ]; then
    log "ERROR: No .fa files found in $SELECTED_DIR or its subdirectories"
    exit 1
fi

log "Concatenated $bin_count selected bins into reference"

# Count sequences in concatenated reference
ref_seqs=$(grep -c "^>" "$CONCAT_REF" 2>/dev/null || echo 0)
log "Total contigs in reference: $ref_seqs"

# Find reads for this treatment
log "Locating input reads..."

READ_LOCATIONS=(
    "${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads"
    "${OUTPUT_DIR}/assembly/${TREATMENT}"
    "${OUTPUT_DIR}/quality_filtering/${TREATMENT}"
    "${OUTPUT_DIR}/validate_repair/${TREATMENT}"
    "${INPUT_DIR}/${TREATMENT}"
)

R1_PATH=""
R2_PATH=""

for location in "${READ_LOCATIONS[@]}"; do
    if [ -d "$location" ]; then
        r1_candidates=$(find "$location" -maxdepth 1 \( -name "*R1*.fastq.gz" -o -name "*_1.fastq.gz" -o -name "filtered_1.fastq.gz" -o -name "merged_R1.fastq.gz" \) 2>/dev/null | head -1)
        r2_candidates=$(find "$location" -maxdepth 1 \( -name "*R2*.fastq.gz" -o -name "*_2.fastq.gz" -o -name "filtered_2.fastq.gz" -o -name "merged_R2.fastq.gz" \) 2>/dev/null | head -1)

        if [ -n "$r1_candidates" ] && [ -n "$r2_candidates" ] && [ -f "$r1_candidates" ] && [ -f "$r2_candidates" ]; then
            R1_PATH="$r1_candidates"
            R2_PATH="$r2_candidates"
            break
        fi
    fi
done

if [ -z "$R1_PATH" ] || [ -z "$R2_PATH" ]; then
    log "ERROR: Could not find R1/R2 reads for treatment $TREATMENT"
    log "Searched locations: ${READ_LOCATIONS[*]}"
    exit 1
fi

log "Found reads:"
log "  R1: $R1_PATH"
log "  R2: $R2_PATH"

# Activate BBMap environment
activate_env bbmap

# Map reads to selected bins
log "Mapping reads to selected bins..."

bbmap.sh \
    in1="$R1_PATH" \
    in2="$R2_PATH" \
    ref="$CONCAT_REF" \
    out="${TREATMENT_DIR}/selected_bins.bam" \
    statsfile="${TREATMENT_DIR}/selected_bins_stats.txt" \
    covstats="${TREATMENT_DIR}/selected_bins_covstats.txt" \
    minid=0.95 \
    ambiguous=random \
    nodisk=true \
    threads=${SLURM_CPUS_PER_TASK:-32} \
    -Xmx${JAVA_MEM}g \
    2>&1 | tee "${TREATMENT_DIR}/selected_bins_bbmap.log"

deactivate_env

# Extract mapping statistics
STATS_FILE="${TREATMENT_DIR}/selected_bins_stats.txt"

# Sum both R1 and R2 mapped reads
mag_mapped=$(grep "^mapped:" "$STATS_FILE" | \
    awk '{gsub(/,/,"",$3); sum+=$3} END {print sum+0}')
mag_mapped=${mag_mapped:-0}

# Get total reads
total_reads=$(grep "^Reads Used:" "$STATS_FILE" | \
    awk '{gsub(/,/,"",$3); print $3}' | tr -d '\n\r')
total_reads=${total_reads:-0}

# Calculate percentage
if [ "$total_reads" -gt 0 ]; then
    mag_pct=$(awk "BEGIN {printf \"%.2f\", ($mag_mapped / $total_reads) * 100}")
else
    mag_pct="0.00"
fi

# Create summary file
cat > "$SUMMARY_FILE" << EOF
Selected Bins (MAG) Mapping Summary
Treatment: $TREATMENT
Generated: $(date)
===========================================

Input Reads:
  R1: $R1_PATH
  R2: $R2_PATH
  Total Reads: $total_reads

Selected Bins: $bin_count
Total Contigs: $ref_seqs

Results:
-----------------------------------

MAG-Mapped:
  Mapped reads: $mag_mapped
  Percentage: ${mag_pct}%
EOF

log "MAG mapping complete for treatment $TREATMENT"
log "  Total reads: $total_reads"
log "  MAG-mapped: $mag_mapped (${mag_pct}%)"
log "Summary saved to: $SUMMARY_FILE"
exit 0
