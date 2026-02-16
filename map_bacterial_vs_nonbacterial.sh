#!/bin/bash
#SBATCH --job-name=bact_map
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel

# map_bacterial_vs_nonbacterial.sh - Map reads to bacterial vs non-bacterial sequences
# Maps to EACH BINNER SEPARATELY to avoid duplicate contig counting

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

log "Starting bacterial vs non-bacterial mapping (Array index: $ARRAY_INDEX)"

# Define directories
EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder_output"
MAPPING_DIR="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping"

# Get list of treatments from EukFinder output
TREATMENTS=($(ls -1 "$EUKFINDER_DIR" 2>/dev/null | sort))

if [ ${#TREATMENTS[@]} -eq 0 ]; then
    log "ERROR: No treatments found in $EUKFINDER_DIR"
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

# Check if already completed
SUMMARY_FILE="${TREATMENT_DIR}/mapping_summary.txt"
if [ -f "$SUMMARY_FILE" ]; then
    log "Mapping already completed for $TREATMENT (found summary file)"
    exit 0
fi

# Identify unique binners for this treatment
log "Identifying binners from EukFinder results..."

EUKFINDER_TREATMENT_DIR="${EUKFINDER_DIR}/${TREATMENT}"

if [ ! -d "$EUKFINDER_TREATMENT_DIR" ]; then
    log "ERROR: EukFinder directory not found: $EUKFINDER_TREATMENT_DIR"
    exit 1
fi

# Extract binner names from directory names (format: sample_binner_binname or treatment_binner_binname)
BINNERS=($(ls -1 "$EUKFINDER_TREATMENT_DIR" | grep -oE '_[^_]+_' | tr -d '_' | sort -u))

if [ ${#BINNERS[@]} -eq 0 ]; then
    log "ERROR: No binners identified for treatment $TREATMENT"
    exit 1
fi

log "Found ${#BINNERS[@]} binners: ${BINNERS[*]}"

# Find reads for this treatment
log "Locating input reads..."

# Try multiple possible locations
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
        # Look for R1/R2 files (try multiple naming patterns)
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

# Map to each binner separately
log "Mapping reads to each binner's sequences..."

declare -A binner_bacterial_reads
declare -A binner_nonbacterial_reads
declare -A binner_unmapped_reads
total_reads=0

for BINNER in "${BINNERS[@]}"; do
    log "Processing binner: $BINNER"

    # Collect sequences for this binner only
    BINNER_BACTERIAL_FA="${TREATMENT_DIR}/${BINNER}_bacterial.fasta"
    BINNER_NONBACTERIAL_FA="${TREATMENT_DIR}/${BINNER}_nonbacterial.fasta"

    > "$BINNER_BACTERIAL_FA"
    > "$BINNER_NONBACTERIAL_FA"

    # Extract sequences from bins belonging to this binner
    for bin_dir in "${EUKFINDER_TREATMENT_DIR}"/*_${BINNER}_*/; do
        if [ ! -d "$bin_dir" ]; then
            continue
        fi

        results_dir="${bin_dir}/Eukfinder_results"
        if [ ! -d "$results_dir" ]; then
            continue
        fi

        bin_name=$(basename "$bin_dir")

        # Collect bacterial sequences (Bact only)
        if [ -f "${results_dir}/Bact.fasta" ] && [ -s "${results_dir}/Bact.fasta" ]; then
            awk -v bin="$bin_name" '/^>/ {print $0"|bin:"bin; next} {print}' \
                "${results_dir}/Bact.fasta" >> "$BINNER_BACTERIAL_FA"
        fi

        # Collect non-bacterial sequences (Arch + Euk + Unk + EUnk + Misc)
        for category in Arch Euk Unk EUnk Misc; do
            if [ -f "${results_dir}/${category}.fasta" ] && [ -s "${results_dir}/${category}.fasta" ]; then
                awk -v bin="$bin_name" -v cat="$category" '/^>/ {print $0"|bin:"bin"|cat:"cat; next} {print}' \
                    "${results_dir}/${category}.fasta" >> "$BINNER_NONBACTERIAL_FA"
            fi
        done
    done

    # Check if we have sequences for this binner
    bact_seqs=$(grep -c "^>" "$BINNER_BACTERIAL_FA" 2>/dev/null || echo "0")
    nonbact_seqs=$(grep -c "^>" "$BINNER_NONBACTERIAL_FA" 2>/dev/null || echo "0")

    if [ "$bact_seqs" -eq 0 ] && [ "$nonbact_seqs" -eq 0 ]; then
        log "  No sequences found for binner $BINNER - skipping"
        rm -f "$BINNER_BACTERIAL_FA" "$BINNER_NONBACTERIAL_FA"
        continue
    fi

    log "  Bacterial sequences: $bact_seqs"
    log "  Non-bacterial sequences: $nonbact_seqs"

    # Map to bacterial sequences
    bact_mapped=0
    if [ "$bact_seqs" -gt 0 ]; then
        log "  Mapping to bacterial sequences..."

        bbmap.sh \
            in1="$R1_PATH" \
            in2="$R2_PATH" \
            ref="$BINNER_BACTERIAL_FA" \
            out="${TREATMENT_DIR}/${BINNER}_bacterial.bam" \
            statsfile="${TREATMENT_DIR}/${BINNER}_bacterial_stats.txt" \
            covstats="${TREATMENT_DIR}/${BINNER}_bacterial_covstats.txt" \
            minid=0.95 \
            ambiguous=random \
            threads=${SLURM_CPUS_PER_TASK:-16} \
            -Xmx${SLURM_MEM_PER_NODE:-64}g \
            2>&1 | tee "${TREATMENT_DIR}/${BINNER}_bacterial_bbmap.log"

        bact_mapped=$(grep "mapped:" "${TREATMENT_DIR}/${BINNER}_bacterial_stats.txt" | \
            head -1 | awk '{print $3}' | sed 's/,//g' | tr -d '\n\r' || echo "0")
        bact_mapped=${bact_mapped:-0}
    fi

    # Map to non-bacterial sequences
    nonbact_mapped=0
    if [ "$nonbact_seqs" -gt 0 ]; then
        log "  Mapping to non-bacterial sequences..."

        bbmap.sh \
            in1="$R1_PATH" \
            in2="$R2_PATH" \
            ref="$BINNER_NONBACTERIAL_FA" \
            out="${TREATMENT_DIR}/${BINNER}_nonbacterial.bam" \
            statsfile="${TREATMENT_DIR}/${BINNER}_nonbacterial_stats.txt" \
            covstats="${TREATMENT_DIR}/${BINNER}_nonbacterial_covstats.txt" \
            minid=0.95 \
            ambiguous=random \
            threads=${SLURM_CPUS_PER_TASK:-16} \
            -Xmx${SLURM_MEM_PER_NODE:-64}g \
            2>&1 | tee "${TREATMENT_DIR}/${BINNER}_nonbacterial_bbmap.log"

        nonbact_mapped=$(grep "mapped:" "${TREATMENT_DIR}/${BINNER}_nonbacterial_stats.txt" | \
            head -1 | awk '{print $3}' | sed 's/,//g' | tr -d '\n\r' || echo "0")
        nonbact_mapped=${nonbact_mapped:-0}
    fi

    # Get total reads (from first binner)
    if [ $total_reads -eq 0 ]; then
        if [ "$bact_seqs" -gt 0 ]; then
            total_reads=$(grep "reads:" "${TREATMENT_DIR}/${BINNER}_bacterial_stats.txt" | \
                head -1 | awk '{print $2}' | sed 's/,//g' | tr -d '\n\r')
        elif [ "$nonbact_seqs" -gt 0 ]; then
            total_reads=$(grep "reads:" "${TREATMENT_DIR}/${BINNER}_nonbacterial_stats.txt" | \
                head -1 | awk '{print $2}' | sed 's/,//g' | tr -d '\n\r')
        fi
        total_reads=${total_reads:-0}
    fi

    # Calculate unmapped for this binner
    unmapped=$((total_reads - bact_mapped - nonbact_mapped))

    # Store results
    binner_bacterial_reads[$BINNER]=$bact_mapped
    binner_nonbacterial_reads[$BINNER]=$nonbact_mapped
    binner_unmapped_reads[$BINNER]=$unmapped

    log "  Results for $BINNER:"
    log "    Bacterial: $bact_mapped reads"
    log "    Non-bacterial: $nonbact_mapped reads"
    log "    Unmapped: $unmapped reads"
done

deactivate_env

# Calculate averages across binners
log "Calculating average results across binners..."

num_binners=${#binner_bacterial_reads[@]}

if [ $num_binners -eq 0 ]; then
    log "ERROR: No binners had mappable sequences"
    exit 1
fi

avg_bact=0
avg_nonbact=0
avg_unmapped=0

for binner in "${!binner_bacterial_reads[@]}"; do
    avg_bact=$((avg_bact + binner_bacterial_reads[$binner]))
    avg_nonbact=$((avg_nonbact + binner_nonbacterial_reads[$binner]))
    avg_unmapped=$((avg_unmapped + binner_unmapped_reads[$binner]))
done

avg_bact=$((avg_bact / num_binners))
avg_nonbact=$((avg_nonbact / num_binners))
avg_unmapped=$((avg_unmapped / num_binners))

# Calculate percentages
bact_pct=$(awk "BEGIN {printf \"%.2f\", ($avg_bact / $total_reads) * 100}")
nonbact_pct=$(awk "BEGIN {printf \"%.2f\", ($avg_nonbact / $total_reads) * 100}")
unmapped_pct=$(awk "BEGIN {printf \"%.2f\", ($avg_unmapped / $total_reads) * 100}")

# Create summary file
cat > "$SUMMARY_FILE" << EOF
Bacterial vs Non-Bacterial Mapping Summary
Treatment: $TREATMENT
Generated: $(date)
===========================================

Input Reads:
  R1: $R1_PATH
  R2: $R2_PATH
  Total Reads: $total_reads

Binners Analyzed: $num_binners
  ${BINNERS[*]}

Results (Averaged Across Binners):
-----------------------------------

Bacterial (Bact):
  Mapped reads: $avg_bact
  Percentage: ${bact_pct}%

Non-Bacterial (Arch + Euk + Unk + EUnk + Misc):
  Mapped reads: $avg_nonbact
  Percentage: ${nonbact_pct}%

Unmapped:
  Unmapped reads: $avg_unmapped
  Percentage: ${unmapped_pct}%

Per-Binner Results:
-------------------
EOF

for binner in "${!binner_bacterial_reads[@]}"; do
    bact=${binner_bacterial_reads[$binner]}
    nonbact=${binner_nonbacterial_reads[$binner]}
    unmap=${binner_unmapped_reads[$binner]}

    bact_p=$(awk "BEGIN {printf \"%.2f\", ($bact / $total_reads) * 100}")
    nonbact_p=$(awk "BEGIN {printf \"%.2f\", ($nonbact / $total_reads) * 100}")
    unmap_p=$(awk "BEGIN {printf \"%.2f\", ($unmap / $total_reads) * 100}")

    cat >> "$SUMMARY_FILE" << EOF

$binner:
  Bacterial: $bact (${bact_p}%)
  Non-bacterial: $nonbact (${nonbact_p}%)
  Unmapped: $unmap (${unmap_p}%)
EOF
done

log "Mapping complete for treatment $TREATMENT"
log "Summary saved to: $SUMMARY_FILE"
exit 0
