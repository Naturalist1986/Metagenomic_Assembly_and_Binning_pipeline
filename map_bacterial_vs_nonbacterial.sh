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

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get treatment from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Get directories
BACTERIAL_DIR="${OUTPUT_DIR}/bacterial_sequences"
NONBACTERIAL_DIR="${OUTPUT_DIR}/nonbacterial_sequences"
MAPPING_DIR="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping"

mkdir -p "$MAPPING_DIR"
mkdir -p "${LOG_DIR}/bacterial_vs_nonbacterial_mapping"

# Get list of treatments (use bacterial files as reference)
TREATMENT_FILES=("${BACTERIAL_DIR}"/*_bacterial.fasta)

if [ ${#TREATMENT_FILES[@]} -eq 0 ] || [ ! -f "${TREATMENT_FILES[0]}" ]; then
    echo "ERROR: No bacterial sequence files found"
    echo "Please run collect_bacterial_sequences.sh first"
    exit 0
fi

# Get treatment for this task
if [ $TASK_ID -ge ${#TREATMENT_FILES[@]} ]; then
    echo "No treatment for array index $TASK_ID"
    exit 0
fi

BACTERIAL_FASTA="${TREATMENT_FILES[$TASK_ID]}"
TREATMENT=$(basename "$BACTERIAL_FASTA" _bacterial.fasta)
NONBACTERIAL_FASTA="${NONBACTERIAL_DIR}/${TREATMENT}_nonbacterial.fasta"

echo "====== Starting Bacterial vs Non-Bacterial Mapping for $TREATMENT ======"
echo "Bacterial sequences: $BACTERIAL_FASTA"
echo "Non-bacterial sequences: $NONBACTERIAL_FASTA"
echo ""

# Verify non-bacterial file exists
if [ ! -f "$NONBACTERIAL_FASTA" ]; then
    echo "ERROR: Non-bacterial sequence file not found: $NONBACTERIAL_FASTA"
    echo "Please run collect_nonbacterial_sequences.sh first"
    exit 1
fi

# Initialize
init_conda

# Setup output directory for this treatment
TREATMENT_MAPPING_DIR="${MAPPING_DIR}/${TREATMENT}"
mkdir -p "$TREATMENT_MAPPING_DIR"

# Check if already completed
CHECKPOINT_FILE="${TREATMENT_MAPPING_DIR}/.mapping_complete"
if [ -f "$CHECKPOINT_FILE" ]; then
    echo "Mapping already completed for $TREATMENT"
    exit 0
fi

# Find reads for this treatment
READS_R1=""
READS_R2=""

# Try multiple possible locations
if [ -f "${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/merged_R1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/merged_R1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/merged_R2.fastq.gz"
elif [ -f "${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/filtered_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/filtered_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/filtered_2.fastq.gz"
elif [ -f "${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_2.fastq.gz"
elif [ -f "${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_2.fastq.gz"
fi

if [ -z "$READS_R1" ] || [ ! -f "$READS_R1" ]; then
    echo "ERROR: Could not find read files for treatment $TREATMENT"
    echo "Searched locations:"
    echo "  ${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/"
    echo "  ${OUTPUT_DIR}/assembly/${TREATMENT}/"
    echo "  ${OUTPUT_DIR}/quality_filtering/${TREATMENT}/"
    exit 1
fi

echo "Found reads:"
echo "  R1: $READS_R1"
echo "  R2: $READS_R2"
echo ""

# Count total reads
TOTAL_READS=$(zcat "$READS_R1" | echo $((`wc -l`/4)))
echo "Total read pairs: $TOTAL_READS"
echo ""

# Activate BBMap environment
activate_env bbtools

# Function to map reads to reference
map_to_reference() {
    local reference="$1"
    local output_prefix="$2"
    local label="$3"

    echo "Mapping to $label sequences..."

    # Create BAM output
    local bam_output="${TREATMENT_MAPPING_DIR}/${output_prefix}.bam"
    local stats_output="${TREATMENT_MAPPING_DIR}/${output_prefix}_stats.txt"

    bbmap.sh \
        in1="$READS_R1" \
        in2="$READS_R2" \
        ref="$reference" \
        out="$bam_output" \
        bamscript="${TREATMENT_MAPPING_DIR}/${output_prefix}_script.sh" \
        statsfile="$stats_output" \
        covstats="${TREATMENT_MAPPING_DIR}/${output_prefix}_covstats.txt" \
        rpkm="${TREATMENT_MAPPING_DIR}/${output_prefix}_rpkm.txt" \
        threads=${SLURM_CPUS_PER_TASK:-16} \
        minid=0.95 \
        ambiguous=random \
        2>&1 | tee "${LOG_DIR}/bacterial_vs_nonbacterial_mapping/${TREATMENT}_${output_prefix}.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: Mapping to $label failed"
        return 1
    fi

    # Parse mapped reads count from stats file
    local mapped_reads=$(grep "mapped:" "$stats_output" | awk '{print $2}')
    echo "$mapped_reads"
    return 0
}

# Map to bacterial sequences
BACTERIAL_MAPPED=$(map_to_reference "$BACTERIAL_FASTA" "bacterial" "bacterial")
if [ $? -ne 0 ]; then
    deactivate_env
    exit 1
fi

echo "Bacterial mapped reads: $BACTERIAL_MAPPED"
echo ""

# Map to non-bacterial sequences
NONBACTERIAL_MAPPED=$(map_to_reference "$NONBACTERIAL_FASTA" "nonbacterial" "non-bacterial")
if [ $? -ne 0 ]; then
    deactivate_env
    exit 1
fi

echo "Non-bacterial mapped reads: $NONBACTERIAL_MAPPED"
echo ""

deactivate_env

# Calculate percentages
BACTERIAL_PCT=$(awk -v mapped="$BACTERIAL_MAPPED" -v total="$TOTAL_READS" 'BEGIN {printf "%.2f", (mapped/total)*100}')
NONBACTERIAL_PCT=$(awk -v mapped="$NONBACTERIAL_MAPPED" -v total="$TOTAL_READS" 'BEGIN {printf "%.2f", (mapped/total)*100}')

# Calculate unmapped (reads that didn't map to either)
TOTAL_MAPPED=$((BACTERIAL_MAPPED + NONBACTERIAL_MAPPED))
UNMAPPED=$((TOTAL_READS - TOTAL_MAPPED))
UNMAPPED_PCT=$(awk -v unmapped="$UNMAPPED" -v total="$TOTAL_READS" 'BEGIN {printf "%.2f", (unmapped/total)*100}')

# Create summary file
SUMMARY_FILE="${TREATMENT_MAPPING_DIR}/mapping_summary.txt"
cat > "$SUMMARY_FILE" << EOF
Bacterial vs Non-Bacterial Mapping Summary
Treatment: $TREATMENT
===========================================

Total Reads: $TOTAL_READS

Bacterial (Bact):
  Mapped reads: $BACTERIAL_MAPPED
  Percentage: ${BACTERIAL_PCT}%

Non-Bacterial (Arch + Euk + Unk + EUnk + Misc):
  Mapped reads: $NONBACTERIAL_MAPPED
  Percentage: ${NONBACTERIAL_PCT}%

Unmapped:
  Unmapped reads: $UNMAPPED
  Percentage: ${UNMAPPED_PCT}%

Date: $(date)
EOF

echo ""
echo "Summary:"
echo "--------"
cat "$SUMMARY_FILE"

# Create checkpoint
touch "$CHECKPOINT_FILE"

echo ""
echo "====== Mapping completed for $TREATMENT ======"
exit 0
