#!/bin/bash
#SBATCH --job-name=nonprok_map
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel

# map_to_nonprokaryotic.sh - Map reads to non-prokaryotic sequences

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get treatment from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Get list of treatments
NONPROK_DIR="${OUTPUT_DIR}/nonprokaryotic_sequences"
MAPPING_DIR="${OUTPUT_DIR}/nonprokaryotic_mapping"

mkdir -p "$MAPPING_DIR"
mkdir -p "${LOG_DIR}/nonprokaryotic_mapping"

# Get list of non-prokaryotic fasta files
TREATMENT_FILES=("${NONPROK_DIR}"/*_nonprokaryotic.fasta)

if [ ${#TREATMENT_FILES[@]} -eq 0 ] || [ ! -f "${TREATMENT_FILES[0]}" ]; then
    echo "ERROR: No non-prokaryotic sequence files found"
    echo "Please run collect_nonprokaryotic_sequences.sh first"
    exit 0
fi

# Get treatment for this task
if [ $TASK_ID -ge ${#TREATMENT_FILES[@]} ]; then
    echo "No treatment for array index $TASK_ID"
    exit 0
fi

NONPROK_FASTA="${TREATMENT_FILES[$TASK_ID]}"
TREATMENT=$(basename "$NONPROK_FASTA" _nonprokaryotic.fasta)

echo "====== Starting Read Mapping for $TREATMENT ======"
echo "Non-prokaryotic sequences: $NONPROK_FASTA"
echo ""

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
# Try multiple possible locations for coassembly reads
READS_R1=""
READS_R2=""

# Location 1: In assembly directory (coassembly mode)
if [ -f "${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/assembly/${TREATMENT}/filtered_2.fastq.gz"
# Location 2: In quality_filtering directory
elif [ -f "${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/filtered_2.fastq.gz"
# Location 3: Try coassembly_reads directory
elif [ -f "${OUTPUT_DIR}/coassembly_reads/${TREATMENT}_1.fastq.gz" ]; then
    READS_R1="${OUTPUT_DIR}/coassembly_reads/${TREATMENT}_1.fastq.gz"
    READS_R2="${OUTPUT_DIR}/coassembly_reads/${TREATMENT}_2.fastq.gz"
fi

if [ -z "$READS_R1" ] || [ ! -f "$READS_R1" ]; then
    echo "ERROR: Could not find read files for treatment $TREATMENT"
    echo "Searched in:"
    echo "  ${OUTPUT_DIR}/assembly/${TREATMENT}/"
    echo "  ${OUTPUT_DIR}/quality_filtering/${TREATMENT}/"
    echo "  ${OUTPUT_DIR}/coassembly_reads/"
    exit 1
fi

echo "Found read files:"
echo "  R1: $READS_R1"
echo "  R2: $READS_R2"
echo ""

# Build bowtie2 index
INDEX_DIR="${TREATMENT_MAPPING_DIR}/bowtie2_index"
mkdir -p "$INDEX_DIR"
INDEX_BASE="${INDEX_DIR}/${TREATMENT}_nonprok"

if [ ! -f "${INDEX_BASE}.1.bt2" ]; then
    echo "Building bowtie2 index..."
    activate_env metawrap-env

    bowtie2-build \
        --threads ${SLURM_CPUS_PER_TASK:-16} \
        "$NONPROK_FASTA" \
        "$INDEX_BASE" \
        2>&1 | tee "${LOG_DIR}/nonprokaryotic_mapping/${TREATMENT}_index.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "ERROR: Failed to build bowtie2 index"
        conda deactivate
        exit 1
    fi

    conda deactivate
    echo "Index built successfully"
    echo ""
else
    echo "Bowtie2 index already exists"
    echo ""
fi

# Map reads with bowtie2
echo "Mapping reads to non-prokaryotic sequences..."
activate_env metawrap-env

SAM_FILE="${TREATMENT_MAPPING_DIR}/${TREATMENT}.sam"
BAM_FILE="${TREATMENT_MAPPING_DIR}/${TREATMENT}.bam"
SORTED_BAM="${TREATMENT_MAPPING_DIR}/${TREATMENT}.sorted.bam"

bowtie2 \
    -x "$INDEX_BASE" \
    -1 "$READS_R1" \
    -2 "$READS_R2" \
    -S "$SAM_FILE" \
    -p ${SLURM_CPUS_PER_TASK:-16} \
    --very-sensitive \
    --no-unal \
    2>&1 | tee "${LOG_DIR}/nonprokaryotic_mapping/${TREATMENT}_mapping.log"

MAPPING_EXIT=$?

if [ $MAPPING_EXIT -ne 0 ]; then
    echo "ERROR: Bowtie2 mapping failed"
    conda deactivate
    exit 1
fi

echo ""
echo "Converting SAM to BAM and sorting..."

# Convert to BAM and sort
samtools view -@ ${SLURM_CPUS_PER_TASK:-16} -bS "$SAM_FILE" > "$BAM_FILE"
samtools sort -@ ${SLURM_CPUS_PER_TASK:-16} "$BAM_FILE" -o "$SORTED_BAM"
samtools index "$SORTED_BAM"

# Clean up intermediate files
rm -f "$SAM_FILE" "$BAM_FILE"

echo "BAM file created and sorted: $SORTED_BAM"
echo ""

# Calculate mapping statistics
echo "Calculating mapping statistics..."

STATS_FILE="${TREATMENT_MAPPING_DIR}/mapping_stats.txt"

# Get total reads from the fastq file
echo "Counting total reads..."
TOTAL_READS=$(zcat "$READS_R1" | wc -l)
TOTAL_READS=$((TOTAL_READS / 4))  # Convert lines to read pairs

# Get mapped reads
MAPPED_READS=$(samtools view -c -F 4 "$SORTED_BAM")

# Calculate percentage
PERCENT_MAPPED=$(awk "BEGIN {printf \"%.4f\", ($MAPPED_READS / $TOTAL_READS) * 100}")

# Create statistics file
cat > "$STATS_FILE" << EOF
Mapping Statistics for $TREATMENT
==================================
Date: $(date)

Input Reads:
  Total read pairs: $TOTAL_READS
  Forward reads: $READS_R1
  Reverse reads: $READS_R2

Target Sequences:
  Non-prokaryotic FASTA: $NONPROK_FASTA
  Total sequences: $(grep -c "^>" "$NONPROK_FASTA")

Mapping Results:
  Mapped read pairs: $MAPPED_READS
  Percentage mapped: $PERCENT_MAPPED%

Output Files:
  Sorted BAM: $SORTED_BAM
  BAM index: ${SORTED_BAM}.bai
EOF

cat "$STATS_FILE"

# Also create a simple TSV for easy parsing
echo -e "${TREATMENT}\t${TOTAL_READS}\t${MAPPED_READS}\t${PERCENT_MAPPED}" > "${TREATMENT_MAPPING_DIR}/mapping_summary.tsv"

conda deactivate

# Create checkpoint
touch "$CHECKPOINT_FILE"

echo ""
echo "====== Mapping completed for $TREATMENT ======"

exit 0
