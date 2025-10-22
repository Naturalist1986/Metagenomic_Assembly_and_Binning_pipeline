#!/bin/bash
#SBATCH --job-name=coassembly
#SBATCH --array=0-9%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=512G
#SBATCH --time=7-0

# 01b_coassembly.sh - MetaSPAdes co-assembly stage (combines all samples per treatment)

# CRITICAL: Ensure OUTPUT_DIR is set and normalize it (remove trailing slash)
if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR environment variable not set!" >&2
    exit 1
fi

# Remove trailing slashes
OUTPUT_DIR="${OUTPUT_DIR%/}"
export OUTPUT_DIR

# Set dependent paths
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"

# Source configuration and utilities FIRST (before setting file paths)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# NOW set these variables AFTER sourcing (so they don't get overwritten by empty initializations)
export TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
export SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"

# Debug: Print all paths
echo "========================================" >&2
echo "Environment Check:" >&2
echo "  OUTPUT_DIR: $OUTPUT_DIR" >&2
echo "  WORK_DIR: $WORK_DIR" >&2
echo "  TREATMENTS_FILE: $TREATMENTS_FILE" >&2
echo "  SAMPLE_INFO_FILE: $SAMPLE_INFO_FILE" >&2
echo "========================================" >&2

# Verify critical files exist
if [ ! -f "$TREATMENTS_FILE" ]; then
    echo "ERROR: Treatments file not found at: $TREATMENTS_FILE" >&2
    echo "Listing WORK_DIR contents:" >&2
    ls -la "$WORK_DIR" 2>&1 >&2 || echo "WORK_DIR doesn't exist!" >&2
    exit 1
fi

if [ ! -f "$SAMPLE_INFO_FILE" ]; then
    echo "ERROR: Sample info file not found at: $SAMPLE_INFO_FILE" >&2
    exit 1
fi

# Get treatment from array task ID
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

# Read treatments into array directly from file
mapfile -t TREATMENTS_ARRAY < "$TREATMENTS_FILE"

# Debug treatment array
echo "Treatment array info:" >&2
echo "  Total treatments: ${#TREATMENTS_ARRAY[@]}" >&2
echo "  Array index: $ARRAY_INDEX" >&2
for i in "${!TREATMENTS_ARRAY[@]}"; do
    echo "    [$i]: '${TREATMENTS_ARRAY[$i]}'" >&2
done

if [ $ARRAY_INDEX -ge ${#TREATMENTS_ARRAY[@]} ]; then
    echo "No treatment found for array index $ARRAY_INDEX" >&2
    exit 0
fi

TREATMENT="${TREATMENTS_ARRAY[$ARRAY_INDEX]}"
echo "Selected treatment: '$TREATMENT'" >&2
export TREATMENT

# Initialize
init_conda
create_treatment_dirs "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Co-Assembly for Treatment: $TREATMENT ======"

# Check if co-assembly already completed
if check_treatment_checkpoint "$TREATMENT" "coassembly"; then
    log "Co-assembly already completed for $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main co-assembly function
stage_coassembly() {
    local treatment="$1"
    
    log "Running metaSPAdes co-assembly for treatment: $treatment"
    
    local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
    mkdir -p "$coassembly_dir"
    
    # Check if already assembled
    if [ -f "${coassembly_dir}/contigs.fasta" ] && [ -s "${coassembly_dir}/contigs.fasta" ]; then
        log "Treatment $treatment already co-assembled, skipping..."
        return 0
    fi
    
    # Collect all samples for this treatment
    local samples=($(get_samples_for_treatment "$treatment"))
    
    if [ ${#samples[@]} -eq 0 ]; then
        log "ERROR: No samples found for treatment $treatment"
        return 1
    fi
    
    log "Found ${#samples[@]} samples for co-assembly: ${samples[*]}"
    
    # Collect all input files (prefer validated, fall back to quality_filtering)
    local r1_files=()
    local r2_files=()
    local singleton_files=()
    local total_r1_reads=0
    local total_r2_reads=0
    local total_singleton_reads=0
    
    for sample in "${samples[@]}"; do
        # Check for validated files first (after stage 0.5)
        local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample}"
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample}"
        
        local input1=""
        local input2=""
        local singletons=""
        
        if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
            input1="${validated_dir}/validated_1.fastq.gz"
            input2="${validated_dir}/validated_2.fastq.gz"
            singletons="${validated_dir}/singletons.fastq.gz"
            log "  Sample $sample: using validated files"
        elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
            input1="${quality_dir}/filtered_1.fastq.gz"
            input2="${quality_dir}/filtered_2.fastq.gz"
            singletons="${quality_dir}/singletons.fastq.gz"
            log "  Sample $sample: using quality-filtered files"
        else
            log "WARNING: No input files found for $sample, skipping"
            continue
        fi
        
        # Validate files exist
        if [ ! -f "$input1" ] || [ ! -f "$input2" ]; then
            log "WARNING: Missing input reads for $sample, skipping"
            continue
        fi
        
        # Validate read counts match
        if ! validate_read_counts "$input1" "$input2" "$sample"; then
            log "ERROR: Read count mismatch for $sample, aborting co-assembly"
            return 1
        fi
        
        r1_files+=("$input1")
        r2_files+=("$input2")
        
        # Count reads
        local r1_reads=$(count_reads "$input1")
        local r2_reads=$(count_reads "$input2")
        total_r1_reads=$((total_r1_reads + r1_reads))
        total_r2_reads=$((total_r2_reads + r2_reads))
        
        log "  Sample $sample: R1=$r1_reads, R2=$r2_reads reads"
        
        # Add singletons if available
        if [ -f "$singletons" ] && [ -s "$singletons" ]; then
            singleton_files+=("$singletons")
            local singleton_count=$(count_reads "$singletons")
            total_singleton_reads=$((total_singleton_reads + singleton_count))
            log "  Sample $sample singletons: $singleton_count reads"
        fi
    done
    
    # Verify we have files to assemble
    if [ ${#r1_files[@]} -eq 0 ]; then
        log "ERROR: No valid input files found for treatment $treatment"
        return 1
    fi
    
    log "Total reads for co-assembly:"
    log "  R1: $total_r1_reads"
    log "  R2: $total_r2_reads"
    log "  Singletons: $total_singleton_reads"
    
    # Create concatenated input files
    local concat_r1="${TEMP_DIR}/${treatment}_combined_R1.fastq.gz"
    local concat_r2="${TEMP_DIR}/${treatment}_combined_R2.fastq.gz"
    local concat_singletons="${TEMP_DIR}/${treatment}_combined_singletons.fastq.gz"
    
    log "Concatenating R1 files..."
    cat "${r1_files[@]}" > "$concat_r1"
    
    log "Concatenating R2 files..."
    cat "${r2_files[@]}" > "$concat_r2"
    
    # Concatenate singletons if any exist
    if [ ${#singleton_files[@]} -gt 0 ]; then
        log "Concatenating singleton files..."
        cat "${singleton_files[@]}" > "$concat_singletons"
    fi
    
    # Final validation of concatenated files
    log "Validating concatenated files..."
    if ! validate_read_counts "$concat_r1" "$concat_r2" "$treatment"; then
        log "ERROR: Concatenated read count mismatch for $treatment"
        return 1
    fi
    
    # Activate SPAdes environment
    activate_env spades
    
    # Prepare SPAdes command
    local spades_cmd="spades.py --meta"
    spades_cmd+=" -1 $concat_r1"
    spades_cmd+=" -2 $concat_r2"
    
    # Add concatenated singletons if available
    if [ -f "$concat_singletons" ] && [ -s "$concat_singletons" ]; then
        spades_cmd+=" -s $concat_singletons"
        log "Including $total_singleton_reads singleton reads in co-assembly"
    fi
    
    spades_cmd+=" -o $coassembly_dir"
    spades_cmd+=" -t $SLURM_CPUS_PER_TASK"
    spades_cmd+=" -m 500"
    spades_cmd+=" --tmp-dir ${TEMP_DIR}/spades_${treatment}"
    spades_cmd+=" --only-assembler"
    
    # Create temporary directory for SPAdes
    mkdir -p "${TEMP_DIR}/spades_${treatment}"
    
    # Run SPAdes
    log "Running metaSPAdes co-assembly command:"
    log "$spades_cmd"
    
    eval $spades_cmd 2>&1 | tee "${LOG_DIR}/${treatment}_coassembly.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ] && [ -f "${coassembly_dir}/contigs.fasta" ]; then
        # Create symlink for compatibility
        ln -sf "${coassembly_dir}/contigs.fasta" "${coassembly_dir}/final_assembly.fasta"

        # Log assembly statistics
        log_coassembly_stats "$treatment" "$coassembly_dir"

        # Clean up SPAdes intermediate files
        log "Cleaning up intermediate files..."
        rm -rf "${coassembly_dir}/corrected" "${coassembly_dir}/misc" "${coassembly_dir}/tmp"
        rm -rf "${coassembly_dir}/K21" "${coassembly_dir}/K33" "${coassembly_dir}/K55"

        # Save merged reads for downstream use (bin reassembly, etc.)
        local merged_reads_dir="${coassembly_dir}/merged_reads"
        mkdir -p "$merged_reads_dir"

        log "Saving merged reads for treatment $treatment..."
        mv "$concat_r1" "${merged_reads_dir}/merged_R1.fastq.gz"
        mv "$concat_r2" "${merged_reads_dir}/merged_R2.fastq.gz"

        if [ -f "$concat_singletons" ] && [ -s "$concat_singletons" ]; then
            mv "$concat_singletons" "${merged_reads_dir}/merged_singletons.fastq.gz"
            log "  Saved merged reads: R1, R2, and singletons"
        else
            log "  Saved merged reads: R1 and R2"
        fi

        log "Co-assembly successful for treatment: $treatment"
        return 0
    else
        log "ERROR: metaSPAdes co-assembly failed for $treatment (exit code: $exit_code)"
        rm -rf "$coassembly_dir"
        return 1
    fi
    
    conda deactivate
}

# Log co-assembly statistics
log_coassembly_stats() {
    local treatment="$1"
    local assembly_dir="$2"
    local contigs_file="${assembly_dir}/contigs.fasta"
    
    if [ ! -f "$contigs_file" ]; then
        return 1
    fi
    
    log "Co-assembly statistics for treatment: $treatment"
    
    # Basic statistics
    local num_contigs=$(grep -c '^>' "$contigs_file")
    local total_length=$(grep -v '^>' "$contigs_file" | tr -d '\n' | wc -c)
    
    log "  Total contigs: $num_contigs"
    log "  Total length: $total_length bp"
    
    # Extract contig lengths
    local lengths_file="${TEMP_DIR}/contig_lengths.txt"
    grep '^>' "$contigs_file" | sed 's/.*length_\([0-9]*\).*/\1/' > "$lengths_file"
    
    if [ -s "$lengths_file" ]; then
        local max_length=$(sort -nr "$lengths_file" | head -1)
        local min_length=$(sort -n "$lengths_file" | head -1)
        
        log "  Maximum contig length: $max_length bp"
        log "  Minimum contig length: $min_length bp"
        
        # Calculate N50
        local n50=$(calculate_n50 "$lengths_file")
        log "  N50: $n50 bp"
        
        # Length distribution
        local long_contigs=$(awk '$1 >= 1000' "$lengths_file" | wc -l)
        local very_long_contigs=$(awk '$1 >= 10000' "$lengths_file" | wc -l)
        
        log "  Contigs ≥ 1kb: $long_contigs"
        log "  Contigs ≥ 10kb: $very_long_contigs"
    fi
    
    rm -f "$lengths_file"
}

# Calculate N50
calculate_n50() {
    local lengths_file="$1"
    
    awk '
        {lengths[NR] = $1; total += $1}
        END {
            n = asort(lengths, sorted_lengths, "@val_num_desc")
            target = total / 2
            cumulative = 0
            for (i = 1; i <= n; i++) {
                cumulative += sorted_lengths[i]
                if (cumulative >= target) {
                    print sorted_lengths[i]
                    break
                }
            }
        }
    ' "$lengths_file"
}

# Run the co-assembly stage
if stage_coassembly "$TREATMENT"; then
    # Create checkpoint
    create_treatment_checkpoint "$TREATMENT" "coassembly"
    
    log "====== Co-assembly completed for treatment: $TREATMENT ======"
else
    log "ERROR: Co-assembly stage failed for treatment: $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"