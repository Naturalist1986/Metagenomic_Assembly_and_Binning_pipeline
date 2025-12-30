#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=3-0

# 01_assembly.sh - MetaSPAdes assembly stage with singletons

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get sample info from array task ID
SAMPLE_INFO=$(get_sample_info_by_index $SLURM_ARRAY_TASK_ID)
if [ -z "$SAMPLE_INFO" ]; then
    log "No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 0
fi

# Parse sample information
IFS='|' read -r SAMPLE_NAME TREATMENT _ _ <<< "$SAMPLE_INFO"
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Assembly for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "assembly"; then
    log "Assembly already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_assembly() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Running metaSPAdes assembly for $sample_name ($treatment)"
    
    local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    mkdir -p "$assembly_dir"
    
    # Check if already assembled
    if [ -f "${assembly_dir}/contigs.fasta" ] && [ -s "${assembly_dir}/contigs.fasta" ]; then
        log "Sample $sample_name already assembled, skipping..."
        return 0
    fi
    
    # Check for validated files first (after stage 0.5), fall back to quality_filtering
    local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    local input_dir=""
    
    if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
        input_dir="$validated_dir"
        local input1="${input_dir}/validated_1.fastq.gz"
        local input2="${input_dir}/validated_2.fastq.gz"
        local singletons="${input_dir}/singletons.fastq.gz"
        log "Using validated files from: $input_dir"
    elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
        input_dir="$quality_dir"
        local input1="${input_dir}/filtered_1.fastq.gz"
        local input2="${input_dir}/filtered_2.fastq.gz"
        local singletons="${input_dir}/singletons.fastq.gz"
        log "Using quality-filtered files from: $input_dir"
    else
        log "ERROR: No input files found for $sample_name"
        log "  Checked validated: $validated_dir"
        log "  Checked quality: $quality_dir"
        return 1
    fi
    
    # Verify input files exist
    if [ ! -f "$input1" ] || [ ! -f "$input2" ]; then
        log "ERROR: Missing input reads for $sample_name"
        log "  Expected: $input1 and $input2"
        return 1
    fi
    
    # Validate read counts match before assembly
    log "Validating input file synchronization..."
    if ! validate_read_counts "$input1" "$input2" "$sample_name"; then
        log "ERROR: Read count validation failed for $sample_name"
        log "Input files are not synchronized - this should not happen after validation stage!"
        return 1
    fi
    
    # Activate SPAdes environment
    activate_env spades
    
    # Log input statistics
    log "Input files for $sample_name:"
    local r1_reads=$(count_reads "$input1")
    local r2_reads=$(count_reads "$input2")
    log "  R1: $input1 ($r1_reads reads)"
    log "  R2: $input2 ($r2_reads reads)"
    
    # Check for singletons
    local singleton_reads=0
    if [ -f "$singletons" ] && [ -s "$singletons" ]; then
        singleton_reads=$(count_reads "$singletons")
        log "  Singletons: $singletons ($singleton_reads reads)"
    else
        log "  No singleton reads found or file is empty"
    fi
    
    # Prepare SPAdes command
    local spades_cmd="spades.py --meta"
    spades_cmd+=" -1 $input1"
    spades_cmd+=" -2 $input2"
    
    # Add singletons if available and non-empty
    if [ -f "$singletons" ] && [ -s "$singletons" ] && [ $singleton_reads -gt 0 ]; then
        spades_cmd+=" -s $singletons"
        log "Including $singleton_reads singleton reads in assembly"
    fi
    
    spades_cmd+=" -o $assembly_dir"
    spades_cmd+=" -t $ASSEMBLY_THREADS"
    spades_cmd+=" -m $ASSEMBLY_MEMORY"
    spades_cmd+=" --tmp-dir ${TEMP_DIR}/spades_${sample_name}"
    spades_cmd+=" --only-assembler"  # Skip read error correction to save time
    
    # Create temporary directory for SPAdes
    mkdir -p "${TEMP_DIR}/spades_${sample_name}"
    
    # Run SPAdes
    log "Running metaSPAdes command:"
    log "$spades_cmd"
    
    eval $spades_cmd 2>&1 | tee "${LOG_DIR}/${treatment}/${sample_name}_assembly.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ] && [ -f "${assembly_dir}/contigs.fasta" ]; then
        # Create symlink for compatibility
        ln -sf "${assembly_dir}/contigs.fasta" "${assembly_dir}/final_assembly.fasta"
        
        # Log assembly statistics
        log_assembly_stats "$sample_name" "$assembly_dir"
        
        # Clean up SPAdes intermediate files to save space
        log "Cleaning up intermediate files..."
        rm -rf "${assembly_dir}/corrected" "${assembly_dir}/misc" "${assembly_dir}/tmp"
        rm -rf "${assembly_dir}/K21" "${assembly_dir}/K33" "${assembly_dir}/K55"
        
        log "Assembly successful for $sample_name"
        return 0
    else
        log "ERROR: metaSPAdes failed for $sample_name (exit code: $exit_code)"
        rm -rf "$assembly_dir"
        return 1
    fi
    
    conda deactivate
}

# Log assembly statistics
log_assembly_stats() {
    local sample_name="$1"
    local assembly_dir="$2"
    local contigs_file="${assembly_dir}/contigs.fasta"

    if [ ! -f "$contigs_file" ]; then
        return 1
    fi

    log "Assembly statistics for $sample_name:"

    # Basic statistics
    local num_contigs=$(grep -c '^>' "$contigs_file")
    local total_length=$(grep -v '^>' "$contigs_file" | tr -d '\n' | wc -c)

    log "  Total contigs: $num_contigs"
    log "  Total length: $total_length bp"

    # Extract contig lengths from headers (SPAdes format)
    local lengths_file="${TEMP_DIR}/contig_lengths.txt"
    grep '^>' "$contigs_file" | sed 's/.*length_\([0-9]*\).*/\1/' > "$lengths_file"

    if [ -s "$lengths_file" ]; then
        local max_length=$(sort -nr "$lengths_file" | head -1)
        local min_length=$(sort -n "$lengths_file" | head -1)
        local median_length=$(sort -n "$lengths_file" | awk '{a[NR]=$1} END {print (NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}')

        log "  Maximum contig length: $max_length bp"
        log "  Minimum contig length: $min_length bp"
        log "  Median contig length: $median_length bp"

        # Calculate N50
        local n50=$(calculate_n50 "$lengths_file")
        log "  N50: $n50 bp"

        # Length distribution
        local long_contigs=$(awk '$1 >= 1000' "$lengths_file" | wc -l)
        local very_long_contigs=$(awk '$1 >= 10000' "$lengths_file" | wc -l)

        log "  Contigs ≥ 1kb: $long_contigs"
        log "  Contigs ≥ 10kb: $very_long_contigs"
    fi

    # GC content
    local gc_content=$(calculate_gc_content "$contigs_file")
    log "  GC content: ${gc_content}%"

    # Calculate assembly success rate
    log "  Calculating assembly success rate..."
    local input_dir=""
    local validated_dir="${OUTPUT_DIR}/validated/${TREATMENT}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/${sample_name}"

    # Get input files (same logic as in stage_assembly)
    if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
        local input_r1="${validated_dir}/validated_1.fastq.gz"
        local input_r2="${validated_dir}/validated_2.fastq.gz"
    elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
        local input_r1="${quality_dir}/filtered_1.fastq.gz"
        local input_r2="${quality_dir}/filtered_2.fastq.gz"
    else
        log "  WARNING: Could not find input files for assembly success rate calculation"
        rm -f "$lengths_file"
        return 0
    fi

    # Calculate success rate
    local assembly_success_rate=$(calculate_assembly_success_rate \
        "$contigs_file" \
        "$input_r1" \
        "$input_r2" \
        "$sample_name" \
        "${assembly_dir}/assembly_mapping")

    log "  Assembly success rate: ${assembly_success_rate}%"

    # Save success rate to file for later reference
    echo "$assembly_success_rate" > "${assembly_dir}/assembly_success_rate.txt"

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

# Calculate GC content
calculate_gc_content() {
    local fasta_file="$1"
    
    grep -v '^>' "$fasta_file" | \
    awk '
        {
            sequence = sequence $0
        }
        END {
            gc_count = gsub(/[GCgc]/, "", sequence)
            total_count = length(sequence)
            if (total_count > 0) {
                gc_percentage = (gc_count * 100) / total_count
                printf "%.2f", gc_percentage
            } else {
                print "0.00"
            }
        }
    '
}

# Generate assembly statistics report
generate_assembly_stats() {
    local sample_name="$1"
    local treatment="$2"
    local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    local stats_file="${assembly_dir}/assembly_statistics.txt"
    
    if [ ! -f "${assembly_dir}/contigs.fasta" ]; then
        return 1
    fi
    
    echo "Assembly Statistics for $sample_name" > "$stats_file"
    echo "===================================" >> "$stats_file"
    echo "" >> "$stats_file"
    echo "Date: $(date)" >> "$stats_file"
    echo "Sample: $sample_name" >> "$stats_file"
    echo "Treatment: $treatment" >> "$stats_file"
    echo "" >> "$stats_file"
    
    # Detailed statistics
    local contigs_file="${assembly_dir}/contigs.fasta"
    local num_contigs=$(grep -c '^>' "$contigs_file")
    local total_length=$(grep -v '^>' "$contigs_file" | tr -d '\n' | wc -c)
    
    echo "Basic Statistics:" >> "$stats_file"
    echo "  Total contigs: $num_contigs" >> "$stats_file"
    echo "  Total length: $total_length bp" >> "$stats_file"
    
    # Contig length statistics
    local lengths_file="${TEMP_DIR}/contig_lengths_report.txt"
    grep '^>' "$contigs_file" | sed 's/.*length_\([0-9]*\).*/\1/' > "$lengths_file"
    
    if [ -s "$lengths_file" ]; then
        local max_length=$(sort -nr "$lengths_file" | head -1)
        local min_length=$(sort -n "$lengths_file" | head -1)
        local n50=$(calculate_n50 "$lengths_file")
        
        echo "  Maximum contig length: $max_length bp" >> "$stats_file"
        echo "  Minimum contig length: $min_length bp" >> "$stats_file"
        echo "  N50: $n50 bp" >> "$stats_file"
        
        # Length distribution
        echo "" >> "$stats_file"
        echo "Length Distribution:" >> "$stats_file"
        echo "  Contigs ≥ 1kb: $(awk '$1 >= 1000' "$lengths_file" | wc -l)" >> "$stats_file"
        echo "  Contigs ≥ 5kb: $(awk '$1 >= 5000' "$lengths_file" | wc -l)" >> "$stats_file"
        echo "  Contigs ≥ 10kb: $(awk '$1 >= 10000' "$lengths_file" | wc -l)" >> "$stats_file"
        echo "  Contigs ≥ 50kb: $(awk '$1 >= 50000' "$lengths_file" | wc -l)" >> "$stats_file"
    fi
    
    # GC content
    local gc_content=$(calculate_gc_content "$contigs_file")
    echo "" >> "$stats_file"
    echo "Composition:" >> "$stats_file"
    echo "  GC content: ${gc_content}%" >> "$stats_file"

    # Add assembly success rate if available
    if [ -f "${assembly_dir}/assembly_success_rate.txt" ]; then
        local success_rate=$(cat "${assembly_dir}/assembly_success_rate.txt")
        echo "" >> "$stats_file"
        echo "Assembly Quality:" >> "$stats_file"
        echo "  Assembly success rate: ${success_rate}%" >> "$stats_file"
        echo "  (Percentage of input reads that mapped to assembled contigs)" >> "$stats_file"
    fi

    rm -f "$lengths_file"

    log "Assembly statistics report created: $stats_file"
}

# Run the assembly stage
if stage_assembly "$SAMPLE_NAME" "$TREATMENT"; then
    # Generate detailed statistics
    generate_assembly_stats "$SAMPLE_NAME" "$TREATMENT"
    
    # Create checkpoint
    create_sample_checkpoint "$SAMPLE_NAME" "assembly"
    
    log "====== Assembly completed for $SAMPLE_NAME ======"
else
    log "ERROR: Assembly stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"