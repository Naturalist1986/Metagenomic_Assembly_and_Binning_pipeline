#!/bin/bash
#SBATCH --job-name=quality_filtering
#SBATCH --account=$SLURM_ACCOUNT
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00

# 00_quality_filtering.sh - Quality filtering using Trimmomatic

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get sample info from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
SAMPLE_INFO=$(get_sample_info_by_index "$TASK_ID")

if [ -z "$SAMPLE_INFO" ]; then
    log "No sample found for array index $TASK_ID"
    exit 0
fi

# Parse sample information
IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"

# Validate parsed variables
if [ -z "$SAMPLE_NAME" ] || [ -z "$TREATMENT" ]; then
    log "ERROR: Failed to parse sample information properly"
    exit 1
fi

# Export variables
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Quality Filtering for ${SAMPLE_NAME} (${TREATMENT}) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "quality_filtering"; then
    log "Quality filtering already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_quality_filtering() {
    local sample_name="$1"
    local treatment="$2"
    local r1_path="$3"
    local r2_path="$4"
    
    # Check for merged files first (from stage -1), fall back to original paths
    local merged_dir="${OUTPUT_DIR}/merged_lanes/${treatment}/${sample_name}"
    if [ -f "${merged_dir}/merged_R1.fq.gz" ]; then
        r1_path="${merged_dir}/merged_R1.fq.gz"
        r2_path="${merged_dir}/merged_R2.fq.gz"
        log "Using merged lane files from: $merged_dir"
    else
        log "Using original input files (single lane or no merge stage)"
    fi
    
    log "Running quality filtering for $sample_name ($treatment)"
    log "  R1: $r1_path"
    log "  R2: $r2_path"
    
    # Create output directories
    local output_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    mkdir -p "$output_dir"
    
    # Check if already processed
    if [ -f "${output_dir}/filtered_1.fastq.gz" ] && [ -f "${output_dir}/filtered_2.fastq.gz" ]; then
        log "Sample $sample_name already processed, skipping..."
        return 0
    fi
    
    # Verify input files exist
    if [ ! -f "$r1_path" ] || [ ! -f "$r2_path" ]; then
        log "ERROR: Input files not found for $sample_name"
        log "  R1: $r1_path"
        log "  R2: $r2_path"
        return 1
    fi
    
    # Validate input file synchronization
    log "Validating input file synchronization..."
    if ! validate_read_counts "$r1_path" "$r2_path" "$sample_name"; then
        log "ERROR: Input files are not synchronized!"
        return 1
    fi
    
    # Activate Trimmomatic environment
    activate_env metawrap-env
    
    # Set up Trimmomatic parameters
    local adapters="${TRIMMOMATIC_ADAPTERS:-TruSeq3-PE-2.fa}"
    local leading="${TRIMMOMATIC_LEADING:-3}"
    local trailing="${TRIMMOMATIC_TRAILING:-3}"
    local slidingwindow="${TRIMMOMATIC_SLIDINGWINDOW:-4:15}"
    local minlen="${TRIMMOMATIC_MINLEN:-36}"
    
    # Run Trimmomatic
    log "Running Trimmomatic with parameters:"
    log "  Adapters: $adapters"
    log "  Leading: $leading"
    log "  Trailing: $trailing"
    log "  Sliding window: $slidingwindow"
    log "  Min length: $minlen"
    
    trimmomatic PE \
        -threads ${SLURM_CPUS_PER_TASK:-4} \
        -phred33 \
        "$r1_path" "$r2_path" \
        "${output_dir}/filtered_1.fastq.gz" \
        "${output_dir}/unpaired_1.fastq.gz" \
        "${output_dir}/filtered_2.fastq.gz" \
        "${output_dir}/unpaired_2.fastq.gz" \
        ILLUMINACLIP:${TRIMMOMATIC_DB}/${adapters}:2:30:10 \
        LEADING:${leading} \
        TRAILING:${trailing} \
        SLIDINGWINDOW:${slidingwindow} \
        MINLEN:${minlen} \
        2>&1 | tee "${LOG_DIR}/${treatment}/${sample_name}_trimmomatic.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        # Validate Trimmomatic output
        log "Validating Trimmomatic output..."
        if ! validate_read_counts "${output_dir}/filtered_1.fastq.gz" \
                                  "${output_dir}/filtered_2.fastq.gz" \
                                  "$sample_name"; then
            log "ERROR: Trimmomatic produced unequal read counts"
            log "This should not happen - Trimmomatic maintains pair synchronization"
            return 1
        fi
        
        # Combine unpaired reads into singletons
        log "Combining unpaired reads into singletons..."
        cat "${output_dir}/unpaired_1.fastq.gz" "${output_dir}/unpaired_2.fastq.gz" > \
            "${output_dir}/singletons.fastq.gz"
        
        # Remove individual unpaired files to save space
        rm -f "${output_dir}/unpaired_1.fastq.gz" "${output_dir}/unpaired_2.fastq.gz"
        
        # Log statistics
        log_filtering_stats "$sample_name" "$output_dir" "$r1_path" "$r2_path"
        
        log "Quality filtering completed successfully for $sample_name"
        return 0
    else
        log "ERROR: Trimmomatic failed for $sample_name (exit code: $exit_code)"
        return 1
    fi
    
    conda deactivate
}

# Log filtering statistics
log_filtering_stats() {
    local sample_name="$1"
    local output_dir="$2"
    local r1_path="$3"
    local r2_path="$4"
    
    log "Quality filtering statistics for $sample_name:"
    
    # Count reads using the count_reads function
    local input_r1_reads=$(count_reads "$r1_path")
    local input_r2_reads=$(count_reads "$r2_path")
    local output_r1_reads=$(count_reads "${output_dir}/filtered_1.fastq.gz")
    local output_r2_reads=$(count_reads "${output_dir}/filtered_2.fastq.gz")
    local singleton_reads=$(count_reads "${output_dir}/singletons.fastq.gz")
    
    log "  Input R1 reads: $input_r1_reads"
    log "  Input R2 reads: $input_r2_reads"
    log "  Output paired R1 reads: $output_r1_reads"
    log "  Output paired R2 reads: $output_r2_reads"
    log "  Singleton reads: $singleton_reads"
    
    # Calculate retention rate
    if [ "$input_r1_reads" -gt 0 ]; then
        local retention_rate=$(echo "scale=2; ($output_r1_reads * 100) / $input_r1_reads" | bc -l 2>/dev/null || echo "N/A")
        log "  Paired read retention rate: ${retention_rate}%"
    fi
}

# Run quality filtering
if stage_quality_filtering "$SAMPLE_NAME" "$TREATMENT" "$R1_PATH" "$R2_PATH"; then
    create_sample_checkpoint "$SAMPLE_NAME" "quality_filtering"
    log "====== Quality filtering completed for $SAMPLE_NAME ======"
else
    log "ERROR: Quality filtering failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"