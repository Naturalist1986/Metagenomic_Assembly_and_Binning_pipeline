#!/bin/bash
#SBATCH --job-name=validate_repair
#SBATCH --array=0-99%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00

# 00b_validate_repair.sh - Validate and repair paired-end reads after QC

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
IFS='|' read -r SAMPLE_NAME TREATMENT _ _ <<< "$SAMPLE_INFO"
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Validation & Repair for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "validate_repair"; then
    log "Validation & repair already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main validation and repair function
validate_and_repair() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Validating paired-end reads for $sample_name ($treatment)"
    
    # Locate quality-filtered files
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    local input_r1="${quality_dir}/filtered_1.fastq.gz"
    local input_r2="${quality_dir}/filtered_2.fastq.gz"
    local input_singletons="${quality_dir}/singletons.fastq.gz"
    
    # Check if input files exist
    if [ ! -f "$input_r1" ] || [ ! -f "$input_r2" ]; then
        log "ERROR: Quality-filtered reads not found for $sample_name"
        log "  Expected: $input_r1 and $input_r2"
        return 1
    fi
    
    # Count reads
    log "Counting reads in quality-filtered files..."
    local r1_count=$(count_reads "$input_r1")
    local r2_count=$(count_reads "$input_r2")
    
    log "  R1 reads: $r1_count"
    log "  R2 reads: $r2_count"
    
    # Check if files are synchronized
    if [ "$r1_count" -eq "$r2_count" ]; then
        log "✅ Files are synchronized, no repair needed"
        
        # Create validated directory with symlinks
        local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
        mkdir -p "$validated_dir"
        
        ln -sf "$(realpath "$input_r1")" "${validated_dir}/validated_1.fastq.gz"
        ln -sf "$(realpath "$input_r2")" "${validated_dir}/validated_2.fastq.gz"
        
        if [ -f "$input_singletons" ] && [ -s "$input_singletons" ]; then
            ln -sf "$(realpath "$input_singletons")" "${validated_dir}/singletons.fastq.gz"
        fi
        
        log "Created symlinks in: $validated_dir"
        return 0
    fi
    
    # Files are not synchronized - repair needed
    log "⚠️  Files are NOT synchronized!"
    log "  R1 has $r1_count reads"
    log "  R2 has $r2_count reads"
    log "  Difference: $((r1_count > r2_count ? r1_count - r2_count : r2_count - r1_count)) reads"
    log ""
    log "Attempting repair with BBMap repair.sh..."
    
    # Setup output directory
    local repair_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
    mkdir -p "$repair_dir"
    
    local output_r1="${repair_dir}/validated_1.fastq.gz"
    local output_r2="${repair_dir}/validated_2.fastq.gz"
    local repair_singletons="${repair_dir}/repair_singletons.fastq.gz"
    
    # Run repair
    if ! repair_paired_reads "$input_r1" "$input_r2" "$output_r1" "$output_r2" "$repair_singletons" "$sample_name"; then
        log "ERROR: Repair failed for $sample_name"
        return 1
    fi
    
    # Count repaired reads
    local repaired_r1=$(count_reads "$output_r1")
    local repaired_r2=$(count_reads "$output_r2")
    local repaired_singletons=$(count_reads "$repair_singletons")
    
    log "Repair statistics:"
    log "  Input paired R1: $r1_count"
    log "  Input paired R2: $r2_count"
    log "  Output paired: $repaired_r1"
    log "  New singletons: $repaired_singletons"
    log "  Total reads recovered: $((repaired_r1 + repaired_singletons))"
    
    # Combine with original singletons if they exist
    if [ -f "$input_singletons" ] && [ -s "$input_singletons" ]; then
        local original_singletons=$(count_reads "$input_singletons")
        log "  Original singletons: $original_singletons"
        
        log "Combining original and repair singletons..."
        cat "$input_singletons" "$repair_singletons" > "${repair_dir}/singletons.fastq.gz"
        rm -f "$repair_singletons"
        
        local total_singletons=$(count_reads "${repair_dir}/singletons.fastq.gz")
        log "  Combined singletons: $total_singletons"
    else
        mv "$repair_singletons" "${repair_dir}/singletons.fastq.gz"
    fi
    
    log "✅ Repair completed successfully for $sample_name"
    return 0
}

# Run validation and repair
if validate_and_repair "$SAMPLE_NAME" "$TREATMENT"; then
    create_sample_checkpoint "$SAMPLE_NAME" "validate_repair"
    log "====== Validation & repair completed for $SAMPLE_NAME ======"
else
    log "ERROR: Validation & repair failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
