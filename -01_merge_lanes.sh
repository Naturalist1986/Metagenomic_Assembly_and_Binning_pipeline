#!/bin/bash
#SBATCH --job-name=merge_lanes
#SBATCH --array=0-99%20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# -01_merge_lanes.sh - Detect and merge multiple sequencing lanes per sample

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
IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATHS R2_PATHS <<< "$SAMPLE_INFO"
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Lane Merge for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "lane_merge"; then
    log "Lane merge already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main merge function
merge_sample_lanes() {
    local sample_name="$1"
    local treatment="$2"
    local r1_paths="$3"  # Comma-separated
    local r2_paths="$4"  # Comma-separated
    
    log "Processing lane merge for $sample_name ($treatment)"
    
    local merge_dir="${OUTPUT_DIR}/merged_lanes/${treatment}/${sample_name}"
    mkdir -p "$merge_dir"
    
    # Output files
    local merged_r1="${merge_dir}/merged_R1.fq.gz"
    local merged_r2="${merge_dir}/merged_R2.fq.gz"
    
    # Check if already merged
    if [ -f "$merged_r1" ] && [ -f "$merged_r2" ]; then
        log "Sample $sample_name already merged, skipping..."
        return 0
    fi
    
    # Split paths into arrays
    IFS=',' read -ra R1_ARRAY <<< "$r1_paths"
    IFS=',' read -ra R2_ARRAY <<< "$r2_paths"
    
    local num_lanes=${#R1_ARRAY[@]}
    
    if [ $num_lanes -ne ${#R2_ARRAY[@]} ]; then
        log "ERROR: Mismatch in number of R1 ($num_lanes) and R2 (${#R2_ARRAY[@]}) files"
        return 1
    fi
    
    log "Found $num_lanes lane(s) for $sample_name:"
    for i in "${!R1_ARRAY[@]}"; do
        log "  Lane $((i+1)):"
        log "    R1: ${R1_ARRAY[$i]}"
        log "    R2: ${R2_ARRAY[$i]}"
    done
    
    # Validate all input files exist
    for file in "${R1_ARRAY[@]}" "${R2_ARRAY[@]}"; do
        if [ ! -f "$file" ]; then
            log "ERROR: Input file not found: $file"
            return 1
        fi
    done
    
    # Count reads in input files
    log "Counting reads in input lanes..."
    local total_r1=0
    local total_r2=0
    
    for i in "${!R1_ARRAY[@]}"; do
        local r1_count=$(count_reads "${R1_ARRAY[$i]}")
        local r2_count=$(count_reads "${R2_ARRAY[$i]}")
        
        log "  Lane $((i+1)): R1=$r1_count, R2=$r2_count reads"
        
        # Validate each lane is synchronized
        if [ "$r1_count" -ne "$r2_count" ]; then
            log "ERROR: Lane $((i+1)) has unequal read counts!"
            return 1
        fi
        
        total_r1=$((total_r1 + r1_count))
        total_r2=$((total_r2 + r2_count))
    done
    
    log "Total input reads: R1=$total_r1, R2=$total_r2"
    
    # Single lane - create symlinks
    if [ $num_lanes -eq 1 ]; then
        log "Single lane detected, creating symlinks..."
        ln -sf "$(realpath "${R1_ARRAY[0]}")" "$merged_r1"
        ln -sf "$(realpath "${R2_ARRAY[0]}")" "$merged_r2"
        log "✅ Symlinked single lane for $sample_name"
        return 0
    fi
    
    # Multiple lanes - merge using cat (safe for gzipped FASTQ)
    log "Merging $num_lanes lanes using cat..."
    
    # Merge R1 files
    log "  Concatenating R1 files..."
    cat "${R1_ARRAY[@]}" > "$merged_r1"
    
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to concatenate R1 files"
        rm -f "$merged_r1"
        return 1
    fi
    
    # Merge R2 files
    log "  Concatenating R2 files..."
    cat "${R2_ARRAY[@]}" > "$merged_r2"
    
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to concatenate R2 files"
        rm -f "$merged_r1" "$merged_r2"
        return 1
    fi
    
    # Validate merged output
    log "Validating merged output..."
    local merged_r1_count=$(count_reads "$merged_r1")
    local merged_r2_count=$(count_reads "$merged_r2")
    
    log "Merged read counts: R1=$merged_r1_count, R2=$merged_r2_count"
    
    # Check for read loss
    if [ "$merged_r1_count" -ne "$total_r1" ]; then
        log "⚠️  WARNING: R1 read count mismatch!"
        log "  Expected: $total_r1, Got: $merged_r1_count"
        log "  Lost reads: $((total_r1 - merged_r1_count))"
    fi
    
    if [ "$merged_r2_count" -ne "$total_r2" ]; then
        log "⚠️  WARNING: R2 read count mismatch!"
        log "  Expected: $total_r2, Got: $merged_r2_count"
        log "  Lost reads: $((total_r2 - merged_r2_count))"
    fi
    
    # Validate synchronization
    if [ "$merged_r1_count" -ne "$merged_r2_count" ]; then
        log "ERROR: Merged files are not synchronized!"
        log "  R1 reads: $merged_r1_count"
        log "  R2 reads: $merged_r2_count"
        return 1
    fi
    
    log "✅ Lane merge successful for $sample_name"
    log "  Merged $num_lanes lanes into $merged_r1_count paired reads"
    
    return 0
}

# Run the merge
if merge_sample_lanes "$SAMPLE_NAME" "$TREATMENT" "$R1_PATHS" "$R2_PATHS"; then
    create_sample_checkpoint "$SAMPLE_NAME" "lane_merge"
    log "====== Lane merge completed for $SAMPLE_NAME ======"
else
    log "ERROR: Lane merge failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"