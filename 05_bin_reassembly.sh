#!/bin/bash
#SBATCH --job-name=bin_reassembly
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5-0

# 05_bin_reassembly.sh - Bin reassembly using MetaWRAP and SPAdes

# Source configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"

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

log "====== Starting Bin Reassembly for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "bin_reassembly"; then
    log "Bin reassembly already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Function to filter reads for a specific bin
filter_reads_for_bin() {
    local bin_file="$1"
    local read1="$2"
    local read2="$3"
    local singletons="$4"
    local output_dir="$5"
    local bin_name="$6"
    
    log "Filtering reads for bin: $bin_name"
    
    # Create bin-specific read files
    local bin_read1="${output_dir}/${bin_name}_R1.fastq.gz"
    local bin_read2="${output_dir}/${bin_name}_R2.fastq.gz"
    local bin_singletons="${output_dir}/${bin_name}_singletons.fastq.gz"
    
    # Activate environment with mapping tools
    activate_env metawrap-env
    
    # Build index for the bin
    local bin_index="${output_dir}/${bin_name}_index"
    bowtie2-build "$bin_file" "$bin_index" 2>/dev/null
    
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to build index for $bin_name"
        conda deactivate
        return 1
    fi
    
    # Map reads to bin
    local mapped_sam="${output_dir}/${bin_name}_mapped.sam"
    bowtie2 -x "$bin_index" \
        -1 "$read1" -2 "$read2" \
        -S "$mapped_sam" \
        -p 4 \
        --no-unal \
        --very-sensitive \
        2>/dev/null
    
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to map reads for $bin_name"
        conda deactivate
        return 1
    fi
    
    # Extract mapped reads
    samtools view -bS "$mapped_sam" | samtools sort -n -@ 4 | \
    samtools fastq -1 "$bin_read1" -2 "$bin_read2" -s "$bin_singletons" -
    
    # Clean up
    rm -f "$mapped_sam" "${bin_index}".* 2>/dev/null
    
    conda deactivate
    
    # Check if we got any reads
    local read_count=0
    if [ -f "$bin_read1" ] && [ -s "$bin_read1" ]; then
        read_count=$(zcat "$bin_read1" | wc -l | awk '{print $1/4}')
    fi
    
    log "  Extracted $read_count read pairs for $bin_name"
    
    if [ $read_count -lt 100 ]; then
        log "  WARNING: Very few reads extracted for $bin_name"
        return 1
    fi
    
    return 0
}

# Function to reassemble a single bin using SPAdes
reassemble_bin_spades() {
    local bin_file="$1"
    local bin_read1="$2"
    local bin_read2="$3"
    local bin_singletons="$4"
    local output_dir="$5"
    local bin_name="$6"
    
    log "Reassembling $bin_name with SPAdes..."
    
    # Activate SPAdes environment
    activate_env spades
    
    # Create output directory for this bin
    local bin_output="${output_dir}/${bin_name}_spades"
    mkdir -p "$bin_output"
    
    # Prepare SPAdes command
    local spades_cmd="spades.py --isolate"  # Use isolate mode for better bin assembly
    spades_cmd+=" -1 $bin_read1"
    spades_cmd+=" -2 $bin_read2"
    
    # Add singletons if available and non-empty
    if [ -f "$bin_singletons" ] && [ -s "$bin_singletons" ]; then
        spades_cmd+=" -s $bin_singletons"
    fi
    
    spades_cmd+=" -o $bin_output"
    spades_cmd+=" -t 8"  # Use fewer threads per bin
    spades_cmd+=" -m 32"  # Use less memory per bin
    spades_cmd+=" --tmp-dir ${TEMP_DIR}/spades_${bin_name}"
    
    # Create temp directory
    mkdir -p "${TEMP_DIR}/spades_${bin_name}"
    
    # Run SPAdes
    eval $spades_cmd 2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_${bin_name}_spades.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${bin_output}/contigs.fasta" ]; then
        # Copy the reassembled bin
        cp "${bin_output}/contigs.fasta" "${output_dir}/${bin_name}.fa"
        
        # Clean up SPAdes output to save space
        rm -rf "${bin_output}" "${TEMP_DIR}/spades_${bin_name}"
        
        log "  Successfully reassembled $bin_name"
        return 0
    else
        log "  ERROR: SPAdes failed for $bin_name"
        rm -rf "${bin_output}" "${TEMP_DIR}/spades_${bin_name}"
        return 1
    fi
}

# Function to run MetaWRAP bin reassembly
run_metawrap_reassembly() {
    local sample_name="$1"
    local treatment="$2"
    local refined_bins_dir="$3"
    local read1="$4"
    local read2="$5"
    local output_dir="$6"
    
    log "Running MetaWRAP bin reassembly for $sample_name..."
    
    # Activate MetaWRAP environment
    activate_env metawrap-env
    
    # Check if MetaWRAP is available
    if ! command -v metawrap &> /dev/null; then
        log "WARNING: MetaWRAP not available for bin reassembly"
        conda deactivate
        return 1
    fi
    
    # Create temporary directory for this sample
    local temp_reassembly="${TEMP_DIR}/reassembly_${sample_name}"
    mkdir -p "$temp_reassembly"
    
    # Run MetaWRAP reassembly
    metawrap reassemble_bins \
        -o "$output_dir" \
        -t $SLURM_CPUS_PER_TASK \
        -1 "$read1" \
        -2 "$read2" \
        -c 50 \
        -x 10 \
        -b "$refined_bins_dir" \
        2>&1 | tee "${LOG_DIR}/${treatment}/${sample_name}_metawrap_reassembly.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/reassembled_bins" ]; then
        local reassembled_count=$(ls -1 "${output_dir}/reassembled_bins"/*.fa 2>/dev/null | wc -l)
        log "MetaWRAP reassembly successful: $reassembled_count reassembled bins"
        return 0
    else
        log "WARNING: MetaWRAP reassembly failed (exit code: $exit_code)"
        return 1
    fi
}

# Function to run manual bin reassembly
run_manual_reassembly() {
    local sample_name="$1"
    local treatment="$2"
    local refined_bins_dir="$3"
    local read1="$4"
    local read2="$5"
    local singletons="$6"
    local output_dir="$7"
    
    log "Running manual bin reassembly for $sample_name..."
    
    local reassembled_dir="${output_dir}/reassembled_bins"
    mkdir -p "$reassembled_dir"
    
    local temp_reads="${TEMP_DIR}/reads_${sample_name}"
    mkdir -p "$temp_reads"
    
    local reassembled_count=0
    local failed_count=0
    local copied_count=0
    
    # Process each refined bin
    for bin_file in "$refined_bins_dir"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        local bin_name=$(basename "$bin_file" .fa)
        log "Processing bin: $bin_name"
        
        # Check bin size - only reassemble larger bins
        local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
        local contig_count=$(grep -c "^>" "$bin_file")
        
        if [ $bin_size -lt 1000000 ]; then  # < 1MB
            log "  Bin too small ($bin_size bp), copying original"
            cp "$bin_file" "${reassembled_dir}/${bin_name}.fa"
            ((copied_count++))
            continue
        fi
        
        # Filter reads for this bin
        if filter_reads_for_bin "$bin_file" "$read1" "$read2" "$singletons" "$temp_reads" "$bin_name"; then
            # Try reassembly
            if reassemble_bin_spades "$bin_file" \
                "${temp_reads}/${bin_name}_R1.fastq.gz" \
                "${temp_reads}/${bin_name}_R2.fastq.gz" \
                "${temp_reads}/${bin_name}_singletons.fastq.gz" \
                "$temp_reads" "$bin_name"; then
                
                # Move reassembled bin to final location
                mv "${temp_reads}/${bin_name}.fa" "${reassembled_dir}/${bin_name}.fa"
                ((reassembled_count++))
            else
                log "  Reassembly failed for $bin_name, copying original"
                cp "$bin_file" "${reassembled_dir}/${bin_name}.fa"
                ((failed_count++))
            fi
        else
            log "  Read filtering failed for $bin_name, copying original"
            cp "$bin_file" "${reassembled_dir}/${bin_name}.fa"
            ((failed_count++))
        fi
        
        # Clean up bin-specific files
        rm -f "${temp_reads}/${bin_name}"_* 2>/dev/null
    done
    
    log "Manual reassembly completed:"
    log "  Reassembled: $reassembled_count bins"
    log "  Failed (copied original): $failed_count bins"
    log "  Too small (copied original): $copied_count bins"
    
    local total_bins=$((reassembled_count + failed_count + copied_count))
    
    if [ $total_bins -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

# Main processing function
stage_bin_reassembly() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Running bin reassembly for $sample_name ($treatment)"
    
    local refinement_dir="${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}"
    local output_dir="${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    
    mkdir -p "$output_dir"
    
    # Check if already processed
    if [ -d "${output_dir}/reassembled_bins" ] && [ "$(ls -A "${output_dir}/reassembled_bins/"*.fa 2>/dev/null)" ]; then
        log "Sample $sample_name already reassembled, skipping..."
        return 0
    fi
    
    # Check for refined bins
    local refined_bins_dir="${refinement_dir}/metawrap_50_10_bins"
    if [ ! -d "$refined_bins_dir" ] || [ ! "$(ls -A "$refined_bins_dir"/*.fa 2>/dev/null)" ]; then
        log "ERROR: No refined bins found for $sample_name"
        log "  Expected: $refined_bins_dir"
        return 1
    fi
    
    # Check for reads
    local read1="${quality_dir}/filtered_1.fastq.gz"
    local read2="${quality_dir}/filtered_2.fastq.gz"
    local singletons="${quality_dir}/singletons.fastq.gz"
    
    if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
        log "ERROR: Missing quality-filtered reads for $sample_name"
        log "  Expected: $read1 and $read2"
        return 1
    fi
    
    # Count bins to reassemble
    local bin_count=$(ls -1 "$refined_bins_dir"/*.fa 2>/dev/null | wc -l)
    log "Reassembling $bin_count refined bins for $sample_name"
    
    if [ $bin_count -eq 0 ]; then
        log "ERROR: No bins to reassemble for $sample_name"
        return 1
    fi
    
    # Try MetaWRAP reassembly first
    if run_metawrap_reassembly "$sample_name" "$treatment" "$refined_bins_dir" "$read1" "$read2" "$output_dir"; then
        log "MetaWRAP reassembly completed successfully"
    else
        log "MetaWRAP reassembly failed, trying manual reassembly..."
        
        if run_manual_reassembly "$sample_name" "$treatment" "$refined_bins_dir" "$read1" "$read2" "$singletons" "$output_dir"; then
            log "Manual reassembly completed successfully"
        else
            log "ERROR: Both MetaWRAP and manual reassembly failed"
            return 1
        fi
    fi
    
    # Generate reassembly statistics
    generate_reassembly_stats "$sample_name" "$treatment" "$refined_bins_dir" "${output_dir}/reassembled_bins"
    
    return 0
}

# Generate reassembly statistics
generate_reassembly_stats() {
    local sample_name="$1"
    local treatment="$2"
    local original_dir="$3"
    local reassembled_dir="$4"
    local stats_file="${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}/reassembly_stats.txt"
    
    cat > "$stats_file" << EOF
Bin Reassembly Statistics for $sample_name
==========================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment

EOF
    
    # Count bins
    local original_count=$(ls -1 "$original_dir"/*.fa 2>/dev/null | wc -l)
    local reassembled_count=$(ls -1 "$reassembled_dir"/*.fa 2>/dev/null | wc -l)
    
    echo "Bin Counts:" >> "$stats_file"
    echo "  Original refined bins: $original_count" >> "$stats_file"
    echo "  Reassembled bins: $reassembled_count" >> "$stats_file"
    echo "" >> "$stats_file"
    
    # Compare bin sizes if both directories exist
    if [ -d "$original_dir" ] && [ -d "$reassembled_dir" ]; then
        echo "Bin Size Comparison:" >> "$stats_file"
        echo "Bin_Name\tOriginal_Size\tReassembled_Size\tSize_Change\tContigs_Before\tContigs_After" >> "$stats_file"
        
        local total_original_size=0
        local total_reassembled_size=0
        local improved_count=0
        local degraded_count=0
        
        for original_bin in "$original_dir"/*.fa; do
            if [ -f "$original_bin" ]; then
                local bin_name=$(basename "$original_bin" .fa)
                local reassembled_bin="${reassembled_dir}/${bin_name}.fa"
                
                if [ -f "$reassembled_bin" ]; then
                    local original_size=$(grep -v '^>' "$original_bin" | tr -d '\n' | wc -c)
                    local reassembled_size=$(grep -v '^>' "$reassembled_bin" | tr -d '\n' | wc -c)
                    local size_change=$((reassembled_size - original_size))
                    
                    local original_contigs=$(grep -c '^>' "$original_bin")
                    local reassembled_contigs=$(grep -c '^>' "$reassembled_bin")
                    
                    echo "$bin_name\t$original_size\t$reassembled_size\t$size_change\t$original_contigs\t$reassembled_contigs" >> "$stats_file"
                    
                    total_original_size=$((total_original_size + original_size))
                    total_reassembled_size=$((total_reassembled_size + reassembled_size))
                    
                    if [ $size_change -gt 0 ]; then
                        ((improved_count++))
                    elif [ $size_change -lt 0 ]; then
                        ((degraded_count++))
                    fi
                fi
            fi
        done
        
        echo "" >> "$stats_file"
        echo "Summary:" >> "$stats_file"
        echo "  Total original size: $total_original_size bp" >> "$stats_file"
        echo "  Total reassembled size: $total_reassembled_size bp" >> "$stats_file"
        echo "  Net size change: $((total_reassembled_size - total_original_size)) bp" >> "$stats_file"
        echo "  Bins improved: $improved_count" >> "$stats_file"
        echo "  Bins degraded: $degraded_count" >> "$stats_file"
        echo "  Bins unchanged: $((reassembled_count - improved_count - degraded_count))" >> "$stats_file"
        
        if [ $reassembled_count -gt 0 ]; then
            local success_rate=$(echo "scale=1; ($reassembled_count * 100) / $original_count" | bc -l)
            echo "  Success rate: ${success_rate}%" >> "$stats_file"
        fi
    fi
    
    log "Reassembly statistics created: $stats_file"
}

# Validation function
validate_bin_reassembly() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}"
    
    # Check if reassembled bins directory exists and has content
    if [ -d "${output_dir}/reassembled_bins" ] && [ "$(ls -A "${output_dir}/reassembled_bins/"*.fa 2>/dev/null)" ]; then
        # Check that reassembled bins are not empty
        for bin in "${output_dir}/reassembled_bins"/*.fa; do
            if [ -f "$bin" ] && [ -s "$bin" ]; then
                return 0  # At least one non-empty bin found
            fi
        done
    fi
    
    return 1
}

# Run the bin reassembly stage
if stage_bin_reassembly "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_bin_reassembly "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "bin_reassembly"
        log "====== Bin reassembly completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Bin reassembly validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: Bin reassembly stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"