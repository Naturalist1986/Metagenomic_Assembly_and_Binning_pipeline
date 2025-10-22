#!/bin/bash
#SBATCH --job-name=magpurify
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# 06_magpurify.sh - MAGpurify contamination removal stage

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

log "====== Starting MAGpurify for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "magpurify"; then
    log "MAGpurify already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Function to clean bins by removing contigs with ambiguous bases
clean_bin_sequences() {
    local input_bin="$1"
    local output_bin="$2"
    
    log "  Pre-cleaning sequences in $(basename "$input_bin")..."
    
    # Remove contigs with ambiguous bases and very short contigs
    awk '
        BEGIN {RS=">"; ORS=""}
        NR>1 {
            header = $0
            gsub(/\n.*/, "", header)
            seq = $0
            gsub(/^[^\n]*\n/, "", seq)
            gsub(/\n/, "", seq)
            
            # Only keep sequences with standard nucleotides and length >= 1000bp
            if (!match(seq, /[^ATGCatgc]/) && length(seq) >= 1000) {
                print ">" header "\n" seq "\n"
            }
        }
    ' "$input_bin" > "$output_bin"
    
    local input_contigs=$(grep -c "^>" "$input_bin" 2>/dev/null || echo 0)
    local output_contigs=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
    
    log "    Contigs: $input_contigs → $output_contigs (removed $((input_contigs - output_contigs)) short/ambiguous contigs)"
    
    if [ $output_contigs -eq 0 ]; then
        log "    WARNING: No contigs remaining after cleaning"
        return 1
    fi
    
    return 0
}

# Function to run MAGpurify modules on a single bin
run_magpurify_on_bin() {
    local input_bin="$1"
    local work_dir="$2"
    local bin_name="$3"
    local output_bin="$4"
    
    log "  Running MAGpurify on $bin_name..."
    
    # Create bin-specific work directory
    local bin_work_dir="${work_dir}/${bin_name}"
    mkdir -p "$bin_work_dir"
    
    # Clean input bin first
    local cleaned_bin="${TEMP_DIR}/${bin_name}_cleaned.fa"
    if ! clean_bin_sequences "$input_bin" "$cleaned_bin"; then
        log "    ERROR: Bin cleaning failed for $bin_name"
        return 1
    fi
    
    # Check if bin is large enough for MAGpurify
    local bin_size=$(grep -v '^>' "$cleaned_bin" | tr -d '\n' | wc -c)
    local contig_count=$(grep -c '^>' "$cleaned_bin")
    
    if [ $bin_size -lt 500000 ] || [ $contig_count -lt 5 ]; then
        log "    Bin too small for MAGpurify ($bin_size bp, $contig_count contigs), copying as-is"
        cp "$cleaned_bin" "$output_bin"
        return 0
    fi
    
    log "    Processing bin: $bin_size bp, $contig_count contigs"
    
    # Run MAGpurify modules
    local modules_run=0
    local modules_success=0
    local flagged_contigs_found=false
    
    # Phylo-markers module
    log "    Running phylo-markers module..."
    if magpurify phylo-markers "$cleaned_bin" "$bin_work_dir" 2>/dev/null; then
        ((modules_run++))
        if [ -f "${bin_work_dir}/phylo-markers/flagged_contigs" ] && [ -s "${bin_work_dir}/phylo-markers/flagged_contigs" ]; then
            flagged_contigs_found=true
            ((modules_success++))
            log "      Found $(wc -l < "${bin_work_dir}/phylo-markers/flagged_contigs") flagged contigs"
        fi
    else
        log "      WARNING: phylo-markers module failed"
    fi
    
    # Clade-markers module
    log "    Running clade-markers module..."
    if magpurify clade-markers "$cleaned_bin" "$bin_work_dir" 2>/dev/null; then
        ((modules_run++))
        if [ -f "${bin_work_dir}/clade-markers/flagged_contigs" ] && [ -s "${bin_work_dir}/clade-markers/flagged_contigs" ]; then
            flagged_contigs_found=true
            ((modules_success++))
            log "      Found $(wc -l < "${bin_work_dir}/clade-markers/flagged_contigs") flagged contigs"
        fi
    else
        log "      WARNING: clade-markers module failed"
    fi
    
    # Tetra-freq module
    log "    Running tetra-freq module..."
    if magpurify tetra-freq "$cleaned_bin" "$bin_work_dir" 2>/dev/null; then
        ((modules_run++))
        if [ -f "${bin_work_dir}/tetra-freq/flagged_contigs" ] && [ -s "${bin_work_dir}/tetra-freq/flagged_contigs" ]; then
            flagged_contigs_found=true
            ((modules_success++))
            log "      Found $(wc -l < "${bin_work_dir}/tetra-freq/flagged_contigs") flagged contigs"
        fi
    else
        log "      WARNING: tetra-freq module failed"
    fi
    
    # GC-content module
    log "    Running gc-content module..."
    if magpurify gc-content "$cleaned_bin" "$bin_work_dir" 2>/dev/null; then
        ((modules_run++))
        if [ -f "${bin_work_dir}/gc-content/flagged_contigs" ] && [ -s "${bin_work_dir}/gc-content/flagged_contigs" ]; then
            flagged_contigs_found=true
            ((modules_success++))
            log "      Found $(wc -l < "${bin_work_dir}/gc-content/flagged_contigs") flagged contigs"
        fi
    else
        log "      WARNING: gc-content module failed"
    fi
    
    log "    MAGpurify modules: $modules_success/$modules_run successful"
    
    # Run clean-bin if any modules flagged contigs
    if [ "$flagged_contigs_found" = true ]; then
        log "    Running clean-bin to remove flagged contigs..."
        if magpurify clean-bin "$cleaned_bin" "$bin_work_dir" "$output_bin" 2>/dev/null; then
            local cleaned_contigs=$(grep -c '^>' "$output_bin" 2>/dev/null || echo 0)
            log "    Successfully cleaned bin: $contig_count → $cleaned_contigs contigs"
            return 0
        else
            log "    WARNING: clean-bin failed, copying cleaned input"
            cp "$cleaned_bin" "$output_bin"
            return 1
        fi
    else
        log "    No contigs flagged for removal, copying cleaned input"
        cp "$cleaned_bin" "$output_bin"
        return 0
    fi
}

# Function to find reassembled bins (handles both treatment-level and sample-level)
get_reassembled_bins_dir() {
    local sample_name="$1"
    local treatment="$2"

    log "Locating reassembled bins for $sample_name ($treatment)..."

    # Check treatment-level bins first (for coassembly/treatment-level binning)
    local treatment_bins="${OUTPUT_DIR}/reassembly/${treatment}/reassembled_bins"
    if [ -d "$treatment_bins" ] && [ "$(ls -A "$treatment_bins"/*.fa 2>/dev/null)" ]; then
        local bin_count=$(ls -1 "$treatment_bins"/*.fa 2>/dev/null | wc -l)
        log "  Found $bin_count treatment-level reassembled bins at: $treatment_bins"
        echo "$treatment_bins"
        return 0
    fi

    # Check sample-level bins (for individual sample binning)
    local sample_bins="${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}/reassembled_bins"
    if [ -d "$sample_bins" ] && [ "$(ls -A "$sample_bins"/*.fa 2>/dev/null)" ]; then
        local bin_count=$(ls -1 "$sample_bins"/*.fa 2>/dev/null | wc -l)
        log "  Found $bin_count sample-level reassembled bins at: $sample_bins"
        echo "$sample_bins"
        return 0
    fi

    log "  ERROR: No reassembled bins found at either location:"
    log "    Treatment-level: $treatment_bins"
    log "    Sample-level: $sample_bins"
    return 1
}

# Main processing function
stage_magpurify() {
    local sample_name="$1"
    local treatment="$2"

    log "Running MAGpurify for $sample_name ($treatment)"

    local output_dir="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"

    mkdir -p "${output_dir}/purified_bins"

    # Check if already processed
    if [ -f "${output_dir}/magpurify_complete.flag" ]; then
        log "Sample $sample_name already processed, skipping..."
        return 0
    fi

    # Find reassembled bins - handles both treatment and sample level
    local input_bins_dir=$(get_reassembled_bins_dir "$sample_name" "$treatment")
    if [ $? -ne 0 ] || [ -z "$input_bins_dir" ]; then
        log "ERROR: No reassembled bins found for $sample_name"
        return 1
    fi
    
    # Activate MAGpurify environment
    activate_env metagenome_assembly
    
    # Check if MAGpurify is available
    if ! command -v magpurify &> /dev/null; then
        log "ERROR: MAGpurify not available in environment"
        conda deactivate
        return 1
    fi
    
    # Set up temporary workspace
    local mag_temp="${TEMP_DIR}/mag_${sample_name}"
    mkdir -p "$mag_temp"
    
    # Get list of bins
    local bins=($(ls "$input_bins_dir"/*.fa))
    local bin_count=${#bins[@]}
    log "Processing $bin_count bins with MAGpurify for $sample_name"
    
    if [ $bin_count -eq 0 ]; then
        log "ERROR: No bins found for $sample_name"
        conda deactivate
        return 1
    fi
    
    # Process each bin
    local processed_count=0
    local cleaned_count=0
    local failed_count=0
    local total_contigs_before=0
    local total_contigs_after=0
    
    for bin_file in "${bins[@]}"; do
        local bin_name=$(basename "$bin_file" .fa)
        local output_bin="${output_dir}/purified_bins/${bin_name}.fa"
        
        log "Processing bin: $bin_name"
        
        # Count input contigs
        local before_contigs=$(grep -c "^>" "$bin_file" 2>/dev/null || echo 0)
        total_contigs_before=$((total_contigs_before + before_contigs))
        
        # Run MAGpurify on this bin
        if run_magpurify_on_bin "$bin_file" "$mag_temp" "$bin_name" "$output_bin"; then
            ((processed_count++))
            
            # Count output contigs
            local after_contigs=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
            total_contigs_after=$((total_contigs_after + after_contigs))
            
            if [ $after_contigs -lt $before_contigs ]; then
                ((cleaned_count++))
                log "  SUCCESS: Cleaned $bin_name ($before_contigs → $after_contigs contigs)"
            else
                log "  SUCCESS: No contamination found in $bin_name"
            fi
        else
            ((failed_count++))
            log "  FAILED: MAGpurify failed for $bin_name"
            
            # Copy original as fallback
            if [ -f "$bin_file" ]; then
                cp "$bin_file" "$output_bin"
                local fallback_contigs=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
                total_contigs_after=$((total_contigs_after + fallback_contigs))
            fi
        fi
        
        # Clean up bin-specific temp files
        rm -rf "${mag_temp}/${bin_name}" 2>/dev/null
    done
    
    conda deactivate
    
    log "MAGpurify completed for $sample_name:"
    log "  Processed: $processed_count/$bin_count bins"
    log "  Cleaned: $cleaned_count bins"
    log "  Failed: $failed_count bins"
    log "  Total contigs: $total_contigs_before → $total_contigs_after"
    
    # Create completion flag
    touch "${output_dir}/magpurify_complete.flag"
    
    # Generate sample statistics
    generate_sample_stats "$sample_name" "$treatment" "$input_bins_dir" "${output_dir}/purified_bins"
    
    if [ $processed_count -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

# Generate statistics for a sample
generate_sample_stats() {
    local sample_name="$1"
    local treatment="$2"
    local input_dir="$3"
    local output_dir="$4"
    local stats_file="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/magpurify_stats.txt"
    
    cat > "$stats_file" << EOF
MAGpurify Statistics for $sample_name
====================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment

Bin Processing Results:
EOF
    
    local total_contigs_before=0
    local total_contigs_after=0
    local total_size_before=0
    local total_size_after=0
    local bins_cleaned=0
    local bins_processed=0
    
    echo "Bin_Name\tContigs_Before\tContigs_After\tContigs_Removed\tSize_Before\tSize_After\tSize_Change" >> "$stats_file"
    
    for input_bin in "$input_dir"/*.fa; do
        if [ -f "$input_bin" ]; then
            local bin_name=$(basename "$input_bin" .fa)
            local output_bin="${output_dir}/${bin_name}.fa"
            
            if [ -f "$output_bin" ]; then
                local before_contigs=$(grep -c "^>" "$input_bin" 2>/dev/null || echo 0)
                local after_contigs=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
                local removed_contigs=$((before_contigs - after_contigs))
                
                local before_size=$(grep -v "^>" "$input_bin" | tr -d '\n' | wc -c)
                local after_size=$(grep -v "^>" "$output_bin" | tr -d '\n' | wc -c)
                local size_change=$((after_size - before_size))
                
                echo "$bin_name\t$before_contigs\t$after_contigs\t$removed_contigs\t$before_size\t$after_size\t$size_change" >> "$stats_file"
                
                total_contigs_before=$((total_contigs_before + before_contigs))
                total_contigs_after=$((total_contigs_after + after_contigs))
                total_size_before=$((total_size_before + before_size))
                total_size_after=$((total_size_after + after_size))
                
                if [ $removed_contigs -gt 0 ]; then
                    ((bins_cleaned++))
                fi
                ((bins_processed++))
            fi
        fi
    done
    
    echo "" >> "$stats_file"
    echo "Summary:" >> "$stats_file"
    echo "  Bins processed: $bins_processed" >> "$stats_file"
    echo "  Bins cleaned: $bins_cleaned" >> "$stats_file"
    echo "  Total contigs: $total_contigs_before → $total_contigs_after ($((total_contigs_after - total_contigs_before)) change)" >> "$stats_file"
    echo "  Total size: $total_size_before → $total_size_after bp ($((total_size_after - total_size_before)) bp change)" >> "$stats_file"
    
    if [ $bins_processed -gt 0 ]; then
        local cleaning_rate=$(echo "scale=1; ($bins_cleaned * 100) / $bins_processed" | bc -l)
        echo "  Cleaning rate: ${cleaning_rate}%" >> "$stats_file"
    fi
    
    if [ $total_contigs_before -gt 0 ]; then
        local contigs_removed_pct=$(echo "scale=1; (($total_contigs_before - $total_contigs_after) * 100) / $total_contigs_before" | bc -l)
        echo "  Contigs removed: ${contigs_removed_pct}%" >> "$stats_file"
    fi
    
    log "MAGpurify statistics created: $stats_file"
}

# Generate overall statistics report
generate_magpurify_summary() {
    local sample_name="$1"
    local treatment="$2"
    local summary_file="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/magpurify_summary.txt"
    
    cat > "$summary_file" << EOF
MAGpurify Summary for $sample_name
=================================

MAGpurify is a tool for removing contaminating contigs from metagenome-assembled genomes (MAGs).

Modules Used:
  1. phylo-markers: Identifies contigs with conflicting phylogenetic signals
  2. clade-markers: Removes contigs from different taxonomic clades
  3. tetra-freq: Flags contigs with different tetranucleotide frequencies
  4. gc-content: Identifies contigs with unusual GC content

Quality Thresholds:
  - Minimum bin size: 500 KB
  - Minimum contigs: 5
  - Minimum contig length: 1000 bp

Processing Steps:
  1. Pre-cleaning: Remove ambiguous bases and short contigs
  2. Module execution: Run all MAGpurify modules
  3. Contig flagging: Identify contaminating contigs
  4. Clean-bin: Remove flagged contigs
  5. Quality check: Validate output bins

Results are saved in: ${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/

EOF
    
    log "MAGpurify summary created: $summary_file"
}

# Validation function
validate_magpurify() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"
    
    # Check completion flag
    if [ ! -f "${output_dir}/magpurify_complete.flag" ]; then
        return 1
    fi
    
    # Check purified bins directory
    if [ ! -d "${output_dir}/purified_bins" ]; then
        return 1
    fi
    
    # Check that we have at least one purified bin
    local purified_count=$(ls -1 "${output_dir}/purified_bins"/*.fa 2>/dev/null | wc -l)
    if [ $purified_count -eq 0 ]; then
        return 1
    fi
    
    # Check that purified bins are not empty
    for bin in "${output_dir}/purified_bins"/*.fa; do
        if [ -f "$bin" ] && [ -s "$bin" ]; then
            return 0  # At least one non-empty bin found
        fi
    done
    
    return 1
}

# Run the MAGpurify stage
if stage_magpurify "$SAMPLE_NAME" "$TREATMENT"; then
    # Generate summary
    generate_magpurify_summary "$SAMPLE_NAME" "$TREATMENT"
    
    # Validate results
    if validate_magpurify "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "magpurify"
        log "====== MAGpurify completed for $SAMPLE_NAME ======"
    else
        log "ERROR: MAGpurify validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: MAGpurify stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"