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

# Determine processing mode: treatment-level or sample-level
if [ "${TREATMENT_LEVEL_BINNING:-false}" = "true" ]; then
    # Treatment-level mode: array index maps to treatment
    ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

    # Check if treatments file exists
    if [ ! -f "$TREATMENTS_FILE" ]; then
        log "ERROR: Treatments file not found: $TREATMENTS_FILE"
        exit 1
    fi

    # Get treatment from treatments file
    TREATMENT=$(sed -n "$((ARRAY_INDEX + 1))p" "$TREATMENTS_FILE" 2>/dev/null | tr -d '\r\n' | xargs)

    if [ -z "$TREATMENT" ]; then
        log "No treatment found for array index $ARRAY_INDEX"
        exit 0
    fi

    export TREATMENT
    PROCESSING_MODE="treatment-level"

    # Initialize for treatment-level
    init_conda
    create_treatment_dirs "$TREATMENT"
    TEMP_DIR=$(setup_temp_dir)

    log "====== Starting MAGpurify for treatment $TREATMENT (treatment-level mode) ======"

    # Check if stage already completed for treatment
    if check_treatment_checkpoint "$TREATMENT" "magpurify"; then
        log "MAGpurify already completed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

else
    # Sample-level mode: array index maps to sample (original behavior)
    SAMPLE_INFO=$(get_sample_info_by_index $SLURM_ARRAY_TASK_ID)
    if [ -z "$SAMPLE_INFO" ]; then
        log "No sample found for array index $SLURM_ARRAY_TASK_ID"
        exit 0
    fi

    # Parse sample information
    IFS='|' read -r SAMPLE_NAME TREATMENT _ _ <<< "$SAMPLE_INFO"
    export SAMPLE_NAME TREATMENT
    PROCESSING_MODE="sample-level"

    # Initialize for sample-level
    init_conda
    create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
    TEMP_DIR=$(setup_temp_dir)

    log "====== Starting MAGpurify for $SAMPLE_NAME ($TREATMENT) ======"

    # Check if stage already completed for sample
    if check_sample_checkpoint "$SAMPLE_NAME" "magpurify"; then
        log "MAGpurify already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
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

# Function to check if a bin has already been successfully purified
is_bin_already_purified() {
    local output_bin="$1"

    # Check if output bin exists
    if [ ! -f "$output_bin" ]; then
        return 1
    fi

    # Check if output bin is not empty
    if [ ! -s "$output_bin" ]; then
        return 1
    fi

    # Check if output bin has valid contigs
    local contig_count=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
    if [ $contig_count -eq 0 ]; then
        return 1
    fi

    # Bin is valid and already purified
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
    local treatment="$1"
    local sample_name="${2:-}"  # Optional for treatment-level mode

    if [ -n "$sample_name" ]; then
        # Sample-level mode
        log "Running MAGpurify for $sample_name ($treatment)"
        local output_dir="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"
        local entity_desc="sample $sample_name"
    else
        # Treatment-level mode
        log "Running MAGpurify for treatment $treatment (all bins)"
        local output_dir="${OUTPUT_DIR}/magpurify/${treatment}"
        local entity_desc="treatment $treatment"
    fi

    mkdir -p "${output_dir}/purified_bins"

    # Note: Removed completion flag check to allow resuming partial processing
    # Individual bins will be checked instead

    # Find reassembled bins - handles both treatment and sample level
    if [ -n "$sample_name" ]; then
        local input_bins_dir=$(get_reassembled_bins_dir "$sample_name" "$treatment")
    else
        # Treatment-level: look for treatment-level bins
        local input_bins_dir="${OUTPUT_DIR}/reassembly/${treatment}/reassembled_bins"
    fi

    if [ ! -d "$input_bins_dir" ] || [ -z "$(ls -A "$input_bins_dir"/*.fa 2>/dev/null)" ]; then
        log "ERROR: No reassembled bins found for $entity_desc at: $input_bins_dir"
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
    local mag_temp="${TEMP_DIR}/mag_work"
    mkdir -p "$mag_temp"

    # Get list of bins
    local bins=($(ls "$input_bins_dir"/*.fa))
    local bin_count=${#bins[@]}
    log "Processing $bin_count bins with MAGpurify for $entity_desc"

    if [ $bin_count -eq 0 ]; then
        log "ERROR: No bins found for $entity_desc"
        conda deactivate
        return 1
    fi
    
    # Process each bin
    local processed_count=0
    local cleaned_count=0
    local failed_count=0
    local skipped_count=0
    local total_contigs_before=0
    local total_contigs_after=0

    for bin_file in "${bins[@]}"; do
        local bin_name=$(basename "$bin_file" .fa)
        local output_bin="${output_dir}/purified_bins/${bin_name}.fa"

        # Check if bin has already been successfully purified
        if is_bin_already_purified "$output_bin"; then
            log "Bin $bin_name already purified, skipping..."
            ((skipped_count++))

            # Count contigs from existing output
            local existing_contigs=$(grep -c "^>" "$output_bin" 2>/dev/null || echo 0)
            total_contigs_after=$((total_contigs_after + existing_contigs))

            # Estimate input contigs for statistics
            local before_contigs=$(grep -c "^>" "$bin_file" 2>/dev/null || echo 0)
            total_contigs_before=$((total_contigs_before + before_contigs))

            continue
        fi

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

    log "MAGpurify completed for $entity_desc:"
    log "  Processed: $processed_count bins"
    log "  Skipped (already purified): $skipped_count bins"
    log "  Cleaned: $cleaned_count bins"
    log "  Failed: $failed_count bins"
    log "  Total bins: $((processed_count + skipped_count))/$bin_count"
    log "  Total contigs: $total_contigs_before → $total_contigs_after"

    # Create completion flag only if all bins have been processed or skipped
    local expected_bins=$bin_count
    local actual_bins=$((processed_count + skipped_count + failed_count))

    if [ $actual_bins -eq $expected_bins ]; then
        touch "${output_dir}/magpurify_complete.flag"
        log "All bins completed - created completion flag"
    else
        log "WARNING: Not all bins processed ($actual_bins/$expected_bins) - completion flag not created"
    fi

    # Generate statistics
    if [ -n "$sample_name" ]; then
        generate_sample_stats "$sample_name" "$treatment" "$input_bins_dir" "${output_dir}/purified_bins"
    else
        generate_treatment_stats "$treatment" "$input_bins_dir" "${output_dir}/purified_bins"
    fi

    # Return success if we have any successful bins (processed or skipped)
    if [ $((processed_count + skipped_count)) -gt 0 ]; then
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

# Generate statistics for a treatment (treatment-level mode)
generate_treatment_stats() {
    local treatment="$1"
    local input_dir="$2"
    local output_dir="$3"
    local stats_file="${OUTPUT_DIR}/magpurify/${treatment}/magpurify_stats.txt"

    cat > "$stats_file" << EOF
MAGpurify Statistics for Treatment $treatment
===========================================

Date: $(date)
Treatment: $treatment
Mode: Treatment-level

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
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        local summary_file="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/magpurify_summary.txt"
        local entity_title="$sample_name"
    else
        local summary_file="${OUTPUT_DIR}/magpurify/${treatment}/magpurify_summary.txt"
        local entity_title="Treatment $treatment"
    fi
    
    if [ -n "$sample_name" ]; then
        local results_path="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/"
    else
        local results_path="${OUTPUT_DIR}/magpurify/${treatment}/purified_bins/"
    fi

    cat > "$summary_file" << EOF
MAGpurify Summary for $entity_title
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

Results are saved in: $results_path

EOF
    
    log "MAGpurify summary created: $summary_file"
}

# Validation function
validate_magpurify() {
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        local output_dir="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"
    else
        local output_dir="${OUTPUT_DIR}/magpurify/${treatment}"
    fi
    
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

# Run the MAGpurify stage based on processing mode
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    # Treatment-level mode: run once for all bins in treatment
    if stage_magpurify "$TREATMENT"; then
        # Generate summary
        generate_magpurify_summary "$TREATMENT"

        # Validate results
        if validate_magpurify "$TREATMENT"; then
            create_treatment_checkpoint "$TREATMENT" "magpurify"
            log "====== MAGpurify completed for treatment $TREATMENT ======"
        else
            log "ERROR: MAGpurify validation failed for treatment $TREATMENT"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: MAGpurify stage failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

else
    # Sample-level mode: run for individual sample
    if stage_magpurify "$TREATMENT" "$SAMPLE_NAME"; then
        # Generate summary
        generate_magpurify_summary "$TREATMENT" "$SAMPLE_NAME"

        # Validate results
        if validate_magpurify "$TREATMENT" "$SAMPLE_NAME"; then
            create_sample_checkpoint "$SAMPLE_NAME" "magpurify"
            log "====== MAGpurify completed for $SAMPLE_NAME ======"
        else
            log "ERROR: MAGpurify validation failed for $SAMPLE_NAME"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: MAGpurify stage failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"