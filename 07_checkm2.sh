#!/bin/bash
#SBATCH --job-name=checkm2_analysis
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00

# 07_checkm2.sh - CheckM2 quality assessment with path length fix

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

# Create a shorter temporary directory to avoid AF_UNIX path length issues
SHORT_TEMP_DIR="/tmp/checkm2_${SAMPLE_NAME}_$$"
mkdir -p "$SHORT_TEMP_DIR"

log "====== Starting CheckM2 Quality Assessment for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "checkm2"; then
    log "CheckM2 analysis already completed for $SAMPLE_NAME"
    rm -rf "$SHORT_TEMP_DIR"
    exit 0
fi

# Function to prepare bins for CheckM2
prepare_bins_for_checkm2() {
    local sample_name="$1"
    local treatment="$2"
    local magpurify_dir="$3"
    local temp_bins_dir="$4"
    
    log "Preparing bins for CheckM2 analysis..."
    
    mkdir -p "$temp_bins_dir"
    
    local bins_prepared=0
    local bins_skipped=0
    
    # Copy MAGpurify cleaned bins - look in purified_bins subdirectory
    local purified_bins_dir="${magpurify_dir}/purified_bins"
    
    if [ -d "$purified_bins_dir" ]; then
        log "  Checking purified bins directory: $purified_bins_dir"
        
        for bin_file in "$purified_bins_dir"/*.fa; do
            if [ -f "$bin_file" ]; then
                local bin_name=$(basename "$bin_file")
                local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
                local contig_count=$(grep -c "^>" "$bin_file")
                
                log "  Found bin: $bin_name ($bin_size bp, $contig_count contigs)"
                
                # Only process bins with reasonable size and contig count
                if [ $bin_size -ge 500000 ] && [ $contig_count -ge 10 ]; then
                    cp "$bin_file" "$temp_bins_dir/"
                    log "  Prepared: $bin_name ($bin_size bp, $contig_count contigs)"
                    ((bins_prepared++))
                else
                    log "  Skipped: $bin_name (too small: $bin_size bp, $contig_count contigs)"
                    ((bins_skipped++))
                fi
            fi
        done
    elif [ -d "$magpurify_dir" ]; then
        # Fallback: look for files directly in main directory
        log "  purified_bins subdirectory not found, checking main directory: $magpurify_dir"
        
        for bin_file in "$magpurify_dir"/*.fa; do
            if [ -f "$bin_file" ]; then
                local bin_name=$(basename "$bin_file")
                local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
                local contig_count=$(grep -c "^>" "$bin_file")
                
                log "  Found bin: $bin_name ($bin_size bp, $contig_count contigs)"
                
                # Only process bins with reasonable size and contig count
                if [ $bin_size -ge 500000 ] && [ $contig_count -ge 10 ]; then
                    cp "$bin_file" "$temp_bins_dir/"
                    log "  Prepared: $bin_name ($bin_size bp, $contig_count contigs)"
                    ((bins_prepared++))
                else
                    log "  Skipped: $bin_name (too small: $bin_size bp, $contig_count contigs)"
                    ((bins_skipped++))
                fi
            fi
        done
    else
        log "  ERROR: MAGpurify directory does not exist: $magpurify_dir"
    fi
    
    log "Bin preparation complete: $bins_prepared valid, $bins_skipped skipped"
    
    if [ $bins_prepared -eq 0 ]; then
        return 1
    fi
    
    return 0
}

# Function to run CheckM2 analysis
run_checkm2_analysis() {
    local sample_name="$1"
    local treatment="$2"
    local temp_bins_dir="$3"
    local output_dir="$4"
    
    log "Running CheckM2 analysis for $sample_name..."
    
    # Activate CheckM2 environment
    activate_env checkm
    
    # Check if CheckM2 is available
    if ! command -v checkm2 &> /dev/null; then
        log "ERROR: CheckM2 not available"
        conda deactivate
        return 1
    fi
    
    local bin_count=$(find "$temp_bins_dir" -name "*.fa" | wc -l)
    log "Running CheckM2 on $bin_count bins..."
    
    # Set environment variables to use shorter paths
    export TMPDIR="$SHORT_TEMP_DIR"
    export TMP="$SHORT_TEMP_DIR"
    export TEMP="$SHORT_TEMP_DIR"
    
    # Create CheckM2 working directory with short path
    local checkm2_work_dir="${SHORT_TEMP_DIR}/checkm2_work"
    mkdir -p "$checkm2_work_dir"
    
    log "Executing CheckM2 predict..."
    log "  Input directory: $temp_bins_dir"
    log "  Output directory: $output_dir"
    log "  Working directory: $checkm2_work_dir"
    log "  Threads: $SLURM_CPUS_PER_TASK"
    
    # Run CheckM2 with explicit working directory
    checkm2 predict \
        --threads $SLURM_CPUS_PER_TASK \
        --input "$temp_bins_dir" \
        --output-directory "$output_dir" \
        --tmpdir "$checkm2_work_dir" \
		-x fa \
        --force \
        2>&1 | tee "${LOG_DIR}/${treatment}/${sample_name}_checkm2.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    # Clean up working directory
    rm -rf "$checkm2_work_dir"
    
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/quality_report.tsv" ]; then
        log "CheckM2 analysis completed successfully"
        return 0
    else
        log "ERROR: CheckM2 analysis failed (exit code: $exit_code)"
        return 1
    fi
}

# Function to create CheckM2 summary
create_checkm2_summary() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="$3"
    local summary_file="${output_dir}/checkm2_summary.txt"
    
    log "Creating CheckM2 summary for $sample_name..."
    
    if [ ! -f "${output_dir}/quality_report.tsv" ]; then
        log "WARNING: CheckM2 quality report not found"
        return 1
    fi
    
    cat > "$summary_file" << EOF
CheckM2 Quality Assessment Summary for $sample_name
=================================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment

Bin Quality Assessment:
EOF
    
    # Parse CheckM2 results
    local high_quality=0
    local medium_quality=0
    local low_quality=0
    local total_bins=0
    
    # Skip header and process each bin
    tail -n +2 "${output_dir}/quality_report.tsv" | while IFS=$'\t' read -r name completeness contamination genome_size gc_content n_contigs n50 score translation_table_used; do
        ((total_bins++))
        
        echo "  ${name}:" >> "$summary_file"
        echo "    Completeness: ${completeness}%" >> "$summary_file"
        echo "    Contamination: ${contamination}%" >> "$summary_file"
        echo "    Genome size: ${genome_size} bp" >> "$summary_file"
        echo "    GC content: ${gc_content}%" >> "$summary_file"
        echo "    Contigs: ${n_contigs}" >> "$summary_file"
        echo "    N50: ${n50} bp" >> "$summary_file"
        echo "    Quality score: ${score}" >> "$summary_file"
        echo "" >> "$summary_file"
        
        # Classify quality (MIMAG standards)
        local comp_float=$(echo "$completeness" | cut -d. -f1)
        local cont_float=$(echo "$contamination" | cut -d. -f1)
        
        if [ "$comp_float" -ge 90 ] && [ "$cont_float" -le 5 ]; then
            ((high_quality++))
        elif [ "$comp_float" -ge 50 ] && [ "$cont_float" -le 10 ]; then
            ((medium_quality++))
        else
            ((low_quality++))
        fi
    done
    
    # Add summary statistics
    echo "Quality Summary:" >> "$summary_file"
    echo "  Total bins analyzed: $total_bins" >> "$summary_file"
    echo "  High quality (≥90% comp, ≤5% cont): $high_quality" >> "$summary_file"
    echo "  Medium quality (≥50% comp, ≤10% cont): $medium_quality" >> "$summary_file"
    echo "  Low quality: $low_quality" >> "$summary_file"
    
    log "CheckM2 summary created: $summary_file"
    return 0
}

# Function to find MAGpurify directory (handles both treatment-level and sample-level)
get_magpurify_dir() {
    local sample_name="$1"
    local treatment="$2"

    log "Locating MAGpurify directory for $sample_name ($treatment)..."

    # Check treatment-level directory first (for coassembly/treatment-level binning)
    local treatment_dir="${OUTPUT_DIR}/magpurify/${treatment}"
    if [ -d "$treatment_dir" ] && [ -d "${treatment_dir}/purified_bins" ]; then
        log "  Found treatment-level MAGpurify directory at: $treatment_dir"
        echo "$treatment_dir"
        return 0
    fi

    # Check sample-level directory (for individual sample binning)
    local sample_dir="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"
    if [ -d "$sample_dir" ]; then
        log "  Found sample-level MAGpurify directory at: $sample_dir"
        echo "$sample_dir"
        return 0
    fi

    log "  ERROR: No MAGpurify directory found at either location:"
    log "    Treatment-level: $treatment_dir"
    log "    Sample-level: $sample_dir"
    return 1
}

# Main processing function
stage_checkm2_analysis() {
    local sample_name="$1"
    local treatment="$2"

    log "Running CheckM2 quality assessment for $sample_name ($treatment)"

    local output_dir="${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}"
    local temp_bins_dir="${SHORT_TEMP_DIR}/bins"

    mkdir -p "$output_dir"

    # Check if already processed
    if [ -f "${output_dir}/quality_report.tsv" ]; then
        log "Sample $sample_name already analyzed, skipping..."
        return 0
    fi

    # Find MAGpurify directory - handles both treatment and sample level
    local magpurify_dir=$(get_magpurify_dir "$sample_name" "$treatment")
    if [ $? -ne 0 ] || [ -z "$magpurify_dir" ]; then
        log "ERROR: MAGpurify directory not found for $sample_name"
        return 1
    fi

    # Prepare bins for CheckM2
    if ! prepare_bins_for_checkm2 "$sample_name" "$treatment" "$magpurify_dir" "$temp_bins_dir"; then
        log "ERROR: No valid bins found for CheckM2 analysis"
        return 1
    fi
    
    # Run CheckM2 analysis
    if run_checkm2_analysis "$sample_name" "$treatment" "$temp_bins_dir" "$output_dir"; then
        log "CheckM2 analysis completed successfully"
        
        # Create summary
        create_checkm2_summary "$sample_name" "$treatment" "$output_dir"
        
        return 0
    else
        log "ERROR: CheckM2 analysis failed for $sample_name"
        return 1
    fi
}

# Validation function
validate_checkm2_analysis() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}"
    
    # Check if CheckM2 output exists and is valid
    if [ -f "${output_dir}/quality_report.tsv" ] && [ -s "${output_dir}/quality_report.tsv" ]; then
        local bin_count=$(tail -n +2 "${output_dir}/quality_report.tsv" | wc -l)
        if [ $bin_count -gt 0 ]; then
            log "Validation successful: Analyzed $bin_count bins"
            return 0
        fi
    fi
    
    log "Validation failed: No valid CheckM2 results found"
    return 1
}

# Run the CheckM2 analysis stage
if stage_checkm2_analysis "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_checkm2_analysis "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "checkm2"
        log "====== CheckM2 analysis completed successfully for $SAMPLE_NAME ======"
    else
        log "ERROR: CheckM2 validation failed for $SAMPLE_NAME"
        rm -rf "$SHORT_TEMP_DIR"
        exit 1
    fi
else
    log "ERROR: CheckM2 stage failed for $SAMPLE_NAME"
    rm -rf "$SHORT_TEMP_DIR"
    exit 1
fi

# Cleanup
rm -rf "$SHORT_TEMP_DIR"