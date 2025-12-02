#!/bin/bash
#SBATCH --job-name=bin_refinement_dastool
#SBATCH --account=$SLURM_ACCOUNT
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=6:00:00

# 04_bin_refinement_dastool.sh - Bin refinement using DAS Tool
# Supports both sample-level and treatment-level binning

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Determine if this should run at treatment or sample level
if [ "${TREATMENT_LEVEL_BINNING:-false}" = "true" ]; then
    # Treatment-level mode - array index maps to treatment
    ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}
    TREATMENT=$(sed -n "$((ARRAY_INDEX + 1))p" "$TREATMENTS_FILE" 2>/dev/null | tr -d '\r\n' | xargs)
    
    if [ -z "$TREATMENT" ]; then
        echo "ERROR: No treatment found at index $ARRAY_INDEX" >&2
        exit 1
    fi
    
    export TREATMENT
    PROCESSING_MODE="treatment-level"
    log "Processing treatment-level refinement for: $TREATMENT"
else
    # Sample-level mode - array index maps to sample
    SAMPLE_INFO=$(get_sample_info_by_index $SLURM_ARRAY_TASK_ID)
    if [ -z "$SAMPLE_INFO" ]; then
        log "No sample found for array index $SLURM_ARRAY_TASK_ID"
        exit 0
    fi
    
    IFS='|' read -r SAMPLE_NAME TREATMENT _ _ <<< "$SAMPLE_INFO"
    export SAMPLE_NAME TREATMENT
    PROCESSING_MODE="sample-level"
    log "Processing sample-level refinement for: $SAMPLE_NAME ($TREATMENT)"
fi

# Initialize
init_conda
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    create_treatment_dirs "$TREATMENT"
else
    create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
fi
TEMP_DIR=$(setup_temp_dir)

if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    log "====== Starting Treatment-Level Bin Refinement (DAS Tool) for $TREATMENT ======"
else
    log "====== Starting Bin Refinement (DAS Tool) for $SAMPLE_NAME ($TREATMENT) ======"
fi

# Check if stage already completed
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    if check_treatment_checkpoint "$TREATMENT" "bin_refinement"; then
        log "Bin refinement already completed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
else
    if check_sample_checkpoint "$SAMPLE_NAME" "bin_refinement"; then
        log "Bin refinement already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
fi

# Function to convert bin FASTA files to DAS Tool format (contig-to-bin TSV)
create_dastool_input() {
    local binner_name="$1"
    local bins_dir="$2"
    local output_tsv="$3"
    
    log "  Converting $binner_name bins to DAS Tool format..."
    
    # DAS Tool expects TSV: contigID<TAB>binID
    > "$output_tsv"  # Clear file
    
    local bin_count=0
    for bin_file in "${bins_dir}"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        local bin_id=$(basename "$bin_file" .fa)
        
        # Extract contig IDs and assign to this bin
        grep "^>" "$bin_file" | sed 's/^>//' | while read -r contig_id; do
            echo -e "${contig_id}\t${bin_id}" >> "$output_tsv"
        done
        
        ((bin_count++))
    done
    
    if [ $bin_count -eq 0 ]; then
        log "    WARNING: No bins found in $bins_dir"
        return 1
    fi
    
    local contig_count=$(wc -l < "$output_tsv")
    log "    Created input for $bin_count bins ($contig_count contigs)"
    
    return 0
}

# Function to get assembly file
get_assembly_file() {
    local treatment="$1"
    local sample_name="$2"
    
    # For treatment-level: look for co-assembly or use first sample's assembly
    if [ "$PROCESSING_MODE" = "treatment-level" ]; then
        # Check for co-assembly
        local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
        
        for possible_file in \
            "${coassembly_dir}/contigs.fasta" \
            "${coassembly_dir}/scaffolds.fasta" \
            "${coassembly_dir}/final_assembly.fasta"; do
            
            if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                echo "$possible_file"
                return 0
            fi
        done
        
        # Fall back to first sample's assembly
        local samples=($(get_treatment_samples "$treatment"))
        if [ ${#samples[@]} -gt 0 ]; then
            local first_sample="${samples[0]}"
            local sample_assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${first_sample}"
            
            for possible_file in \
                "${sample_assembly_dir}/contigs.fasta" \
                "${sample_assembly_dir}/scaffolds.fasta"; do
                
                if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                    echo "$possible_file"
                    return 0
                fi
            done
        fi
    else
        # Sample-level: use sample's assembly
        local sample_assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
        
        for possible_file in \
            "${sample_assembly_dir}/contigs.fasta" \
            "${sample_assembly_dir}/scaffolds.fasta"; do
            
            if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                echo "$possible_file"
                return 0
            fi
        done
    fi
    
    return 1
}

# Function to run DAS Tool
run_dastool() {
    local binning_dir="$1"
    local assembly_file="$2"
    local output_dir="$3"
    local entity_name="$4"
    
    log "Running DAS Tool for $entity_name..."
    
    # Create temporary directory for DAS Tool inputs
    local dastool_input_dir="${TEMP_DIR}/dastool_inputs"
    mkdir -p "$dastool_input_dir"
    
    # Convert each binner's output to DAS Tool format
    local available_binners=()
    local binner_files=()
    local binner_labels=()
    
    # Check and convert MetaBAT2 bins
    if [ -d "${binning_dir}/metabat2_bins" ] && [ "$(find "${binning_dir}/metabat2_bins" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
        local metabat2_tsv="${dastool_input_dir}/metabat2.tsv"
        if create_dastool_input "MetaBAT2" "${binning_dir}/metabat2_bins" "$metabat2_tsv"; then
            available_binners+=("metabat2")
            binner_files+=("$metabat2_tsv")
            binner_labels+=("metabat2")
        fi
    fi
    
    # Check and convert MaxBin2 bins
    if [ -d "${binning_dir}/maxbin2_bins" ] && [ "$(find "${binning_dir}/maxbin2_bins" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
        local maxbin2_tsv="${dastool_input_dir}/maxbin2.tsv"
        if create_dastool_input "MaxBin2" "${binning_dir}/maxbin2_bins" "$maxbin2_tsv"; then
            available_binners+=("maxbin2")
            binner_files+=("$maxbin2_tsv")
            binner_labels+=("maxbin2")
        fi
    fi
    
    # Check and convert CONCOCT bins
    if [ -d "${binning_dir}/concoct_bins" ] && [ "$(find "${binning_dir}/concoct_bins" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
        local concoct_tsv="${dastool_input_dir}/concoct.tsv"
        if create_dastool_input "CONCOCT" "${binning_dir}/concoct_bins" "$concoct_tsv"; then
            available_binners+=("concoct")
            binner_files+=("$concoct_tsv")
            binner_labels+=("concoct")
        fi
    fi
    
    # Need at least 2 binners
    if [ "${#available_binners[@]}" -lt 2 ]; then
        log "ERROR: Need at least 2 binners for DAS Tool (found ${#available_binners[@]})"
        log "  Available: ${available_binners[*]}"
        return 1
    fi
    
    log "Running DAS Tool with ${#available_binners[@]} binners: ${available_binners[*]}"
    
    # Join binner files and labels with commas
    local binner_files_str=$(IFS=,; echo "${binner_files[*]}")
    local binner_labels_str=$(IFS=,; echo "${binner_labels[*]}")
    
    # Activate DAS Tool environment
    activate_env dastool
    
    # Check if DAS Tool is available
    if ! command -v DAS_Tool &> /dev/null; then
        log "ERROR: DAS Tool not available"
        conda deactivate
        return 1
    fi
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Run DAS Tool
    log "DAS Tool command:"
    log "  Input files: $binner_files_str"
    log "  Labels: $binner_labels_str"
    log "  Assembly: $assembly_file"
    log "  Output: ${output_dir}/dastool"
    log "  Threads: ${SLURM_CPUS_PER_TASK:-16}"
    log "  Score threshold: 0.5"
    log "  Search engine: diamond"
    
    DAS_Tool \
        -i "$binner_files_str" \
        -l "$binner_labels_str" \
        -c "$assembly_file" \
        -o "${output_dir}/dastool" \
        --threads ${SLURM_CPUS_PER_TASK:-16} \
        --write_bins \
        --score_threshold 0.5 \
        --search_engine diamond \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${entity_name}_dastool.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    # Check if DAS Tool was successful
    if [ $exit_code -eq 0 ]; then
        # Check for output bins
        local dastool_bins_dir="${output_dir}/dastool_DASTool_bins"
        
        if [ -d "$dastool_bins_dir" ] && [ "$(find "$dastool_bins_dir" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
            local refined_count=$(find "$dastool_bins_dir" -name "*.fa" | wc -l)
            log "DAS Tool successful: $refined_count refined bins"
            return 0
        else
            log "WARNING: DAS Tool completed but produced no bins"
            # This might happen if quality thresholds are too strict
            log "  Check ${output_dir}/dastool_DASTool_summary.tsv for details"
            return 1
        fi
    else
        log "ERROR: DAS Tool failed with exit code: $exit_code"
        return 1
    fi
}

# Generate refinement statistics
generate_refinement_stats() {
    local entity_name="$1"
    local binning_dir="$2"
    local output_dir="$3"
    local mode="$4"
    local stats_file="${output_dir}/dastool_refinement_stats.txt"
    
    log "Generating DAS Tool refinement statistics for $entity_name..."
    
    cat > "$stats_file" << EOF
DAS Tool Bin Refinement Statistics for $entity_name
====================================================

Date: $(date)
Mode: $mode
Entity: $entity_name
Treatment: $TREATMENT

Input Bins:
EOF
    
    local total_input_bins=0
    for binner in metabat2 maxbin2 concoct; do
        local binner_dir="${binning_dir}/${binner}_bins"
        if [ -d "$binner_dir" ]; then
            local bin_count=$(find "$binner_dir" -name "*.fa" 2>/dev/null | wc -l)
            if [ $bin_count -gt 0 ]; then
                echo "  $binner: $bin_count bins" >> "$stats_file"
                total_input_bins=$((total_input_bins + bin_count))
            fi
        fi
    done
    
    echo "  Total input bins: $total_input_bins" >> "$stats_file"
    echo "" >> "$stats_file"
    
    # DAS Tool output bins
    local dastool_bins_dir="${output_dir}/dastool_DASTool_bins"
    local refined_count=0
    
    if [ -d "$dastool_bins_dir" ]; then
        refined_count=$(find "$dastool_bins_dir" -name "*.fa" 2>/dev/null | wc -l)
        echo "DAS Tool Refined Bins: $refined_count" >> "$stats_file"
        echo "" >> "$stats_file"
        
        if [ $refined_count -gt 0 ]; then
            echo "Refined Bin Details:" >> "$stats_file"
            for bin in "$dastool_bins_dir"/*.fa; do
                if [ -f "$bin" ]; then
                    local bin_name=$(basename "$bin")
                    local bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    local contig_count=$(grep -c "^>" "$bin")
                    echo "  $bin_name: $bin_size bp, $contig_count contigs" >> "$stats_file"
                fi
            done
        fi
    fi
    
    echo "" >> "$stats_file"
    
    # Add DAS Tool summary if available
    local dastool_summary="${output_dir}/dastool_DASTool_summary.tsv"
    if [ -f "$dastool_summary" ]; then
        echo "DAS Tool Quality Summary:" >> "$stats_file"
        echo "  (see dastool_DASTool_summary.tsv for full details)" >> "$stats_file"
        echo "" >> "$stats_file"
        
        # Extract quality metrics
        local high_quality=$(awk -F'\t' 'NR>1 && $2>=90 && $3<=5 {count++} END {print count+0}' "$dastool_summary")
        local medium_quality=$(awk -F'\t' 'NR>1 && $2>=50 && $2<90 && $3<=10 {count++} END {print count+0}' "$dastool_summary")
        
        echo "  High quality bins (≥90% complete, ≤5% contamination): $high_quality" >> "$stats_file"
        echo "  Medium quality bins (≥50% complete, ≤10% contamination): $medium_quality" >> "$stats_file"
    fi
    
    echo "" >> "$stats_file"
    
    # Calculate statistics
    if [ $total_input_bins -gt 0 ] && command -v bc &> /dev/null; then
        local reduction=$(echo "scale=1; (($total_input_bins - $refined_count) * 100) / $total_input_bins" | bc -l)
        
        echo "Bin Dereplication:" >> "$stats_file"
        echo "  Reduction: ${reduction}% ($total_input_bins → $refined_count bins)" >> "$stats_file"
        echo "  DAS Tool eliminated $(($total_input_bins - $refined_count)) redundant/low-quality bins" >> "$stats_file"
    fi
    
    log "DAS Tool refinement statistics created: $stats_file"
}

# Validation function
validate_bin_refinement() {
    local output_dir="$1"
    
    # Check if DAS Tool bins directory exists and has content
    local dastool_bins_dir="${output_dir}/dastool_DASTool_bins"
    
    if [ -d "$dastool_bins_dir" ]; then
        local refined_count=$(find "$dastool_bins_dir" -name "*.fa" 2>/dev/null | wc -l)
        if [ $refined_count -gt 0 ]; then
            # Check that at least one bin is not empty
            for bin in "$dastool_bins_dir"/*.fa; do
                if [ -f "$bin" ] && [ -s "$bin" ]; then
                    log "Validation successful: Found $refined_count refined bins"
                    return 0
                fi
            done
        fi
    fi
    
    log "Validation failed: No valid refined bins found"
    return 1
}

# Main processing
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    # Treatment-level bin refinement
    log "Running treatment-level bin refinement (DAS Tool) for $TREATMENT"
    
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    output_dir="${OUTPUT_DIR}/bin_refinement/${TREATMENT}"
    
    # Check if already processed
    if [ -d "${output_dir}/dastool_DASTool_bins" ] && [ "$(find "${output_dir}/dastool_DASTool_bins" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
        log "Treatment $TREATMENT already has refined bins, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
    
    # Get assembly file
    assembly_file=$(get_assembly_file "$TREATMENT" "")
    if [ -z "$assembly_file" ]; then
        log "ERROR: No assembly file found for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    log "Using assembly: $assembly_file"
    
    # Validate input bins exist
    if [ ! -d "$binning_dir" ]; then
        log "ERROR: Binning directory not found for treatment $TREATMENT: $binning_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    
    # Run DAS Tool
    if run_dastool "$binning_dir" "$assembly_file" "$output_dir" "$TREATMENT"; then
        log "DAS Tool refinement completed successfully"
        generate_refinement_stats "$TREATMENT" "$binning_dir" "$output_dir" "treatment-level"
        
        if validate_bin_refinement "$output_dir"; then
            create_treatment_checkpoint "$TREATMENT" "bin_refinement"
            log "====== Treatment-level bin refinement (DAS Tool) completed for $TREATMENT ======"
        else
            log "ERROR: Bin refinement validation failed for treatment $TREATMENT"
            exit 1
        fi
    else
        log "ERROR: DAS Tool bin refinement failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    
else
    # Sample-level bin refinement
    log "Running sample-level bin refinement (DAS Tool) for $SAMPLE_NAME ($TREATMENT)"
    
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    output_dir="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/${SAMPLE_NAME}"
    
    # Check if already processed
    if [ -d "${output_dir}/dastool_DASTool_bins" ] && [ "$(find "${output_dir}/dastool_DASTool_bins" -name "*.fa" 2>/dev/null | wc -l)" -gt 0 ]; then
        log "Sample $SAMPLE_NAME already has refined bins, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
    
    # Get assembly file
    assembly_file=$(get_assembly_file "$TREATMENT" "$SAMPLE_NAME")
    if [ -z "$assembly_file" ]; then
        log "ERROR: No assembly file found for sample $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    log "Using assembly: $assembly_file"
    
    # Validate input bins exist
    if [ ! -d "$binning_dir" ]; then
        log "ERROR: Binning directory not found for $SAMPLE_NAME: $binning_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    
    # Run DAS Tool
    if run_dastool "$binning_dir" "$assembly_file" "$output_dir" "$SAMPLE_NAME"; then
        log "DAS Tool refinement completed successfully"
        generate_refinement_stats "$SAMPLE_NAME" "$binning_dir" "$output_dir" "sample-level"
        
        if validate_bin_refinement "$output_dir"; then
            create_sample_checkpoint "$SAMPLE_NAME" "bin_refinement"
            log "====== Bin refinement (DAS Tool) completed successfully for $SAMPLE_NAME ======"
        else
            log "ERROR: Bin refinement validation failed for $SAMPLE_NAME"
            exit 1
        fi
    else
        log "ERROR: DAS Tool bin refinement failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
