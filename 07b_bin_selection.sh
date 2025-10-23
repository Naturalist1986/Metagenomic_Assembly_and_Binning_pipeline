#!/bin/bash
#SBATCH --job-name=bin_selection
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

# 07b_bin_selection.sh - Select best bin version based on CheckM2 quality

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

    log "====== Starting Bin Selection for treatment $TREATMENT (treatment-level mode) ======"

    # Check if stage already completed for treatment
    if check_treatment_checkpoint "$TREATMENT" "bin_selection"; then
        log "Bin selection already completed for treatment $TREATMENT"
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

    log "====== Starting Bin Selection for $SAMPLE_NAME ($TREATMENT) ======"

    # Check if stage already completed for sample
    if check_sample_checkpoint "$SAMPLE_NAME" "bin_selection"; then
        log "Bin selection already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
fi

# Function to extract quality metrics for a bin from CheckM2 report
get_bin_quality() {
    local checkm2_report="$1"
    local bin_name="$2"

    # Look for this bin in the CheckM2 report
    local bin_line=$(grep "^${bin_name}" "$checkm2_report" | head -1)

    if [ -z "$bin_line" ]; then
        echo ""
        return 1
    fi

    # Extract completeness (column 2) and contamination (column 3)
    local completeness=$(echo "$bin_line" | awk '{print $2}')
    local contamination=$(echo "$bin_line" | awk '{print $3}')

    echo "${completeness}|${contamination}"
    return 0
}

# Function to compare two bin versions and select the better one
# Returns: 0 if first is better, 1 if second is better
compare_bin_quality() {
    local comp1="$1"
    local cont1="$2"
    local comp2="$3"
    local cont2="$4"

    # Calculate quality score: completeness - 5*contamination (MIMAG standard)
    local score1=$(echo "scale=2; $comp1 - (5 * $cont1)" | bc -l)
    local score2=$(echo "scale=2; $comp2 - (5 * $cont2)" | bc -l)

    # Compare scores
    local comparison=$(echo "$score1 > $score2" | bc -l)

    if [ "$comparison" -eq 1 ]; then
        return 0  # First is better
    else
        return 1  # Second is better
    fi
}

# Main processing function
stage_bin_selection() {
    local treatment="$1"
    local sample_name="${2:-}"  # Optional for treatment-level mode

    if [ -n "$sample_name" ]; then
        # Sample-level mode
        log "Selecting best bin versions for $sample_name ($treatment)"
        local checkm2_dir="${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}"
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}/${sample_name}"
        local entity_desc="sample $sample_name"
    else
        # Treatment-level mode
        log "Selecting best bin versions for treatment $treatment"
        local checkm2_dir="${OUTPUT_DIR}/checkm2/${treatment}"
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}"
        local entity_desc="treatment $treatment"
    fi

    mkdir -p "$output_dir"

    # Check if already processed
    if [ -f "${output_dir}/bin_selection_complete.flag" ]; then
        log "Bin selection already completed for $entity_desc, skipping..."
        return 0
    fi

    # Check for CheckM2 quality report
    local checkm2_report="${checkm2_dir}/quality_report.tsv"
    if [ ! -f "$checkm2_report" ]; then
        log "ERROR: CheckM2 quality report not found: $checkm2_report"
        return 1
    fi

    # Get list of unique base bin names (without .orig/.strict/.permissive suffixes)
    log "Analyzing CheckM2 results to identify bin versions..."

    local base_bins=()
    while IFS=$'\t' read -r name rest; do
        # Skip header
        if [ "$name" = "Name" ]; then
            continue
        fi

        # Extract base name by removing version suffixes
        local base_name="$name"
        base_name="${base_name%.orig}"
        base_name="${base_name%.strict}"
        base_name="${base_name%.permissive}"

        # Add to array if not already present
        if [[ ! " ${base_bins[@]} " =~ " ${base_name} " ]]; then
            base_bins+=("$base_name")
        fi
    done < "$checkm2_report"

    local total_bins=${#base_bins[@]}
    log "Found $total_bins unique bins with reassembly versions"

    if [ $total_bins -eq 0 ]; then
        log "ERROR: No bins found in CheckM2 report"
        return 1
    fi

    # Create selection report
    local selection_report="${output_dir}/bin_selection_report.txt"
    cat > "$selection_report" << EOF
Bin Selection Report for $entity_desc
=====================================

Date: $(date)
Treatment: $treatment
$([ -n "$sample_name" ] && echo "Sample: $sample_name")

Selection Criteria:
  Quality Score = Completeness - (5 × Contamination)
  Higher score indicates better quality

Bin Selection Results:
Base_Bin\tSelected_Version\tCompleteness\tContamination\tQuality_Score\tReason
EOF

    # Process each base bin
    local selected_count=0
    local orig_selected=0
    local strict_selected=0
    local permissive_selected=0

    for base_bin in "${base_bins[@]}"; do
        log "Processing bin: $base_bin"

        # Check which versions exist and get their quality metrics
        local orig_quality=$(get_bin_quality "$checkm2_report" "${base_bin}.orig")
        local orig_exists=$?

        # Note: bin without suffix is the same as .orig
        if [ $orig_exists -ne 0 ]; then
            orig_quality=$(get_bin_quality "$checkm2_report" "${base_bin}")
            orig_exists=$?
        fi

        local strict_quality=$(get_bin_quality "$checkm2_report" "${base_bin}.strict")
        local strict_exists=$?

        local permissive_quality=$(get_bin_quality "$checkm2_report" "${base_bin}.permissive")
        local permissive_exists=$?

        # Ensure at least one version exists
        if [ $orig_exists -ne 0 ] && [ $strict_exists -ne 0 ] && [ $permissive_exists -ne 0 ]; then
            log "  WARNING: No versions found for $base_bin, skipping..."
            continue
        fi

        # Parse quality metrics
        local best_version=""
        local best_comp=0
        local best_cont=100
        local best_score=-500
        local reason=""

        # Evaluate orig version
        if [ $orig_exists -eq 0 ]; then
            IFS='|' read -r comp cont <<< "$orig_quality"
            local score=$(echo "scale=2; $comp - (5 * $cont)" | bc -l)
            log "  orig: completeness=$comp%, contamination=$cont%, score=$score"

            best_version="orig"
            best_comp="$comp"
            best_cont="$cont"
            best_score="$score"
            reason="Only version available"
        fi

        # Evaluate strict version
        if [ $strict_exists -eq 0 ]; then
            IFS='|' read -r comp cont <<< "$strict_quality"
            local score=$(echo "scale=2; $comp - (5 * $cont)" | bc -l)
            log "  strict: completeness=$comp%, contamination=$cont%, score=$score"

            if [ -z "$best_version" ] || compare_bin_quality "$comp" "$cont" "$best_comp" "$best_cont"; then
                best_version="strict"
                best_comp="$comp"
                best_cont="$cont"
                best_score="$score"
                reason="Better quality score than other versions"
            fi
        fi

        # Evaluate permissive version
        if [ $permissive_exists -eq 0 ]; then
            IFS='|' read -r comp cont <<< "$permissive_quality"
            local score=$(echo "scale=2; $comp - (5 * $cont)" | bc -l)
            log "  permissive: completeness=$comp%, contamination=$cont%, score=$score"

            if [ -z "$best_version" ] || compare_bin_quality "$comp" "$cont" "$best_comp" "$best_cont"; then
                best_version="permissive"
                best_comp="$comp"
                best_cont="$cont"
                best_score="$score"
                reason="Better quality score than other versions"
            fi
        fi

        log "  SELECTED: $best_version (score=$best_score)"

        # Record selection
        echo -e "${base_bin}\t${best_version}\t${best_comp}\t${best_cont}\t${best_score}\t${reason}" >> "$selection_report"

        # Copy selected bin to output directory
        local source_bin=""
        if [ "$best_version" = "orig" ]; then
            # Check for .orig suffix first, then without suffix
            if [ -n "$sample_name" ]; then
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/${base_bin}.orig.fa"
                if [ ! -f "$source_bin" ]; then
                    source_bin="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/${base_bin}.fa"
                fi
            else
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/purified_bins/${base_bin}.orig.fa"
                if [ ! -f "$source_bin" ]; then
                    source_bin="${OUTPUT_DIR}/magpurify/${treatment}/purified_bins/${base_bin}.fa"
                fi
            fi
            ((orig_selected++))
        elif [ "$best_version" = "strict" ]; then
            if [ -n "$sample_name" ]; then
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/${base_bin}.strict.fa"
            else
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/purified_bins/${base_bin}.strict.fa"
            fi
            ((strict_selected++))
        else
            if [ -n "$sample_name" ]; then
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins/${base_bin}.permissive.fa"
            else
                source_bin="${OUTPUT_DIR}/magpurify/${treatment}/purified_bins/${base_bin}.permissive.fa"
            fi
            ((permissive_selected++))
        fi

        if [ -f "$source_bin" ]; then
            cp "$source_bin" "${output_dir}/${base_bin}.fa"
            ((selected_count++))
            log "  Copied to: ${output_dir}/${base_bin}.fa"
        else
            log "  ERROR: Source bin not found: $source_bin"
        fi
    done

    # Add summary to report
    cat >> "$selection_report" << EOF

Selection Summary:
  Total bins processed: $total_bins
  Bins successfully selected: $selected_count

Version Distribution:
  Original versions selected: $orig_selected
  Strict versions selected: $strict_selected
  Permissive versions selected: $permissive_selected

Output Directory: $output_dir
EOF

    log "Bin selection completed for $entity_desc:"
    log "  Total processed: $total_bins bins"
    log "  Successfully selected: $selected_count bins"
    log "  Version distribution: orig=$orig_selected, strict=$strict_selected, permissive=$permissive_selected"

    # Create completion flag
    touch "${output_dir}/bin_selection_complete.flag"

    if [ $selected_count -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

# Validation function
validate_bin_selection() {
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}/${sample_name}"
    else
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}"
    fi

    # Check completion flag
    if [ ! -f "${output_dir}/bin_selection_complete.flag" ]; then
        return 1
    fi

    # Check that we have selected bins
    local selected_count=$(ls -1 "${output_dir}"/*.fa 2>/dev/null | wc -l)
    if [ $selected_count -eq 0 ]; then
        return 1
    fi

    # Check that selection report exists
    if [ ! -f "${output_dir}/bin_selection_report.txt" ]; then
        return 1
    fi

    log "Validation successful: Selected $selected_count bins"
    return 0
}

# Run the bin selection stage based on processing mode
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    # Treatment-level mode: run once for all bins in treatment
    if stage_bin_selection "$TREATMENT"; then
        # Validate results
        if validate_bin_selection "$TREATMENT"; then
            create_treatment_checkpoint "$TREATMENT" "bin_selection"
            log "====== Bin selection completed successfully for treatment $TREATMENT ======"
        else
            log "ERROR: Bin selection validation failed for treatment $TREATMENT"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: Bin selection stage failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

else
    # Sample-level mode: run for individual sample
    if stage_bin_selection "$TREATMENT" "$SAMPLE_NAME"; then
        # Validate results
        if validate_bin_selection "$TREATMENT" "$SAMPLE_NAME"; then
            create_sample_checkpoint "$SAMPLE_NAME" "bin_selection"
            log "====== Bin selection completed successfully for $SAMPLE_NAME ======"
        else
            log "ERROR: Bin selection validation failed for $SAMPLE_NAME"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: Bin selection stage failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
