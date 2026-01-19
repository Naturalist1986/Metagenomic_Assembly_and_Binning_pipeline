#!/bin/bash
#SBATCH --job-name=bin_selection
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

# 08_bin_selection.sh - Select high-quality bins from Binette consensus binning
# Note: Binette already runs CheckM2, so we just read its quality reports

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Determine processing mode: treatment-level or sample-level
# Treatment-level if TREATMENT_LEVEL_BINNING=true OR ASSEMBLY_MODE=coassembly
if [ "${TREATMENT_LEVEL_BINNING:-false}" = "true" ] || [ "${ASSEMBLY_MODE:-individual}" = "coassembly" ]; then
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

    # Check if stage already completed for treatment (skip if --force flag used)
    if [ "${FORCE_RUN:-false}" != "true" ] && check_treatment_checkpoint "$TREATMENT" "bin_selection"; then
        log "Bin selection already completed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    if [ "${FORCE_RUN:-false}" = "true" ]; then
        log "Force mode enabled: Re-running bin selection even if previously completed"
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

    # Check if stage already completed for sample (skip if --force flag used)
    if [ "${FORCE_RUN:-false}" != "true" ] && check_sample_checkpoint "$SAMPLE_NAME" "bin_selection"; then
        log "Bin selection already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    if [ "${FORCE_RUN:-false}" = "true" ]; then
        log "Force mode enabled: Re-running bin selection even if previously completed"
    fi
fi

# Quality thresholds for bin selection
MIN_COMPLETENESS=50
MAX_CONTAMINATION=10

# Main processing function
stage_bin_selection() {
    local treatment="$1"
    local sample_name="${2:-}"  # Optional for treatment-level mode

    if [ -n "$sample_name" ]; then
        # Sample-level mode
        log "Selecting high-quality bins for $sample_name ($treatment)"
        local binette_dir="${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}/binette"
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}/${sample_name}"
        local entity_desc="sample $sample_name"
    else
        # Treatment-level mode
        log "Selecting high-quality bins for treatment $treatment"
        local binette_dir="${OUTPUT_DIR}/bin_refinement/${treatment}/binette"
        local output_dir="${OUTPUT_DIR}/selected_bins/${treatment}"
        local entity_desc="treatment $treatment"
    fi

    mkdir -p "$output_dir"

    # Check if already processed (skip if --force flag used)
    if [ "${FORCE_RUN:-false}" != "true" ] && [ -f "${output_dir}/bin_selection_complete.flag" ]; then
        log "Bin selection already completed for $entity_desc, skipping..."
        return 0
    fi

    # Check for Binette quality report
    local binette_quality_report="${binette_dir}/final_bins_quality_reports.tsv"
    if [ ! -f "$binette_quality_report" ]; then
        log "ERROR: Binette quality report not found: $binette_quality_report"
        log "Binette must complete successfully before bin selection"
        return 1
    fi

    # Check for Binette final bins directory
    local binette_bins_dir="${binette_dir}/final_bins"
    if [ ! -d "$binette_bins_dir" ]; then
        log "ERROR: Binette bins directory not found: $binette_bins_dir"
        return 1
    fi

    log "Reading Binette quality report: $binette_quality_report"
    log "Source bins directory: $binette_bins_dir"

    # Create selection report
    local selection_report="${output_dir}/bin_selection_report.txt"
    cat > "$selection_report" << EOF
Bin Selection Report for $entity_desc
=====================================

Date: $(date)
Treatment: $treatment
$([ -n "$sample_name" ] && echo "Sample: $sample_name")

Source: Binette consensus binning with CheckM2 quality assessment

Selection Criteria:
  Minimum Completeness: ${MIN_COMPLETENESS}%
  Maximum Contamination: ${MAX_CONTAMINATION}%

Bin Selection Results:
Name\tCompleteness\tContamination\tN50\tSize\tStatus
EOF

    # Parse Binette quality report and select bins
    local total_bins=0
    local selected_bins=0
    local high_quality=0
    local medium_quality=0
    local low_quality=0

    # Skip header line and process each bin
    # Binette columns: name origin is_original original_name completeness contamination score checkm2_size N50 coding_density contig_count
    tail -n +2 "$binette_quality_report" | while IFS=$'\t' read -r name origin is_original original_name completeness contamination score checkm2_size n50 coding_density contig_count; do
        ((total_bins++))

        # Use checkm2_size as size for reporting
        size="$checkm2_size"
        contigs="$contig_count"

        log "Processing bin: $name (completeness=${completeness}%, contamination=${contamination}%)"

        # Determine quality tier
        local quality_tier=""
        local status=""

        if (( $(echo "$completeness >= 90" | bc -l) )) && (( $(echo "$contamination < 5" | bc -l) )); then
            quality_tier="high"
            status="SELECTED (High Quality)"
            ((high_quality++))
            ((selected_bins++))
        elif (( $(echo "$completeness >= ${MIN_COMPLETENESS}" | bc -l) )) && (( $(echo "$contamination <= ${MAX_CONTAMINATION}" | bc -l) )); then
            quality_tier="medium"
            status="SELECTED (Medium Quality)"
            ((medium_quality++))
            ((selected_bins++))
        else
            quality_tier="low"
            status="REJECTED (Below Quality Threshold)"
            ((low_quality++))
        fi

        log "  $status"

        # Record in report
        echo -e "${name}\t${completeness}\t${contamination}\t${n50}\t${size}\t${status}" >> "$selection_report"

        # Copy selected bin to output directory
        if [ "$quality_tier" != "low" ]; then
            local source_bin="${binette_bins_dir}/${name}.fa"
            local dest_bin="${output_dir}/${name}.fa"

            if [ -f "$source_bin" ]; then
                cp "$source_bin" "$dest_bin"
                log "  Copied: $source_bin -> $dest_bin"
            else
                log "  WARNING: Bin file not found: $source_bin"
            fi
        fi
    done

    # Count actual bins copied
    local bins_copied=$(find "$output_dir" -name "*.fa" 2>/dev/null | wc -l)

    # Add summary to report
    cat >> "$selection_report" << EOF

=====================================
Summary:
  Total bins from Binette: $total_bins
  High Quality (>=90% comp, <5% cont): $high_quality
  Medium Quality (>=${MIN_COMPLETENESS}% comp, <=${MAX_CONTAMINATION}% cont): $medium_quality
  Low Quality (rejected): $low_quality
  Total selected: $selected_bins
  Bins copied: $bins_copied
=====================================
EOF

    log "====== Bin Selection Summary ======"
    log "Total bins from Binette: $total_bins"
    log "High Quality: $high_quality"
    log "Medium Quality: $medium_quality"
    log "Low Quality (rejected): $low_quality"
    log "Total selected: $selected_bins"
    log "Bins copied: $bins_copied"
    log "Selection report: $selection_report"

    # Create completion flag
    touch "${output_dir}/bin_selection_complete.flag"

    log "✓ Bin selection completed successfully for $entity_desc"
    return 0
}

# Execute main function
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    stage_bin_selection "$TREATMENT"
    status=$?

    if [ $status -eq 0 ]; then
        create_treatment_checkpoint "$TREATMENT" "bin_selection"
    fi

    cleanup_temp_dir "$TEMP_DIR"
    exit $status
else
    stage_bin_selection "$TREATMENT" "$SAMPLE_NAME"
    status=$?

    if [ $status -eq 0 ]; then
        create_sample_checkpoint "$SAMPLE_NAME" "bin_selection"
    fi

    cleanup_temp_dir "$TEMP_DIR"
    exit $status
fi
