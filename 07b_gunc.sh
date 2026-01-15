#!/bin/bash
#SBATCH --job-name=gunc
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel

# 07b_gunc.sh - GUNC chimerism detection
# Detects chimeric/contaminated bins after quality assessment
# Runs after CheckM2 (stage 7)

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# GUNC database path
GUNC_DB="/sci/backup/ofinkel/moshea/Databases/gunc_db/gunc_db_progenomes2.1.dmnd"

# ===== FUNCTION DEFINITIONS =====

# Run GUNC chimerism detection
run_gunc() {
    local bins_dir="$1"
    local output_dir="$2"

    log "Running GUNC chimerism detection..."
    log "Bins directory: $bins_dir"
    log "Output directory: $output_dir"

    # Activate GUNC environment
    activate_env gunc

    # Check if GUNC is available
    if ! command -v gunc &> /dev/null; then
        log "ERROR: GUNC not available in conda environment"
        conda deactivate
        return 1
    fi

    # Check database
    if [ ! -f "$GUNC_DB" ]; then
        log "ERROR: GUNC database not found at: $GUNC_DB"
        conda deactivate
        return 1
    fi

    mkdir -p "$output_dir"

    # Run GUNC
    gunc run \
        --input_dir "$bins_dir" \
        --db_file "$GUNC_DB" \
        --out_dir "$output_dir" \
        --threads ${SLURM_CPUS_PER_TASK:-32} \
        --detailed_output \
        --file_suffix .fa \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_gunc.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv" ]; then
        log "GUNC completed successfully"
        return 0
    else
        log "GUNC failed with exit code: $exit_code"
        return 1
    fi
}

# Create GUNC summary report
create_gunc_summary() {
    local gunc_dir="$1"
    local summary_file="${gunc_dir}/gunc_summary.txt"
    local gunc_results="${gunc_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv"

    if [ ! -f "$gunc_results" ]; then
        log "ERROR: GUNC results file not found"
        return 1
    fi

    cat > "$summary_file" << EOF
GUNC Chimerism Detection Summary
=================================

Date: $(date)
Sample: ${SAMPLE_NAME:-$TREATMENT}
Treatment: $TREATMENT
Tool: GUNC
Database: ProGenomes 2.1

Chimerism Detection Criteria:
  - clade_separation_score (CSS) > 0.45 = likely chimeric
  - CSS 0.30-0.45 = potential chimerism (review recommended)
  - CSS < 0.30 = likely not chimeric

Results:
EOF

    local total_bins=$(awk 'NR>1' "$gunc_results" | wc -l)
    local chimeric_bins=$(awk -F'\t' 'NR>1 && $NF > 0.45' "$gunc_results" | wc -l)
    local suspect_bins=$(awk -F'\t' 'NR>1 && $NF > 0.30 && $NF <= 0.45' "$gunc_results" | wc -l)
    local clean_bins=$(awk -F'\t' 'NR>1 && $NF <= 0.30' "$gunc_results" | wc -l)

    echo "  Total bins analyzed: $total_bins" >> "$summary_file"
    echo "  Likely chimeric (CSS > 0.45): $chimeric_bins" >> "$summary_file"
    echo "  Potential chimerism (CSS 0.30-0.45): $suspect_bins" >> "$summary_file"
    echo "  Likely clean (CSS ≤ 0.30): $clean_bins" >> "$summary_file"
    echo "" >> "$summary_file"

    # List chimeric bins
    if [ $chimeric_bins -gt 0 ]; then
        echo "Likely Chimeric Bins:" >> "$summary_file"
        awk -F'\t' 'NR>1 && $NF > 0.45 {
            printf "  %s: CSS = %.3f\n", $1, $NF
        }' "$gunc_results" >> "$summary_file"
        echo "" >> "$summary_file"
    fi

    # List suspect bins
    if [ $suspect_bins -gt 0 ]; then
        echo "Potential Chimeric Bins (Review Recommended):" >> "$summary_file"
        awk -F'\t' 'NR>1 && $NF > 0.30 && $NF <= 0.45 {
            printf "  %s: CSS = %.3f\n", $1, $NF
        }' "$gunc_results" >> "$summary_file"
        echo "" >> "$summary_file"
    fi

    echo "Full results available in: $gunc_results" >> "$summary_file"
    echo "" >> "$summary_file"
    echo "Note: Review detailed output for per-contig taxonomic assignments" >> "$summary_file"
    echo "      to understand the source of chimerism (e.g., plasmid contamination," >> "$summary_file"
    echo "      strain mixing, or cross-species contamination)" >> "$summary_file"

    log "GUNC summary created: $summary_file"
}

# ===== MAIN EXECUTION =====

# Initialize
init_conda
TEMP_DIR=$(setup_temp_dir)

# Determine if we're in treatment-level or sample-level mode
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ "${ASSEMBLY_MODE}" = "coassembly" ]; then
    # ===== TREATMENT-LEVEL MODE (coassembly) =====
    log "Running in TREATMENT-LEVEL mode (coassembly)"

    if [ ! -f "$TREATMENTS_FILE" ]; then
        echo "ERROR: Treatments file not found at: $TREATMENTS_FILE"
        exit 1
    fi

    mapfile -t TREATMENTS_ARRAY < "$TREATMENTS_FILE"

    if [ $TASK_ID -ge ${#TREATMENTS_ARRAY[@]} ]; then
        echo "No treatment found for array index $TASK_ID"
        exit 0
    fi

    TREATMENT="${TREATMENTS_ARRAY[$TASK_ID]}"
    export TREATMENT

    log "====== Starting GUNC Chimerism Detection for Treatment: $TREATMENT ======"

    # Set up directories
    local binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    local gunc_dir="${binning_dir}/gunc"

    # Check if already processed
    if [ -f "${gunc_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv" ]; then
        log "GUNC already completed for treatment $TREATMENT, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find bins to analyze - use best available set
    local bins_dir=""

    # Priority: Binette > BinSPreader > Bin Refinement > Initial Binning
    if [ "$USE_BINETTE" = "true" ] && [ -d "${binning_dir}/binette/final_bins" ]; then
        bins_dir="${binning_dir}/binette/final_bins"
    elif [ "$USE_BINSPREADER" = "true" ] && [ -d "${binning_dir}/binspreader/bins" ]; then
        bins_dir="${binning_dir}/binspreader/bins"
    elif [ -d "${binning_dir}/bin_refinement/metawrap_bins" ]; then
        bins_dir="${binning_dir}/bin_refinement/metawrap_bins"
    elif [ "$USE_COMEBIN" = "true" ] && [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        bins_dir="${binning_dir}/comebin/comebin_res_bins"
    elif [ "$USE_SEMIBIN" = "true" ] && [ -d "${binning_dir}/semibin/output_bins" ]; then
        bins_dir="${binning_dir}/semibin/output_bins"
    elif [ -d "${binning_dir}/metawrap_bins" ]; then
        bins_dir="${binning_dir}/metawrap_bins"
    else
        log "ERROR: No bins found for GUNC analysis"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Using bins from: $bins_dir"

    # Run GUNC
    if run_gunc "$bins_dir" "$gunc_dir"; then
        create_gunc_summary "$gunc_dir"
        log "====== GUNC completed for treatment $TREATMENT ======"
    else
        log "ERROR: GUNC failed"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

else
    # ===== SAMPLE-LEVEL MODE (individual assembly) =====
    log "Running in SAMPLE-LEVEL mode (individual assembly)"

    SAMPLE_INFO=$(get_sample_info_by_index "$TASK_ID")
    if [ -z "$SAMPLE_INFO" ]; then
        echo "ERROR: No sample found for array index $TASK_ID"
        exit 1
    fi

    IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"
    export SAMPLE_NAME TREATMENT

    create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"

    log "====== Starting GUNC Chimerism Detection for $SAMPLE_NAME ($TREATMENT) ======"

    # Set up directories
    local binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    local gunc_dir="${binning_dir}/gunc"

    # Check if already processed
    if [ -f "${gunc_dir}/GUNC.progenomes_2.1.maxCSS_level.tsv" ]; then
        log "GUNC already completed for $SAMPLE_NAME, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find bins to analyze - use best available set
    local bins_dir=""

    # Priority: Binette > BinSPreader > Bin Refinement > Initial Binning
    if [ "$USE_BINETTE" = "true" ] && [ -d "${binning_dir}/binette/final_bins" ]; then
        bins_dir="${binning_dir}/binette/final_bins"
    elif [ "$USE_BINSPREADER" = "true" ] && [ -d "${binning_dir}/binspreader/bins" ]; then
        bins_dir="${binning_dir}/binspreader/bins"
    elif [ -d "${binning_dir}/bin_refinement/metawrap_bins" ]; then
        bins_dir="${binning_dir}/bin_refinement/metawrap_bins"
    elif [ "$USE_COMEBIN" = "true" ] && [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        bins_dir="${binning_dir}/comebin/comebin_res_bins"
    elif [ "$USE_SEMIBIN" = "true" ] && [ -d "${binning_dir}/semibin/output_bins" ]; then
        bins_dir="${binning_dir}/semibin/output_bins"
    elif [ -d "${binning_dir}/metawrap_50_10_bins" ]; then
        bins_dir="${binning_dir}/metawrap_50_10_bins"
    else
        log "ERROR: No bins found for GUNC analysis"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Using bins from: $bins_dir"

    # Run GUNC
    if run_gunc "$bins_dir" "$gunc_dir"; then
        create_gunc_summary "$gunc_dir"
        log "====== GUNC completed for $SAMPLE_NAME ======"
    else
        log "ERROR: GUNC failed"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
