#!/bin/bash
#SBATCH --job-name=binette
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --account=ofinkel

# 04b_binette.sh - Consensus binning using Binette
# Combines multiple binning results to produce high-quality consensus bins
# Replaces DAS Tool when --binette flag is used

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# ===== FUNCTION DEFINITIONS =====

# Run Binette consensus binning
run_binette() {
    local contigs_file="$1"
    shift
    local bin_dirs=("$@")  # Array of bin directories
    local output_dir="${BINETTE_OUTPUT_DIR}"  # Use global variable for output directory
    local checkm2_db="${CHECKM2_DB:-}"  # Optional: specify CheckM2 database

    log "Running Binette consensus binning..."
    log "Contigs file: $contigs_file"
    log "Number of bin sets: ${#bin_dirs[@]}"

    # Activate binette environment
    activate_env binette

    # Check if Binette is available
    if ! command -v binette &> /dev/null; then
        log "ERROR: Binette not available in conda environment"
        conda deactivate
        return 1
    fi

    # Build bin_dirs arguments
    local bin_dirs_args=""
    for bin_dir in "${bin_dirs[@]}"; do
        if [ -d "$bin_dir" ] && [ "$(ls -A ${bin_dir}/*.fa 2>/dev/null)" ]; then
            log "  Including bin set: $bin_dir"
            bin_dirs_args="$bin_dirs_args --bin_dirs $bin_dir"
        else
            log "  WARNING: Skipping empty/missing bin directory: $bin_dir"
        fi
    done

    if [ -z "$bin_dirs_args" ]; then
        log "ERROR: No valid bin directories found for Binette"
        conda deactivate
        return 1
    fi

    mkdir -p "$output_dir"

    # Build Binette command
    local binette_cmd="binette \
        $bin_dirs_args \
        --contigs $contigs_file \
        --outdir $output_dir \
        --threads ${SLURM_CPUS_PER_TASK:-32} \
        --min_completeness 50 \
        --max_contamination 10"

    # Add CheckM2 database if specified
    if [ -n "$checkm2_db" ] && [ -f "$checkm2_db" ]; then
        binette_cmd="$binette_cmd --checkm2_db $checkm2_db"
    fi

    # Add low memory flag if needed (for >500 bins)
    # binette_cmd="$binette_cmd --low_mem"

    log "Running Binette command:"
    log "$binette_cmd"

    $binette_cmd 2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_binette.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/final_bins" ]; then
        local bin_count=$(ls -1 "${output_dir}/final_bins"/*.fa 2>/dev/null | wc -l)
        log "Binette completed successfully with $bin_count consensus bins"
        return 0
    else
        log "Binette failed with exit code: $exit_code"
        return 1
    fi
}

# Create Binette summary
create_binette_summary() {
    local binette_dir="$1"
    local summary_file="${binette_dir}/binette_summary.txt"

    cat > "$summary_file" << EOF
Binette Consensus Binning Summary
==================================

Date: $(date)
Sample: ${SAMPLE_NAME:-$TREATMENT}
Treatment: $TREATMENT
Tool: Binette

Results:
EOF

    local total_bins=0
    if [ -d "${binette_dir}/final_bins" ]; then
        total_bins=$(ls -1 "${binette_dir}/final_bins"/*.fa 2>/dev/null | wc -l)
        echo "  Total consensus bins: $total_bins" >> "$summary_file"
        echo "" >> "$summary_file"

        # Add bin quality info if CheckM2 results exist
        if [ -f "${binette_dir}/checkm2_results/quality_report.tsv" ]; then
            echo "Quality Metrics:" >> "$summary_file"
            awk 'NR>1 {
                comp_sum+=$2; cont_sum+=$3; count++
            } END {
                if (count>0) {
                    printf "  Average Completeness: %.1f%%\n", comp_sum/count
                    printf "  Average Contamination: %.1f%%\n", cont_sum/count
                }
            }' "${binette_dir}/checkm2_results/quality_report.tsv" >> "$summary_file"

            # Count high-quality bins
            local hq_bins=$(awk 'NR>1 && $2>=90 && $3<5' "${binette_dir}/checkm2_results/quality_report.tsv" | wc -l)
            local mq_bins=$(awk 'NR>1 && $2>=50 && $3<10 && ($2<90 || $3>=5)' "${binette_dir}/checkm2_results/quality_report.tsv" | wc -l)
            echo "  High Quality (≥90% complete, <5% contam): $hq_bins" >> "$summary_file"
            echo "  Medium Quality (≥50% complete, <10% contam): $mq_bins" >> "$summary_file"
        fi
    else
        echo "  No consensus bins produced" >> "$summary_file"
    fi

    log "Binette summary created: $summary_file"
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

    log "====== Starting Binette Consensus Binning for Treatment: $TREATMENT ======"

    # Set up directories
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"

    # Set output directory for Binette
    BINETTE_OUTPUT_DIR="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/binette"

    # Check if already processed
    if [ -d "${BINETTE_OUTPUT_DIR}/final_bins" ] && [ "$(ls -A ${BINETTE_OUTPUT_DIR}/final_bins/*.fa 2>/dev/null)" ]; then
        log "Binette already completed for treatment $TREATMENT, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly FASTA
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/assembly.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            log "Found assembly: $assembly_fasta"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly FASTA file found"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Collect bin directories from all sources
    bin_dirs_array=()

    # Decide whether to use refined or original bins
    if [ "$BINETTE_USE_REFINED" = "true" ]; then
        log "Using BinSPreader-refined bins for consensus binning"
        log "USE_COMEBIN=${USE_COMEBIN}, USE_SEMIBIN=${USE_SEMIBIN}, BINETTE_USE_REFINED=${BINETTE_USE_REFINED}"

        # Use refined COMEBin bins if available
        if [ "$USE_COMEBIN" = "true" ]; then
            if [ -d "${binning_dir}/binspreader/comebin_refined" ]; then
                log "Adding BinSPreader-refined COMEBin bins: ${binning_dir}/binspreader/comebin_refined"
                bin_dirs_array+=("${binning_dir}/binspreader/comebin_refined")
            elif [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
                log "Refined bins not found, using original COMEBin bins: ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
            elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
                log "Refined bins not found, using original COMEBin bins (legacy path): ${binning_dir}/comebin/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
            else
                log "COMEBin enabled but bins not found"
            fi
        fi

        # Use refined SemiBin bins if available
        if [ "$USE_SEMIBIN" = "true" ]; then
            if [ -d "${binning_dir}/binspreader/semibin_refined" ]; then
                log "Adding BinSPreader-refined SemiBin bins: ${binning_dir}/binspreader/semibin_refined"
                bin_dirs_array+=("${binning_dir}/binspreader/semibin_refined")
            elif [ -d "${binning_dir}/semibin/output_bins" ]; then
                log "Refined bins not found, using original SemiBin bins: ${binning_dir}/semibin/output_bins"
                bin_dirs_array+=("${binning_dir}/semibin/output_bins")
            else
                log "SemiBin enabled but bins not found"
            fi
        fi
    else
        log "Using original bins for consensus binning"
        log "USE_COMEBIN=${USE_COMEBIN}, USE_SEMIBIN=${USE_SEMIBIN}, BINETTE_USE_REFINED=${BINETTE_USE_REFINED}"

        if [ "$USE_COMEBIN" = "true" ]; then
            if [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
                log "Adding COMEBin bins: ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
            elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
                log "Adding COMEBin bins (legacy path): ${binning_dir}/comebin/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
            else
                log "COMEBin enabled but bins not found at:"
                log "  - ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                log "  - ${binning_dir}/comebin/comebin_res_bins"
            fi
        else
            log "COMEBin not enabled (USE_COMEBIN=${USE_COMEBIN})"
        fi

        if [ "$USE_SEMIBIN" = "true" ]; then
            if [ -d "${binning_dir}/semibin/output_bins" ]; then
                log "Adding SemiBin bins: ${binning_dir}/semibin/output_bins"
                bin_dirs_array+=("${binning_dir}/semibin/output_bins")
            else
                log "SemiBin enabled but bins not found at: ${binning_dir}/semibin/output_bins"
            fi
        else
            log "SemiBin not enabled (USE_SEMIBIN=${USE_SEMIBIN})"
        fi
    fi

    # Only use MetaWRAP if neither COMEBin nor SemiBin are enabled
    if [ "$USE_COMEBIN" != "true" ] && [ "$USE_SEMIBIN" != "true" ]; then
        if [ -d "${binning_dir}/metawrap_bins" ]; then
            log "Adding MetaWRAP bins"
            bin_dirs_array+=("${binning_dir}/metawrap_bins")
        fi
    fi

    if [ ${#bin_dirs_array[@]} -eq 0 ]; then
        log "ERROR: No bin directories found for consensus binning"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Run Binette
    if run_binette "$assembly_fasta" "${bin_dirs_array[@]}"; then
        create_binette_summary "${BINETTE_OUTPUT_DIR}"
        log "====== Binette completed for treatment $TREATMENT ======"
    else
        log "ERROR: Binette failed"
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

    log "====== Starting Binette Consensus Binning for $SAMPLE_NAME ($TREATMENT) ======"

    # Set up directories
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"

    # Set output directory for Binette
    BINETTE_OUTPUT_DIR="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/${SAMPLE_NAME}/binette"

    # Check if already processed
    if [ -d "${BINETTE_OUTPUT_DIR}/final_bins" ] && [ "$(ls -A ${BINETTE_OUTPUT_DIR}/final_bins/*.fa 2>/dev/null)" ]; then
        log "Binette already completed for $SAMPLE_NAME, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly FASTA
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/assembly.fasta" \
        "${assembly_dir}/${SAMPLE_NAME}_contigs.fasta" \
        "${assembly_dir}/${SAMPLE_NAME}_scaffolds.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            log "Found assembly: $assembly_fasta"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly FASTA file found"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Collect bin directories from all sources
    bin_dirs_array=()

    # Decide whether to use refined or original bins
    if [ "$BINETTE_USE_REFINED" = "true" ]; then
        log "Using BinSPreader-refined bins for consensus binning"
        log "USE_COMEBIN=${USE_COMEBIN}, USE_SEMIBIN=${USE_SEMIBIN}, BINETTE_USE_REFINED=${BINETTE_USE_REFINED}"

        # Use refined COMEBin bins if available
        if [ "$USE_COMEBIN" = "true" ]; then
            if [ -d "${binning_dir}/binspreader/comebin_refined" ]; then
                log "Adding BinSPreader-refined COMEBin bins: ${binning_dir}/binspreader/comebin_refined"
                bin_dirs_array+=("${binning_dir}/binspreader/comebin_refined")
            elif [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
                log "Refined bins not found, using original COMEBin bins: ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
            elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
                log "Refined bins not found, using original COMEBin bins (legacy path): ${binning_dir}/comebin/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
            else
                log "COMEBin enabled but bins not found"
            fi
        fi

        # Use refined SemiBin bins if available
        if [ "$USE_SEMIBIN" = "true" ]; then
            if [ -d "${binning_dir}/binspreader/semibin_refined" ]; then
                log "Adding BinSPreader-refined SemiBin bins: ${binning_dir}/binspreader/semibin_refined"
                bin_dirs_array+=("${binning_dir}/binspreader/semibin_refined")
            elif [ -d "${binning_dir}/semibin/output_bins" ]; then
                log "Refined bins not found, using original SemiBin bins: ${binning_dir}/semibin/output_bins"
                bin_dirs_array+=("${binning_dir}/semibin/output_bins")
            else
                log "SemiBin enabled but bins not found"
            fi
        fi
    else
        log "Using original bins for consensus binning"
        log "USE_COMEBIN=${USE_COMEBIN}, USE_SEMIBIN=${USE_SEMIBIN}, BINETTE_USE_REFINED=${BINETTE_USE_REFINED}"

        if [ "$USE_COMEBIN" = "true" ]; then
            if [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
                log "Adding COMEBin bins: ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
            elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
                log "Adding COMEBin bins (legacy path): ${binning_dir}/comebin/comebin_res_bins"
                bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
            else
                log "COMEBin enabled but bins not found at:"
                log "  - ${binning_dir}/comebin/comebin_res/comebin_res_bins"
                log "  - ${binning_dir}/comebin/comebin_res_bins"
            fi
        else
            log "COMEBin not enabled (USE_COMEBIN=${USE_COMEBIN})"
        fi

        if [ "$USE_SEMIBIN" = "true" ]; then
            if [ -d "${binning_dir}/semibin/output_bins" ]; then
                log "Adding SemiBin bins: ${binning_dir}/semibin/output_bins"
                bin_dirs_array+=("${binning_dir}/semibin/output_bins")
            else
                log "SemiBin enabled but bins not found at: ${binning_dir}/semibin/output_bins"
            fi
        else
            log "SemiBin not enabled (USE_SEMIBIN=${USE_SEMIBIN})"
        fi
    fi

    # Only use MetaWRAP if neither COMEBin nor SemiBin are enabled
    if [ "$USE_COMEBIN" != "true" ] && [ "$USE_SEMIBIN" != "true" ]; then
        if [ -d "${binning_dir}/metawrap_50_10_bins" ]; then
            log "Adding MetaWRAP bins"
            bin_dirs_array+=("${binning_dir}/metawrap_50_10_bins")
        fi
    fi

    if [ ${#bin_dirs_array[@]} -eq 0 ]; then
        log "ERROR: No bin directories found for consensus binning"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Run Binette
    if run_binette "$assembly_fasta" "${bin_dirs_array[@]}"; then
        create_binette_summary "${BINETTE_OUTPUT_DIR}"
        log "====== Binette completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Binette failed"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
