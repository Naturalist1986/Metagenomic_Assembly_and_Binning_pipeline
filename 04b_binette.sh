#!/usr/bin/env bash
#SBATCH --job-name=binette
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --account=ofinkel

# 04b_binette.sh - Consensus binning using Binette
# Combines ALL 5 binners' results to produce high-quality consensus bins
# Always uses: MetaBAT2, MaxBin2, CONCOCT, COMEBin, and SemiBin
# Supports both treatment-level (coassembly) and sample-level (individual assembly) modes

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Set up temporary directory
TEMP_DIR=$(setup_temp_dir "binette")

# ===== FUNCTION DEFINITIONS =====

# Run Binette consensus binning
run_binette() {
    local contigs_file="$1"
    local output_dir="$2"
    shift 2
    local bin_dirs=("$@")  # Array of bin directories
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
    local valid_count=0
    for bin_dir in "${bin_dirs[@]}"; do
        if [ ! -d "$bin_dir" ]; then
            log "  WARNING: Skipping missing bin directory: $bin_dir"
            continue
        fi

        # Check for both .fa and .fa.gz files
        local fa_count=$(ls -1 "${bin_dir}"/*.fa 2>/dev/null | wc -l)
        local fagz_count=$(ls -1 "${bin_dir}"/*.fa.gz 2>/dev/null | wc -l)
        local total_bins=$((fa_count + fagz_count))

        if [ $total_bins -eq 0 ]; then
            log "  WARNING: Skipping empty bin directory: $bin_dir"
            continue
        fi

        # If there are .fa.gz files, decompress them first
        if [ $fagz_count -gt 0 ]; then
            log "  Decompressing $fagz_count gzipped bins in $bin_dir..."
            for gz_file in "${bin_dir}"/*.fa.gz; do
                if [ -f "$gz_file" ]; then
                    gunzip -f "$gz_file" 2>&1 | head -5
                fi
            done
            # Update counts after decompression
            fa_count=$(ls -1 "${bin_dir}"/*.fa 2>/dev/null | wc -l)
        fi

        log "  Including bin set: $bin_dir ($fa_count bins)"
        bin_dirs_args="$bin_dirs_args --bin_dirs $bin_dir"
        ((valid_count++))
    done

    if [ -z "$bin_dirs_args" ]; then
        log "ERROR: No valid bin directories found for Binette"
        conda deactivate
        return 1
    fi

    log "Using $valid_count bin sets for consensus binning"

    mkdir -p "$output_dir"

    # Validate contig IDs match between bins and assembly
    log "Validating contig IDs between bins and assembly..."

    # Extract contig IDs from assembly (first field only, no descriptions)
    grep "^>" "$contigs_file" | sed 's/^>//' | cut -f1 -d' ' | sort -u > "${TEMP_DIR}/assembly_contigs.txt"
    local assembly_contig_count=$(wc -l < "${TEMP_DIR}/assembly_contigs.txt")
    log "  Assembly has $assembly_contig_count unique contigs"

    # Check a sample bin for contig ID format
    local sample_bin=$(find "${bin_dirs[@]}" -name "*.fa" -type f 2>/dev/null | head -1)
    if [ -f "$sample_bin" ]; then
        log "  Checking contig format in: $(basename $sample_bin)"
        grep "^>" "$sample_bin" | head -3 | sed 's/^>//' > "${TEMP_DIR}/sample_bin_contigs.txt"

        # Check if any of the sample contigs are in the assembly
        local matches=0
        while IFS= read -r contig_id; do
            # Try exact match first
            if grep -qFx "$contig_id" "${TEMP_DIR}/assembly_contigs.txt"; then
                ((matches++))
            else
                # Try matching first field only (strip after first space)
                local contig_base=$(echo "$contig_id" | cut -f1 -d' ')
                if grep -qFx "$contig_base" "${TEMP_DIR}/assembly_contigs.txt"; then
                    ((matches++))
                fi
            fi
        done < "${TEMP_DIR}/sample_bin_contigs.txt"

        log "  Sample check: $matches/3 contigs found in assembly"

        if [ $matches -eq 0 ]; then
            log "  WARNING: Contig ID format mismatch detected!"
            log "  Sample bin contig:"
            head -1 "${TEMP_DIR}/sample_bin_contigs.txt"
            log "  Sample assembly contig:"
            head -1 "${TEMP_DIR}/assembly_contigs.txt"
            log "  Attempting to normalize contig IDs in bins..."

            # Create normalized bin directories
            local normalized_dirs=()
            for bin_dir in "${bin_dirs[@]}"; do
                local norm_dir="${TEMP_DIR}/normalized_bins/$(basename $bin_dir)"
                mkdir -p "$norm_dir"

                # Normalize each bin file (keep only first field of contig IDs)
                for bin_file in "${bin_dir}"/*.fa; do
                    if [ -f "$bin_file" ]; then
                        local norm_file="${norm_dir}/$(basename $bin_file)"

                        # Use sed to normalize contig IDs (faster and simpler than Python)
                        sed '/^>/ s/\s.*//' "$bin_file" > "$norm_file"
                    fi
                done

                normalized_dirs+=("$norm_dir")
            done

            # Use normalized directories instead
            log "  Created normalized bins in: ${TEMP_DIR}/normalized_bins/"
            bin_dirs=("${normalized_dirs[@]}")

            # Rebuild bin_dirs_args with normalized directories
            bin_dirs_args=""
            for norm_dir in "${normalized_dirs[@]}"; do
                local bin_count=$(ls -1 "${norm_dir}"/*.fa 2>/dev/null | wc -l)
                if [ $bin_count -gt 0 ]; then
                    log "    Normalized: $(basename $norm_dir) ($bin_count bins)"
                    bin_dirs_args="$bin_dirs_args --bin_dirs $norm_dir"
                fi
            done
        fi
    fi

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

    log "Running Binette command:"
    log "$binette_cmd"

    # Run Binette and capture output (tee shows output in real-time and saves to log)
    local binette_log="${LOG_DIR}/${TREATMENT}_binette.log"
    mkdir -p "${LOG_DIR}"

    $binette_cmd 2>&1 | tee "$binette_log"

    local exit_code=${PIPESTATUS[0]}

    # Show summary of output for quick diagnosis
    log ""
    log "=== Binette Execution Summary ==="
    if [ -f "$binette_log" ]; then
        local log_lines=$(wc -l < "$binette_log")
        log "Full log: $binette_log ($log_lines lines)"

        # Show key errors or warnings
        if grep -qi "error\|traceback\|exception" "$binette_log"; then
            log "Errors detected in log:"
            grep -i "error\|traceback\|exception" "$binette_log" | tail -10 | while IFS= read -r line; do
                log "  ERROR: $line"
            done
        fi
    fi

    conda deactivate

    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/final_bins" ]; then
        local bin_count=$(ls -1 "${output_dir}/final_bins"/*.fa 2>/dev/null | wc -l)
        log "✓ Binette completed: $bin_count consensus bins"
        return 0
    else
        log "✗ Binette failed with exit code: $exit_code"

        # Show last 50 lines for full context
        if [ -f "$binette_log" ]; then
            log ""
            log "=== Last 50 lines of Binette output ==="
            tail -50 "$binette_log" | while IFS= read -r line; do
                log "  $line"
            done
        else
            log "WARNING: Binette log file not found at: $binette_log"
        fi

        return 1
    fi
}

# Create Binette summary
create_binette_summary() {
    local binette_dir="$1"
    local identifier="$2"  # treatment or sample name
    local summary_file="${binette_dir}/binette_summary.txt"

    cat > "$summary_file" << EOF
Binette Consensus Binning Summary
==================================

Date: $(date)
Identifier: $identifier
Tool: Binette
Bin sets used: MetaBAT2, MaxBin2, CONCOCT, COMEBin, SemiBin

Results:
EOF

    if [ -d "${binette_dir}/final_bins" ]; then
        local total_bins=$(ls -1 "${binette_dir}/final_bins"/*.fa 2>/dev/null | wc -l)
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

# ===== DETERMINE MODE =====

if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
    # ===== TREATMENT-LEVEL MODE (co-assembly) =====
    log "Running in TREATMENT-LEVEL mode (co-assembly)"

    TREATMENTS_ARRAY=($(get_treatments))
    TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

    if [ $TASK_ID -ge ${#TREATMENTS_ARRAY[@]} ]; then
        log "Array task $TASK_ID has no treatment to process"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    TREATMENT="${TREATMENTS_ARRAY[$TASK_ID]}"
    export TREATMENT

    log "====== Starting Binette Consensus Binning for Treatment: $TREATMENT ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"
    binette_dir="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/binette"

    # Check if already processed
    if [ -d "${binette_dir}/final_bins" ] && [ "$(ls -A ${binette_dir}/final_bins/*.fa 2>/dev/null)" ]; then
        log "Binette already completed for treatment $TREATMENT, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly found in $assembly_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Assembly: $assembly_fasta"

    # Collect ALL 5 bin directories
    bin_dirs_array=()

    log "Collecting bin directories from all 5 binners..."

    # 1. MetaBAT2 bins
    if [ -d "${binning_dir}/metabat2_bins" ]; then
        bin_dirs_array+=("${binning_dir}/metabat2_bins")
        log "  Found MetaBAT2 bins"
    else
        log "  WARNING: MetaBAT2 bins not found at ${binning_dir}/metabat2_bins"
    fi

    # 2. MaxBin2 bins
    if [ -d "${binning_dir}/maxbin2_bins" ]; then
        bin_dirs_array+=("${binning_dir}/maxbin2_bins")
        log "  Found MaxBin2 bins"
    else
        log "  WARNING: MaxBin2 bins not found at ${binning_dir}/maxbin2_bins"
    fi

    # 3. CONCOCT bins
    if [ -d "${binning_dir}/concoct_bins" ]; then
        bin_dirs_array+=("${binning_dir}/concoct_bins")
        log "  Found CONCOCT bins"
    else
        log "  WARNING: CONCOCT bins not found at ${binning_dir}/concoct_bins"
    fi

    # 4. COMEBin bins
    if [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
        bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
        log "  Found COMEBin bins"
    elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
        log "  Found COMEBin bins (legacy path)"
    else
        log "  WARNING: COMEBin bins not found"
    fi

    # 5. SemiBin bins
    if [ -d "${binning_dir}/semibin/output_bins" ]; then
        bin_dirs_array+=("${binning_dir}/semibin/output_bins")
        log "  Found SemiBin bins"
    else
        log "  WARNING: SemiBin bins not found at ${binning_dir}/semibin/output_bins"
    fi

    if [ ${#bin_dirs_array[@]} -eq 0 ]; then
        log "ERROR: No bin directories found from any of the 5 binners"
        log "Expected locations:"
        log "  - ${binning_dir}/metabat2_bins"
        log "  - ${binning_dir}/maxbin2_bins"
        log "  - ${binning_dir}/concoct_bins"
        log "  - ${binning_dir}/comebin/comebin_res/comebin_res_bins"
        log "  - ${binning_dir}/semibin/output_bins"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Found ${#bin_dirs_array[@]} bin sets for consensus binning"

    # Run Binette with all available bin sets
    if run_binette "$assembly_fasta" "$binette_dir" "${bin_dirs_array[@]}"; then
        create_binette_summary "$binette_dir" "$TREATMENT"
        log "====== Binette completed for treatment $TREATMENT ======"
    else
        log "ERROR: Binette failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

else
    # ===== SAMPLE-LEVEL MODE (individual assembly) =====
    log "Running in SAMPLE-LEVEL mode (individual assembly)"

    SAMPLE_INFO=$(get_sample_info_by_index "$SLURM_ARRAY_TASK_ID")
    if [ -z "$SAMPLE_INFO" ]; then
        log "ERROR: No sample found for array index $SLURM_ARRAY_TASK_ID"
        exit 1
    fi

    IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"
    export SAMPLE_NAME TREATMENT

    create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"

    log "====== Starting Binette Consensus Binning for $SAMPLE_NAME ($TREATMENT) ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"
    binette_dir="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/${SAMPLE_NAME}/binette"

    # Check if already processed
    if [ -d "${binette_dir}/final_bins" ] && [ "$(ls -A ${binette_dir}/final_bins/*.fa 2>/dev/null)" ]; then
        log "Binette already completed for $SAMPLE_NAME, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/${SAMPLE_NAME}_contigs.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly found in $assembly_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Assembly: $assembly_fasta"

    # Collect ALL 5 bin directories
    bin_dirs_array=()

    log "Collecting bin directories from all 5 binners..."

    # 1. MetaBAT2 bins
    if [ -d "${binning_dir}/metabat2_bins" ]; then
        bin_dirs_array+=("${binning_dir}/metabat2_bins")
        log "  Found MetaBAT2 bins"
    else
        log "  WARNING: MetaBAT2 bins not found at ${binning_dir}/metabat2_bins"
    fi

    # 2. MaxBin2 bins
    if [ -d "${binning_dir}/maxbin2_bins" ]; then
        bin_dirs_array+=("${binning_dir}/maxbin2_bins")
        log "  Found MaxBin2 bins"
    else
        log "  WARNING: MaxBin2 bins not found at ${binning_dir}/maxbin2_bins"
    fi

    # 3. CONCOCT bins
    if [ -d "${binning_dir}/concoct_bins" ]; then
        bin_dirs_array+=("${binning_dir}/concoct_bins")
        log "  Found CONCOCT bins"
    else
        log "  WARNING: CONCOCT bins not found at ${binning_dir}/concoct_bins"
    fi

    # 4. COMEBin bins
    if [ -d "${binning_dir}/comebin/comebin_res/comebin_res_bins" ]; then
        bin_dirs_array+=("${binning_dir}/comebin/comebin_res/comebin_res_bins")
        log "  Found COMEBin bins"
    elif [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        bin_dirs_array+=("${binning_dir}/comebin/comebin_res_bins")
        log "  Found COMEBin bins (legacy path)"
    else
        log "  WARNING: COMEBin bins not found"
    fi

    # 5. SemiBin bins
    if [ -d "${binning_dir}/semibin/output_bins" ]; then
        bin_dirs_array+=("${binning_dir}/semibin/output_bins")
        log "  Found SemiBin bins"
    else
        log "  WARNING: SemiBin bins not found at ${binning_dir}/semibin/output_bins"
    fi

    if [ ${#bin_dirs_array[@]} -eq 0 ]; then
        log "ERROR: No bin directories found from any of the 5 binners"
        log "Expected locations:"
        log "  - ${binning_dir}/metabat2_bins"
        log "  - ${binning_dir}/maxbin2_bins"
        log "  - ${binning_dir}/concoct_bins"
        log "  - ${binning_dir}/comebin/comebin_res/comebin_res_bins"
        log "  - ${binning_dir}/semibin/output_bins"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Found ${#bin_dirs_array[@]} bin sets for consensus binning"

    # Run Binette with all available bin sets
    if run_binette "$assembly_fasta" "$binette_dir" "${bin_dirs_array[@]}"; then
        create_binette_summary "$binette_dir" "$SAMPLE_NAME"
        log "====== Binette completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Binette failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

cleanup_temp_dir "$TEMP_DIR"
log "====== Binette completed successfully ======"
