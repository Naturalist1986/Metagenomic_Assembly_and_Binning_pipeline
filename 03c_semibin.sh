#!/usr/bin/env bash
#SBATCH --job-name=semibin
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --account=ofinkel

# 03c_semibin.sh - Metagenomic binning using SemiBin2
# Uses shared BAM files created by 03a_create_shared_bams.sh
# Supports both treatment-level (coassembly) and sample-level (individual assembly) modes

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Set up temporary directory
TEMP_DIR=$(setup_temp_dir "semibin")

# ===== FUNCTION DEFINITIONS =====

# Run SemiBin2 binning
run_semibin() {
    local assembly_file="$1"
    local bam_dir="$2"
    local output_dir="$3"

    log "Running SemiBin2 binning module..."

    # Activate SemiBin environment
    activate_env semibin

    # Check if SemiBin2 is available
    if ! command -v SemiBin2 &> /dev/null; then
        log "ERROR: SemiBin2 not available in conda environment"
        conda deactivate
        return 1
    fi

    # Collect all BAM files
    local bam_files=("${bam_dir}"/*.sorted.bam)
    local bam_args=""
    for bam in "${bam_files[@]}"; do
        if [ -f "$bam" ]; then
            bam_args="$bam_args -b $bam"
        fi
    done

    if [ -z "$bam_args" ]; then
        log "ERROR: No BAM files found in $bam_dir"
        conda deactivate
        return 1
    fi

    log "Running SemiBin2 command:"
    log "SemiBin2 single_easy_bin -i $assembly_file $bam_args -o $output_dir -t ${SLURM_CPUS_PER_TASK:-50}"

    # Run SemiBin2 in single-sample mode with multiple BAM files
    SemiBin2 single_easy_bin \
        -i "$assembly_file" \
        $bam_args \
        -o "$output_dir" \
        -t ${SLURM_CPUS_PER_TASK:-50} \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_semibin.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    # Check results
    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/output_bins" ]; then
        local bin_count=$(ls -1 "${output_dir}/output_bins"/*.fa 2>/dev/null | wc -l)
        log "✓ SemiBin2 completed: $bin_count bins"
        return 0
    else
        log "✗ SemiBin2 failed with exit code: $exit_code"
        return 1
    fi
}

# Create summary
create_semibin_summary() {
    local semibin_dir="$1"
    local identifier="$2"  # treatment name or sample name
    local sample_count="$3"
    local summary_file="${semibin_dir}/semibin_summary.txt"

    cat > "$summary_file" << EOF
SemiBin2 Binning Summary
========================

Date: $(date)
Identifier: $identifier
Samples used: $sample_count
Tool: SemiBin2

Results:
EOF

    if [ -d "${semibin_dir}/output_bins" ]; then
        local total_bins=$(ls -1 "${semibin_dir}/output_bins"/*.fa 2>/dev/null | wc -l)
        echo "  Total bins: $total_bins" >> "$summary_file"
        echo "" >> "$summary_file"

        # Add bin size information
        echo "Bin Size Information:" >> "$summary_file"
        for bin in "${semibin_dir}/output_bins"/*.fa; do
            if [ -f "$bin" ]; then
                local bin_name=$(basename "$bin")
                local bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                local contig_count=$(grep -c "^>" "$bin")
                echo "  $bin_name: $bin_size bp, $contig_count contigs" >> "$summary_file"
            fi
        done
    else
        echo "  No bins directory found" >> "$summary_file"
    fi

    log "SemiBin2 summary created: $summary_file"
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

    log "====== Starting SemiBin2 Binning for Treatment: $TREATMENT ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"
    semibin_dir="${binning_dir}/semibin"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    mkdir -p "$semibin_dir"

    # Check if shared BAM directory exists
    if [ ! -d "$shared_bam_dir" ]; then
        log "ERROR: Shared BAM directory not found: $shared_bam_dir"
        log "Please run 03a_create_shared_bams.sh first"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Count BAM files
    bam_count=$(ls -1 "${shared_bam_dir}"/*.sorted.bam 2>/dev/null | wc -l)
    if [ $bam_count -eq 0 ]; then
        log "ERROR: No BAM files found in $shared_bam_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Found $bam_count BAM files in shared directory"

    # Check if already processed
    if [ -d "${semibin_dir}/output_bins" ] && [ "$(ls -A ${semibin_dir}/output_bins/*.fa 2>/dev/null)" ]; then
        log "SemiBin2 already completed for $TREATMENT, skipping..."
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
    log "Shared BAM directory: $shared_bam_dir"

    # Filter contigs by minimum length (1000 bp recommended for SemiBin)
    filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
    log "Filtering contigs to minimum length 1000 bp..."

    awk '/^>/ {if(seq) {if(length(seq) >= 1000) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1000) print header "\n" seq}' "$assembly_fasta" > "$filtered_assembly"

    filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
    log "After filtering: $filtered_contigs contigs >= 1000 bp"

    if [ $filtered_contigs -lt 10 ]; then
        log "ERROR: Too few contigs >= 1000 bp ($filtered_contigs) for meaningful binning"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Run SemiBin2 using shared BAM files
    if run_semibin "$filtered_assembly" "$shared_bam_dir" "$semibin_dir"; then
        create_semibin_summary "$semibin_dir" "$TREATMENT" "$bam_count"
        log "====== SemiBin2 binning completed for $TREATMENT ======"
    else
        log "ERROR: SemiBin2 binning failed for $TREATMENT"
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

    log "====== Starting SemiBin2 Binning for $SAMPLE_NAME ($TREATMENT) ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"
    semibin_dir="${binning_dir}/semibin"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    mkdir -p "$semibin_dir"

    # Check if shared BAM directory exists
    if [ ! -d "$shared_bam_dir" ]; then
        log "ERROR: Shared BAM directory not found: $shared_bam_dir"
        log "Please run 03a_create_shared_bams.sh first"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Count BAM files
    bam_count=$(ls -1 "${shared_bam_dir}"/*.sorted.bam 2>/dev/null | wc -l)
    if [ $bam_count -eq 0 ]; then
        log "ERROR: No BAM files found in $shared_bam_dir"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "Found $bam_count BAM files in shared directory"

    # Check if already processed
    if [ -d "${semibin_dir}/output_bins" ] && [ "$(ls -A ${semibin_dir}/output_bins/*.fa 2>/dev/null)" ]; then
        log "SemiBin2 already completed for $SAMPLE_NAME, skipping..."
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
    log "Shared BAM directory: $shared_bam_dir"

    # Filter contigs by minimum length (1000 bp recommended for SemiBin)
    filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
    log "Filtering contigs to minimum length 1000 bp..."

    awk '/^>/ {if(seq) {if(length(seq) >= 1000) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1000) print header "\n" seq}' "$assembly_fasta" > "$filtered_assembly"

    filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
    log "After filtering: $filtered_contigs contigs >= 1000 bp"

    if [ $filtered_contigs -lt 10 ]; then
        log "ERROR: Too few contigs >= 1000 bp ($filtered_contigs) for meaningful binning"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Run SemiBin2 using shared BAM files
    if run_semibin "$filtered_assembly" "$shared_bam_dir" "$semibin_dir"; then
        create_semibin_summary "$semibin_dir" "$SAMPLE_NAME" "$bam_count"
        log "====== SemiBin2 binning completed for $SAMPLE_NAME ======"
    else
        log "ERROR: SemiBin2 binning failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

cleanup_temp_dir "$TEMP_DIR"
log "====== SemiBin2 completed successfully ======"
