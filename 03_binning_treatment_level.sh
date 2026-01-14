#!/bin/bash
#SBATCH --job-name=binning_treatment
#SBATCH --array=0-9%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# 03_binning_treatment_level.sh - Treatment-level metagenomic binning using MetaWRAP
# Uses coassembly contigs and merged reads for treatment-level binning

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Verify treatments file exists
if [ ! -f "$TREATMENTS_FILE" ]; then
    echo "ERROR: Treatments file not found at: $TREATMENTS_FILE"
    exit 1
fi

# Get treatment from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
mapfile -t TREATMENTS_ARRAY < "$TREATMENTS_FILE"

if [ $TASK_ID -ge ${#TREATMENTS_ARRAY[@]} ]; then
    echo "No treatment found for array index $TASK_ID"
    exit 0
fi

TREATMENT="${TREATMENTS_ARRAY[$TASK_ID]}"
export TREATMENT

# Use treatment name for checkpoints
SAMPLE_NAME="${TREATMENT}"
export SAMPLE_NAME

# Initialize
init_conda
create_treatment_dirs "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Treatment-Level Binning for $TREATMENT ======"

# Check if stage already completed (using treatment-level checkpoint)
if check_sample_checkpoint "$TREATMENT" "binning"; then
    log "Binning already completed for treatment $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_binning_treatment() {
    local treatment="$1"

    log "Running treatment-level binning for $treatment"

    local binning_dir="${OUTPUT_DIR}/binning/${treatment}"
    local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
    local merged_reads_dir="${coassembly_dir}/merged_reads"

    mkdir -p "$binning_dir"

    # Check if already processed
    if ls "${binning_dir}"/*_bins/*.fa &>/dev/null 2>&1; then
        log "Treatment $treatment already binned, skipping..."
        return 0
    fi

    # Check for coassembly contigs
    local assembly_file="${coassembly_dir}/contigs.fasta"

    if [ ! -f "$assembly_file" ] || [ ! -s "$assembly_file" ]; then
        log "ERROR: No coassembly found for treatment $treatment"
        log "  Expected: $assembly_file"
        return 1
    fi

    log "Found coassembly file: $assembly_file"

    # Check for merged reads from coassembly
    local read1=""
    local read2=""

    # Try different naming patterns for merged reads
    for r1_pattern in "merged_R1.fastq.gz" "merged_1.fastq.gz" "all_R1.fastq.gz"; do
        local r1_path="${merged_reads_dir}/${r1_pattern}"
        if [ -f "$r1_path" ] && [ -s "$r1_path" ]; then
            read1="$r1_path"
            break
        fi
    done

    for r2_pattern in "merged_R2.fastq.gz" "merged_2.fastq.gz" "all_R2.fastq.gz"; do
        local r2_path="${merged_reads_dir}/${r2_pattern}"
        if [ -f "$r2_path" ] && [ -s "$r2_path" ]; then
            read2="$r2_path"
            break
        fi
    done

    if [ -z "$read1" ] || [ -z "$read2" ]; then
        log "ERROR: Missing merged reads for treatment $treatment"
        log "  Expected location: $merged_reads_dir"
        log "  Checked patterns: merged_R1.fastq.gz, merged_1.fastq.gz, all_R1.fastq.gz"
        return 1
    fi

    log "Found merged reads:"
    log "  R1: $read1"
    log "  R2: $read2"

    local num_contigs=$(grep -c "^>" "$assembly_file")
    log "Binning $num_contigs contigs for treatment $treatment"

    # Check minimum contig requirement
    if [ $num_contigs -lt 100 ]; then
        log "WARNING: Very few contigs ($num_contigs) for binning. Consider increasing assembly quality."
    fi

    # Filter contigs by minimum length (1500 bp) for binning
    local filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
    log "Filtering contigs to minimum length 1500 bp..."

    awk '/^>/ {if(seq) {if(length(seq) >= 1500) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1500) print header "\n" seq}' "$assembly_file" > "$filtered_assembly"

    local filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
    log "After filtering: $filtered_contigs contigs >= 1500 bp"

    if [ $filtered_contigs -lt 10 ]; then
        log "ERROR: Too few contigs >= 1500 bp ($filtered_contigs) for meaningful binning"
        return 1
    fi

    # Prepare reads for MetaWRAP (proper naming convention)
    local metawrap_reads_dir="${TEMP_DIR}/metawrap_reads"
    mkdir -p "$metawrap_reads_dir"

    # MetaWRAP expects files ending in _1.fastq and _2.fastq (uncompressed)
    local metawrap_r1="${metawrap_reads_dir}/${treatment}_1.fastq"
    local metawrap_r2="${metawrap_reads_dir}/${treatment}_2.fastq"

    log "Creating uncompressed reads for MetaWRAP..."
    gunzip -c "$read1" > "$metawrap_r1"
    gunzip -c "$read2" > "$metawrap_r2"

    log "Created MetaWRAP-compatible read files:"
    log "  R1: $metawrap_r1"
    log "  R2: $metawrap_r2"

    # Run MetaWRAP binning
    if run_metawrap_binning "$filtered_assembly" "$metawrap_r1" "$metawrap_r2" "$binning_dir"; then
        log "MetaWRAP binning completed successfully"
        create_binning_summary "$treatment" "$binning_dir"
        return 0
    else
        log "ERROR: MetaWRAP binning failed for treatment $treatment"
        return 1
    fi
}

# Run MetaWRAP binning module
run_metawrap_binning() {
    local assembly_file="$1"
    local read1="$2"
    local read2="$3"
    local output_dir="$4"

    log "Running MetaWRAP binning module..."

    # Activate MetaWRAP environment
    activate_env metawrap-env

    # Check if MetaWRAP is available
    if ! command -v metawrap &> /dev/null; then
        log "ERROR: MetaWRAP not available"
        conda deactivate
        return 1
    fi

    # Create temporary output directory for MetaWRAP
    local metawrap_output="${TEMP_DIR}/metawrap_binning"
    mkdir -p "$metawrap_output"

    # Run MetaWRAP binning
    log "Running MetaWRAP binning command:"
    log "metawrap binning -o $metawrap_output -t $SLURM_CPUS_PER_TASK -a $assembly_file --metabat2 --maxbin2 --concoct $read1 $read2"

    metawrap binning \
        -o "$metawrap_output" \
        -t $SLURM_CPUS_PER_TASK \
        -m $((${SLURM_MEM_PER_NODE:-128000} / 1000)) \
        -l 1500 \
        -a "$assembly_file" \
        --metabat2 --maxbin2 --concoct \
        "$read1" "$read2" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${TREATMENT}_metawrap_binning.log"

    local exit_code=${PIPESTATUS[0]}

    # Check results and copy to final output directory
    if [ $exit_code -eq 0 ] && [ -d "$metawrap_output" ]; then
        log "MetaWRAP binning completed successfully"

        # Copy results to final output directory
        copy_metawrap_results "$metawrap_output" "$output_dir"

        conda deactivate
        return 0
    else
        log "MetaWRAP binning failed with exit code: $exit_code"
        conda deactivate
        return 1
    fi
}

# Copy MetaWRAP results to final output directory
copy_metawrap_results() {
    local metawrap_output="$1"
    local final_output="$2"

    log "Copying MetaWRAP results to final output directory..."

    # Copy binner outputs
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${metawrap_output}/${binner}" ]; then
            cp -r "${metawrap_output}/${binner}" "$final_output/"
            log "Copied ${binner} results"
        else
            log "No results found for ${binner}"
        fi
    done

    # Copy other important files
    for file in binning_results.txt bin_abundance_table.tab; do
        if [ -f "${metawrap_output}/${file}" ]; then
            cp "${metawrap_output}/${file}" "$final_output/"
            log "Copied ${file}"
        fi
    done

    # Copy log files
    if [ -d "${metawrap_output}/work_files" ]; then
        cp -r "${metawrap_output}/work_files" "$final_output/"
    fi
}

# Create binning summary
create_binning_summary() {
    local treatment="$1"
    local binning_dir="$2"
    local summary_file="${binning_dir}/binning_summary.txt"

    cat > "$summary_file" << EOF
Treatment-Level Binning Summary for $treatment
===============================================

Date: $(date)
Treatment: $treatment
Mode: Treatment-level (coassembly)
Tool: MetaWRAP binning module

Results:
EOF

    local total_bins=0
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${binning_dir}/${binner}" ]; then
            local bin_count=$(ls -1 "${binning_dir}/${binner}"/*.fa 2>/dev/null | wc -l)
            echo "  $(echo $binner | sed 's/_bins$//'): $bin_count bins" >> "$summary_file"
            total_bins=$((total_bins + bin_count))
        else
            echo "  $(echo $binner | sed 's/_bins$//'): No results" >> "$summary_file"
        fi
    done

    echo "  Total bins: $total_bins" >> "$summary_file"
    echo "" >> "$summary_file"

    # Add bin size information
    echo "Bin Size Information:" >> "$summary_file"
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${binning_dir}/${binner}" ]; then
            echo "  $(echo $binner | sed 's/_bins$//) bins:" >> "$summary_file"
            for bin in "${binning_dir}/${binner}"/*.fa; do
                if [ -f "$bin" ]; then
                    local bin_name=$(basename "$bin")
                    local bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    local contig_count=$(grep -c "^>" "$bin")
                    echo "    $bin_name: $bin_size bp, $contig_count contigs" >> "$summary_file"
                fi
            done
        fi
    done

    log "Binning summary created: $summary_file"
}

# Validation function
validate_binning() {
    local treatment="$1"
    local binning_dir="${OUTPUT_DIR}/binning/${treatment}"

    # Check if at least one binner produced bins
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${binning_dir}/${binner}" ] && [ "$(ls -A "${binning_dir}/${binner}/"*.fa 2>/dev/null)" ]; then
            return 0
        fi
    done

    return 1
}

# Run the binning stage
if stage_binning_treatment "$TREATMENT"; then
    # Validate results
    if validate_binning "$TREATMENT"; then
        create_sample_checkpoint "$TREATMENT" "binning"
        log "====== Treatment-level binning completed for $TREATMENT ======"
    else
        log "ERROR: Binning validation failed for treatment $TREATMENT"
        exit 1
    fi
else
    log "ERROR: Binning stage failed for treatment $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
