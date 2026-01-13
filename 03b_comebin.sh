#!/bin/bash
#SBATCH --job-name=comebin
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --account=ofinkel

# 03b_comebin.sh - Metagenomic binning using COMEBin

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get sample info from array task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
SAMPLE_INFO=$(get_sample_info_by_index "$TASK_ID")
if [ -z "$SAMPLE_INFO" ]; then
    echo "ERROR: No sample found for array index $TASK_ID"
    exit 1
fi

# Parse sample information
IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting COMEBin Binning for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "comebin"; then
    log "COMEBin binning already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_comebin() {
    local sample_name="$1"
    local treatment="$2"

    log "Running COMEBin binning for $sample_name ($treatment)"

    local comebin_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}/comebin"
    local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    local metawrap_binning_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}"

    mkdir -p "$comebin_dir"

    # Check if already processed
    if [ -d "${comebin_dir}/comebin_res_bins" ] && [ "$(ls -A ${comebin_dir}/comebin_res_bins 2>/dev/null)" ]; then
        log "Sample $sample_name already binned with COMEBin, skipping..."
        return 0
    fi

    # Check for input files - assembly
    local assembly_file=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/assembly.fasta" \
        "${assembly_dir}/${sample_name}_contigs.fasta" \
        "${assembly_dir}/${sample_name}_scaffolds.fasta"; do

        if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
            assembly_file="$possible_file"
            log "Found assembly file: $assembly_file"
            break
        fi
    done

    if [ -z "$assembly_file" ]; then
        log "ERROR: No assembly file found for $sample_name"
        return 1
    fi

    # Check for quality-filtered reads
    local read1="${quality_dir}/filtered_1.fastq.gz"
    local read2="${quality_dir}/filtered_2.fastq.gz"

    if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
        log "ERROR: Missing quality-filtered reads for $sample_name"
        log "  Expected: $read1 and $read2"
        return 1
    fi

    local num_contigs=$(grep -c "^>" "$assembly_file")
    log "Processing $num_contigs contigs for COMEBin binning"

    # Filter contigs by minimum length (1000 bp recommended for COMEBin)
    local filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
    log "Filtering contigs to minimum length 1000 bp..."

    awk '/^>/ {if(seq) {if(length(seq) >= 1000) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1000) print header "\n" seq}' "$assembly_file" > "$filtered_assembly"

    local filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
    log "After filtering: $filtered_contigs contigs >= 1000 bp"

    if [ $filtered_contigs -lt 10 ]; then
        log "ERROR: Too few contigs >= 1000 bp ($filtered_contigs) for meaningful binning"
        return 1
    fi

    # Setup BAM files directory
    local bam_dir="${TEMP_DIR}/bam_files"
    mkdir -p "$bam_dir"

    # Check if BAM files exist from MetaWRAP
    local metawrap_bam=""
    if [ -d "${metawrap_binning_dir}/work_files" ]; then
        # Look for sorted BAM file from MetaWRAP
        for bam_file in "${metawrap_binning_dir}/work_files"/*.bam; do
            if [ -f "$bam_file" ] && [[ "$bam_file" == *"sorted"* ]]; then
                metawrap_bam="$bam_file"
                log "Found MetaWRAP BAM file: $metawrap_bam"
                break
            fi
        done
    fi

    if [ -n "$metawrap_bam" ] && [ -f "$metawrap_bam" ]; then
        # Use existing BAM from MetaWRAP
        log "Using existing BAM file from MetaWRAP"
        ln -s "$metawrap_bam" "${bam_dir}/${sample_name}.sorted.bam"
    else
        # Create BAM files using bowtie2
        log "Creating BAM files for COMEBin..."

        if ! create_bam_files "$filtered_assembly" "$read1" "$read2" "$bam_dir" "$sample_name"; then
            log "ERROR: Failed to create BAM files"
            return 1
        fi
    fi

    # Run COMEBin
    if run_comebin "$filtered_assembly" "$bam_dir" "$comebin_dir" "$sample_name"; then
        log "COMEBin binning completed successfully"
        create_comebin_summary "$sample_name" "$treatment" "$comebin_dir"
        return 0
    else
        log "ERROR: COMEBin binning failed for $sample_name"
        return 1
    fi
}

# Create BAM files using bowtie2
create_bam_files() {
    local assembly="$1"
    local read1="$2"
    local read2="$3"
    local output_dir="$4"
    local sample_name="$5"

    log "Building bowtie2 index..."

    # Activate environment with bowtie2 and samtools
    activate_env metawrap-env

    local index_base="${output_dir}/${sample_name}_index"

    bowtie2-build \
        --threads ${SLURM_CPUS_PER_TASK:-50} \
        "$assembly" \
        "$index_base" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_comebin_index.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Failed to build bowtie2 index"
        conda deactivate
        return 1
    fi

    log "Mapping reads with bowtie2..."

    local sam_file="${output_dir}/${sample_name}.sam"
    local bam_file="${output_dir}/${sample_name}.bam"
    local sorted_bam="${output_dir}/${sample_name}.sorted.bam"

    bowtie2 \
        -x "$index_base" \
        -1 "$read1" \
        -2 "$read2" \
        -S "$sam_file" \
        -p ${SLURM_CPUS_PER_TASK:-50} \
        --sensitive \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_comebin_mapping.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Bowtie2 mapping failed"
        conda deactivate
        return 1
    fi

    log "Converting SAM to BAM and sorting..."

    samtools view -@ ${SLURM_CPUS_PER_TASK:-50} -bS "$sam_file" > "$bam_file"
    samtools sort -@ ${SLURM_CPUS_PER_TASK:-50} "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"

    # Clean up intermediate files
    rm -f "$sam_file" "$bam_file"
    rm -f "${index_base}".*

    conda deactivate

    log "BAM file created: $sorted_bam"
    return 0
}

# Run COMEBin binning
run_comebin() {
    local assembly_file="$1"
    local bam_dir="$2"
    local output_dir="$3"
    local sample_name="$4"

    log "Running COMEBin binning module..."

    # Activate COMEBin environment
    activate_env comebin

    # Check if COMEBin is available
    if [ ! -f "${CONDA_BASE}/envs/comebin/bin/run_comebin.sh" ]; then
        log "ERROR: COMEBin not available in conda environment"
        log "Expected: ${CONDA_BASE}/envs/comebin/bin/run_comebin.sh"
        conda deactivate
        return 1
    fi

    log "Running COMEBin command:"
    log "bash run_comebin.sh -a $assembly_file -o $output_dir -p $bam_dir -t ${SLURM_CPUS_PER_TASK:-50}"

    # Run COMEBin
    bash "${CONDA_BASE}/envs/comebin/bin/run_comebin.sh" \
        -a "$assembly_file" \
        -o "$output_dir" \
        -p "$bam_dir" \
        -t ${SLURM_CPUS_PER_TASK:-50} \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_comebin.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    # Check results
    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/comebin_res_bins" ]; then
        local bin_count=$(ls -1 "${output_dir}/comebin_res_bins"/*.fa 2>/dev/null | wc -l)
        log "COMEBin completed successfully with $bin_count bins"
        return 0
    else
        log "COMEBin failed with exit code: $exit_code"
        return 1
    fi
}

# Create COMEBin summary
create_comebin_summary() {
    local sample_name="$1"
    local treatment="$2"
    local comebin_dir="$3"
    local summary_file="${comebin_dir}/comebin_summary.txt"

    cat > "$summary_file" << EOF
COMEBin Binning Summary for $sample_name
=========================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment
Tool: COMEBin

Results:
EOF

    local total_bins=0
    if [ -d "${comebin_dir}/comebin_res_bins" ]; then
        total_bins=$(ls -1 "${comebin_dir}/comebin_res_bins"/*.fa 2>/dev/null | wc -l)
        echo "  Total bins: $total_bins" >> "$summary_file"
        echo "" >> "$summary_file"

        # Add bin size information
        echo "Bin Size Information:" >> "$summary_file"
        for bin in "${comebin_dir}/comebin_res_bins"/*.fa; do
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

    log "COMEBin summary created: $summary_file"
}

# Validation function
validate_comebin() {
    local sample_name="$1"
    local treatment="$2"
    local comebin_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}/comebin"

    # Check if bins were produced
    if [ -d "${comebin_dir}/comebin_res_bins" ] && [ "$(ls -A ${comebin_dir}/comebin_res_bins/*.fa 2>/dev/null)" ]; then
        return 0
    fi

    return 1
}

# Run the COMEBin binning stage
if stage_comebin "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_comebin "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "comebin"
        log "====== COMEBin binning completed for $SAMPLE_NAME ======"
    else
        log "ERROR: COMEBin validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: COMEBin binning stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
