#!/bin/bash
#SBATCH --job-name=semibin
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=256G
#SBATCH --time=48:00:00
#SBATCH --account=ofinkel

# 03c_semibin.sh - Metagenomic binning using SemiBin2
# Supports both sample-level (individual assembly) and treatment-level (coassembly) modes

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# ===== FUNCTION DEFINITIONS (must be defined before use) =====

# Create BAM files for treatment-level (used in coassembly mode)
create_bam_files_treatment() {
    local assembly="$1"
    local read1="$2"
    local read2="$3"
    local output_dir="$4"
    local sample_name="$5"

    log "Creating BAM for sample $sample_name..."

    # Activate environment with bowtie2 and samtools
    activate_env metawrap-env

    # Use shared index if it exists, otherwise create it
    local index_base="${TEMP_DIR}/shared_index"
    if [ ! -f "${index_base}.1.bt2" ]; then
        log "Building shared bowtie2 index..."
        bowtie2-build \
            --threads ${SLURM_CPUS_PER_TASK:-50} \
            "$assembly" \
            "$index_base" \
            2>&1 | tee "${LOG_DIR}/${TREATMENT}_semibin_index.log"

        if [ ${PIPESTATUS[0]} -ne 0 ]; then
            log "ERROR: Failed to build bowtie2 index"
            conda deactivate
            return 1
        fi
    fi

    log "Mapping reads with bowtie2..."

    local sam_file="${TEMP_DIR}/${sample_name}.sam"
    local bam_file="${TEMP_DIR}/${sample_name}.bam"
    local sorted_bam="${output_dir}/${sample_name}.sorted.bam"

    bowtie2 \
        -x "$index_base" \
        -1 "$read1" \
        -2 "$read2" \
        -S "$sam_file" \
        -p ${SLURM_CPUS_PER_TASK:-50} \
        --sensitive \
        2>&1 | tee -a "${LOG_DIR}/${TREATMENT}_semibin_mapping.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Bowtie2 mapping failed for $sample_name"
        conda deactivate
        return 1
    fi

    log "Converting SAM to BAM and sorting..."

    samtools view -@ ${SLURM_CPUS_PER_TASK:-50} -bS "$sam_file" > "$bam_file"
    samtools sort -@ ${SLURM_CPUS_PER_TASK:-50} "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"

    # Clean up intermediate files
    rm -f "$sam_file" "$bam_file"

    conda deactivate

    log "BAM file created: $sorted_bam"
    return 0
}

# Create BAM files for sample-level (used in individual assembly mode)
create_bam_files_sample() {
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
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_semibin_index.log"

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
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_semibin_mapping.log"

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

# Run SemiBin2 binning
run_semibin() {
    local assembly_file="$1"
    local bam_dir="$2"
    local output_dir="$3"
    local identifier="$4"  # sample name or treatment name

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
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${identifier}_semibin.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    # Check results
    if [ $exit_code -eq 0 ] && [ -d "${output_dir}/output_bins" ]; then
        # Decompress any gzipped bins
        local gz_count=$(ls -1 "${output_dir}/output_bins"/*.fa.gz 2>/dev/null | wc -l)
        if [ $gz_count -gt 0 ]; then
            log "Decompressing $gz_count gzipped bin files..."
            gunzip "${output_dir}/output_bins"/*.fa.gz
            if [ $? -eq 0 ]; then
                log "Successfully decompressed all bins"
            else
                log "WARNING: Some bins failed to decompress"
            fi
        fi

        local bin_count=$(ls -1 "${output_dir}/output_bins"/*.fa 2>/dev/null | wc -l)
        log "SemiBin2 completed successfully with $bin_count bins"
        return 0
    else
        log "SemiBin2 failed with exit code: $exit_code"
        return 1
    fi
}

# Create summary for treatment-level
create_semibin_summary_treatment() {
    local treatment="$1"
    local semibin_dir="$2"
    local sample_count="$3"
    local summary_file="${semibin_dir}/semibin_summary.txt"

    cat > "$summary_file" << EOF
SemiBin2 Binning Summary for Treatment: $treatment
==================================================

Date: $(date)
Treatment: $treatment
Samples used: $sample_count
Mode: Treatment-level (coassembly)
Tool: SemiBin2

Results:
EOF

    local total_bins=0
    if [ -d "${semibin_dir}/output_bins" ]; then
        total_bins=$(ls -1 "${semibin_dir}/output_bins"/*.fa 2>/dev/null | wc -l)
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

# Create summary for sample-level
create_semibin_summary_sample() {
    local sample_name="$1"
    local treatment="$2"
    local semibin_dir="$3"
    local sample_count="${4:-1}"  # Optional: number of samples mapped (default 1)
    local summary_file="${semibin_dir}/semibin_summary.txt"

    local mode_desc="Sample-level (individual assembly)"
    if [ "$sample_count" -gt 1 ]; then
        mode_desc="Sample-level (individual assembly with cross-sample mapping)"
    fi

    cat > "$summary_file" << EOF
SemiBin2 Binning Summary for $sample_name
=========================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment
Mode: $mode_desc
Samples mapped: $sample_count
Tool: SemiBin2

Results:
EOF

    local total_bins=0
    if [ -d "${semibin_dir}/output_bins" ]; then
        total_bins=$(ls -1 "${semibin_dir}/output_bins"/*.fa 2>/dev/null | wc -l)
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

# Validation for treatment-level
validate_semibin_treatment() {
    local treatment="$1"
    local semibin_dir="${OUTPUT_DIR}/binning/${treatment}/semibin"

    # Check if bins were produced (both .fa and .fa.gz)
    if [ -d "${semibin_dir}/output_bins" ]; then
        local bin_count=$(ls -1 "${semibin_dir}/output_bins"/*.fa "${semibin_dir}/output_bins"/*.fa.gz 2>/dev/null | wc -l)
        if [ $bin_count -gt 0 ]; then
            return 0
        fi
    fi

    return 1
}

# Validation for sample-level
validate_semibin_sample() {
    local sample_name="$1"
    local treatment="$2"
    local semibin_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}/semibin"

    # Check if bins were produced (both .fa and .fa.gz)
    if [ -d "${semibin_dir}/output_bins" ]; then
        local bin_count=$(ls -1 "${semibin_dir}/output_bins"/*.fa "${semibin_dir}/output_bins"/*.fa.gz 2>/dev/null | wc -l)
        if [ $bin_count -gt 0 ]; then
            return 0
        fi
    fi

    return 1
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

    log "====== Starting SemiBin2 Binning for Treatment: $TREATMENT ======"

    # Check if stage already completed for this treatment
    if check_treatment_checkpoint "$TREATMENT" "semibin"; then
        log "SemiBin2 binning already completed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Treatment-level binning function
    stage_semibin_treatment() {
        local treatment="$1"

        log "Running SemiBin2 binning for treatment $treatment"

        local semibin_dir="${OUTPUT_DIR}/binning/${treatment}/semibin"
        local assembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"

        mkdir -p "$semibin_dir"

        # Check if already processed
        if [ -d "${semibin_dir}/output_bins" ] && [ "$(ls -A ${semibin_dir}/output_bins 2>/dev/null)" ]; then
            log "Treatment $treatment already binned with SemiBin2, skipping..."
            return 0
        fi

        # Clean up partial/failed SemiBin2 runs
        if [ -d "${semibin_dir}/data" ]; then
            log "Removing existing data directory from previous run..."
            rm -rf "${semibin_dir}/data"
        fi
        if [ -d "${semibin_dir}/model" ]; then
            log "Removing existing model directory from previous run..."
            rm -rf "${semibin_dir}/model"
        fi

        # Check for assembly
        local assembly_file=""
        for possible_file in \
            "${assembly_dir}/contigs.fasta" \
            "${assembly_dir}/scaffolds.fasta" \
            "${assembly_dir}/final_contigs.fasta" \
            "${assembly_dir}/assembly.fasta"; do

            if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                assembly_file="$possible_file"
                log "Found assembly file: $assembly_file"
                break
            fi
        done

        if [ -z "$assembly_file" ]; then
            log "ERROR: No assembly file found for treatment $treatment"
            log "  Checked locations in: $assembly_dir"
            return 1
        fi

        local num_contigs=$(grep -c "^>" "$assembly_file")
        log "Processing $num_contigs contigs for SemiBin2 binning"

        # Filter contigs by minimum length (1000 bp recommended)
        local filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
        log "Filtering contigs to minimum length 1000 bp..."

        awk '/^>/ {if(seq) {if(length(seq) >= 1000) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1000) print header "\n" seq}' "$assembly_file" > "$filtered_assembly"

        local filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
        log "After filtering: $filtered_contigs contigs >= 1000 bp"

        if [ $filtered_contigs -lt 10 ]; then
            log "ERROR: Too few contigs >= 1000 bp ($filtered_contigs) for meaningful binning"
            return 1
        fi

        # Setup BAM files directory - use shared location for all binning tools
        local bam_dir="${TEMP_DIR}/bam_files"
        local shared_bam_dir="${binning_dir}/shared_bam_files"  # Shared across COMEBin/SemiBin/etc
        local persistent_bam_dir="${semibin_dir}/bam_files"      # SemiBin-specific backup
        mkdir -p "$bam_dir"
        mkdir -p "$shared_bam_dir"
        mkdir -p "$persistent_bam_dir"

        # Get all samples in this treatment and create BAM files for each
        log "Finding all samples in treatment $treatment..."
        local total_samples=$(get_total_samples)
        local sample_count=0

        for i in $(seq 0 $((total_samples - 1))); do
            local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
            if [ -n "$sample_info" ]; then
                IFS='|' read -r sample_name sample_treatment _ _ <<< "$sample_info"

                if [ "$sample_treatment" = "$treatment" ]; then
                    log "Processing sample: $sample_name"

                    # Check if BAM already exists (shared location first, then tool-specific)
                    if [ -f "${shared_bam_dir}/${sample_name}.sorted.bam" ] && \
                       [ -f "${shared_bam_dir}/${sample_name}.sorted.bam.bai" ]; then
                        log "Found existing BAM for $sample_name in shared directory, reusing..."
                        ln -sf "${shared_bam_dir}/${sample_name}.sorted.bam" "${bam_dir}/${sample_name}.sorted.bam"
                        ln -sf "${shared_bam_dir}/${sample_name}.sorted.bam.bai" "${bam_dir}/${sample_name}.sorted.bam.bai"
                        ((sample_count++))
                        continue
                    elif [ -f "${persistent_bam_dir}/${sample_name}.sorted.bam" ] && \
                         [ -f "${persistent_bam_dir}/${sample_name}.sorted.bam.bai" ]; then
                        log "Found existing BAM for $sample_name in SemiBin directory, reusing..."
                        ln -sf "${persistent_bam_dir}/${sample_name}.sorted.bam" "${bam_dir}/${sample_name}.sorted.bam"
                        ln -sf "${persistent_bam_dir}/${sample_name}.sorted.bam.bai" "${bam_dir}/${sample_name}.sorted.bam.bai"
                        ((sample_count++))
                        continue
                    fi

                    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
                    local read1="${quality_dir}/filtered_1.fastq.gz"
                    local read2="${quality_dir}/filtered_2.fastq.gz"

                    if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
                        log "WARNING: Missing reads for $sample_name, skipping"
                        continue
                    fi

                    # Create BAM file for this sample
                    if create_bam_files_treatment "$filtered_assembly" "$read1" "$read2" "$bam_dir" "$sample_name"; then
                        # Copy to shared directory (primary) and tool-specific directory (backup)
                        cp "${bam_dir}/${sample_name}.sorted.bam" "${shared_bam_dir}/"
                        cp "${bam_dir}/${sample_name}.sorted.bam.bai" "${shared_bam_dir}/"
                        cp "${bam_dir}/${sample_name}.sorted.bam" "${persistent_bam_dir}/"
                        cp "${bam_dir}/${sample_name}.sorted.bam.bai" "${persistent_bam_dir}/"
                        ((sample_count++))
                        log "Successfully created BAM for $sample_name ($sample_count total)"
                    else
                        log "WARNING: Failed to create BAM for $sample_name"
                    fi
                fi
            fi
        done

        if [ $sample_count -eq 0 ]; then
            log "ERROR: No BAM files created for treatment $treatment"
            return 1
        fi

        log "Created BAM files for $sample_count samples in treatment $treatment"

        # Run SemiBin2
        if run_semibin "$filtered_assembly" "$bam_dir" "$semibin_dir" "$treatment"; then
            log "SemiBin2 binning completed successfully for treatment $treatment"
            create_semibin_summary_treatment "$treatment" "$semibin_dir" "$sample_count"
            return 0
        else
            log "ERROR: SemiBin2 binning failed for treatment $treatment"
            return 1
        fi
    }

    # Run treatment-level binning
    if stage_semibin_treatment "$TREATMENT"; then
        if validate_semibin_treatment "$TREATMENT"; then
            create_treatment_checkpoint "$TREATMENT" "semibin"
            log "====== SemiBin2 binning completed for treatment $TREATMENT ======"
        else
            log "ERROR: SemiBin2 validation failed for treatment $TREATMENT"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: SemiBin2 binning stage failed for treatment $TREATMENT"
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

    log "====== Starting SemiBin2 Binning for $SAMPLE_NAME ($TREATMENT) ======"

    # Check if stage already completed
    if check_sample_checkpoint "$SAMPLE_NAME" "semibin"; then
        log "SemiBin2 binning already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Sample-level binning function
    stage_semibin_sample() {
        local sample_name="$1"
        local treatment="$2"

        log "Running SemiBin2 binning for $sample_name ($treatment)"

        local semibin_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}/semibin"
        local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
        local metawrap_binning_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}"

        mkdir -p "$semibin_dir"

        # Check if already processed
        if [ -d "${semibin_dir}/output_bins" ] && [ "$(ls -A ${semibin_dir}/output_bins 2>/dev/null)" ]; then
            log "Sample $sample_name already binned with SemiBin2, skipping..."
            return 0
        fi

        # Clean up partial/failed SemiBin2 runs
        if [ -d "${semibin_dir}/data" ]; then
            log "Removing existing data directory from previous run..."
            rm -rf "${semibin_dir}/data"
        fi
        if [ -d "${semibin_dir}/model" ]; then
            log "Removing existing model directory from previous run..."
            rm -rf "${semibin_dir}/model"
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
            log "  Checked locations in: $assembly_dir"
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
        log "Processing $num_contigs contigs for SemiBin2 binning"

        # Filter contigs by minimum length (1000 bp recommended)
        local filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
        log "Filtering contigs to minimum length 1000 bp..."

        awk '/^>/ {if(seq) {if(length(seq) >= 1000) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1000) print header "\n" seq}' "$assembly_file" > "$filtered_assembly"

        local filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
        log "After filtering: $filtered_contigs contigs >= 1000 bp"

        if [ $filtered_contigs -lt 10 ]; then
            log "ERROR: Too few contigs >= 1000 bp ($filtered_contigs) for meaningful binning"
            return 1
        fi

        # Setup BAM files directory - use shared location for all binning tools
        local bam_dir="${TEMP_DIR}/bam_files"
        local shared_bam_dir="${binning_dir}/shared_bam_files"  # Shared across COMEBin/SemiBin/etc
        local persistent_bam_dir="${semibin_dir}/bam_files"      # SemiBin-specific backup
        mkdir -p "$bam_dir"
        mkdir -p "$shared_bam_dir"
        mkdir -p "$persistent_bam_dir"

        # Check if per-sample cross-mapping is enabled
        if [ "${PER_SAMPLE_CROSS_MAPPING}" = "true" ]; then
            log "Per-sample cross-mapping enabled: mapping all samples in treatment to this assembly"

            # Get all samples in this treatment and create BAM files for each
            log "Finding all samples in treatment $treatment..."
            local total_samples=$(get_total_samples)
            local sample_count=0

            for i in $(seq 0 $((total_samples - 1))); do
                local other_sample_info=$(get_sample_info_by_index $i 2>/dev/null)
                if [ -n "$other_sample_info" ]; then
                    IFS='|' read -r other_sample_name other_treatment _ _ <<< "$other_sample_info"

                    if [ "$other_treatment" = "$treatment" ]; then
                        log "Processing sample: $other_sample_name"

                        # Check if BAM already exists (shared location first, then tool-specific)
                        if [ -f "${shared_bam_dir}/${other_sample_name}.sorted.bam" ] && \
                           [ -f "${shared_bam_dir}/${other_sample_name}.sorted.bam.bai" ]; then
                            log "Found existing BAM for $other_sample_name in shared directory, reusing..."
                            ln -sf "${shared_bam_dir}/${other_sample_name}.sorted.bam" "${bam_dir}/${other_sample_name}.sorted.bam"
                            ln -sf "${shared_bam_dir}/${other_sample_name}.sorted.bam.bai" "${bam_dir}/${other_sample_name}.sorted.bam.bai"
                            ((sample_count++))
                            continue
                        elif [ -f "${persistent_bam_dir}/${other_sample_name}.sorted.bam" ] && \
                             [ -f "${persistent_bam_dir}/${other_sample_name}.sorted.bam.bai" ]; then
                            log "Found existing BAM for $other_sample_name in SemiBin directory, reusing..."
                            ln -sf "${persistent_bam_dir}/${other_sample_name}.sorted.bam" "${bam_dir}/${other_sample_name}.sorted.bam"
                            ln -sf "${persistent_bam_dir}/${other_sample_name}.sorted.bam.bai" "${bam_dir}/${other_sample_name}.sorted.bam.bai"
                            ((sample_count++))
                            continue
                        fi

                        local other_quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${other_sample_name}"
                        local other_read1="${other_quality_dir}/filtered_1.fastq.gz"
                        local other_read2="${other_quality_dir}/filtered_2.fastq.gz"

                        if [ ! -f "$other_read1" ] || [ ! -f "$other_read2" ]; then
                            log "WARNING: Missing reads for $other_sample_name, skipping"
                            continue
                        fi

                        # Create BAM file for this sample mapping to current sample's assembly
                        if create_bam_files_sample "$filtered_assembly" "$other_read1" "$other_read2" "$bam_dir" "$other_sample_name"; then
                            # Copy to shared directory (primary) and tool-specific directory (backup)
                            cp "${bam_dir}/${other_sample_name}.sorted.bam" "${shared_bam_dir}/"
                            cp "${bam_dir}/${other_sample_name}.sorted.bam.bai" "${shared_bam_dir}/"
                            cp "${bam_dir}/${other_sample_name}.sorted.bam" "${persistent_bam_dir}/"
                            cp "${bam_dir}/${other_sample_name}.sorted.bam.bai" "${persistent_bam_dir}/"
                            ((sample_count++))
                            log "Successfully created BAM for $other_sample_name ($sample_count total)"
                        else
                            log "WARNING: Failed to create BAM for $other_sample_name"
                        fi
                    fi
                fi
            done

            if [ $sample_count -eq 0 ]; then
                log "ERROR: No BAM files created for treatment $treatment"
                return 1
            fi

            log "Created BAM files for $sample_count samples in treatment $treatment"

        else
            # Original behavior: single-sample mapping
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
                sample_count=1
            else
                # Create BAM files using bowtie2
                log "Creating BAM files for SemiBin2..."

                if ! create_bam_files_sample "$filtered_assembly" "$read1" "$read2" "$bam_dir" "$sample_name"; then
                    log "ERROR: Failed to create BAM files"
                    return 1
                fi
                sample_count=1
            fi
        fi

        # Run SemiBin2
        if run_semibin "$filtered_assembly" "$bam_dir" "$semibin_dir" "$sample_name"; then
            log "SemiBin2 binning completed successfully"
            # Pass sample count if cross-mapping was used
            if [ "${PER_SAMPLE_CROSS_MAPPING}" = "true" ]; then
                create_semibin_summary_sample "$sample_name" "$treatment" "$semibin_dir" "$sample_count"
            else
                create_semibin_summary_sample "$sample_name" "$treatment" "$semibin_dir"
            fi
            return 0
        else
            log "ERROR: SemiBin2 binning failed for $sample_name"
            return 1
        fi
    }

    # Run sample-level binning
    if stage_semibin_sample "$SAMPLE_NAME" "$TREATMENT"; then
        if validate_semibin_sample "$SAMPLE_NAME" "$TREATMENT"; then
            create_sample_checkpoint "$SAMPLE_NAME" "semibin"
            log "====== SemiBin2 binning completed for $SAMPLE_NAME ======"
        else
            log "ERROR: SemiBin2 validation failed for $SAMPLE_NAME"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: SemiBin2 binning stage failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
