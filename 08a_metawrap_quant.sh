#!/bin/bash
#SBATCH --job-name=metawrap_quant
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# 08a_metawrap_quant.sh - MetaWRAP quant_bins for refined bins (before reassembly)

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

    log "====== Starting MetaWRAP Quant_bins for treatment $TREATMENT (treatment-level mode) ======"

    # Check if stage already completed for treatment
    if check_treatment_checkpoint "$TREATMENT" "metawrap_quant"; then
        log "MetaWRAP quant_bins already completed for treatment $TREATMENT"
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

    log "====== Starting MetaWRAP Quant_bins for $SAMPLE_NAME ($TREATMENT) ======"

    # Check if stage already completed for sample
    if check_sample_checkpoint "$SAMPLE_NAME" "metawrap_quant"; then
        log "MetaWRAP quant_bins already completed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi
fi

# Function to get refined bins directory (from stage 04 - DAS Tool output)
get_refined_bins_dir() {
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        local entity_desc="sample $sample_name ($treatment)"
    else
        local entity_desc="treatment $treatment"
    fi

    log "Locating refined bins for $entity_desc..."

    # Check treatment-level bins first (for coassembly/treatment-level binning)
    local treatment_bins="${OUTPUT_DIR}/bin_refinement/${treatment}/dastool_DASTool_bins"
    if [ -d "$treatment_bins" ] && [ "$(ls -A "$treatment_bins"/*.fa 2>/dev/null)" ]; then
        local bin_count=$(ls -1 "$treatment_bins"/*.fa 2>/dev/null | wc -l)
        log "  Found $bin_count treatment-level refined bins at: $treatment_bins"
        echo "$treatment_bins"
        return 0
    fi

    # Check sample-level bins (for individual sample binning)
    if [ -n "$sample_name" ]; then
        local sample_bins="${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}/dastool_DASTool_bins"
        if [ -d "$sample_bins" ] && [ "$(ls -A "$sample_bins"/*.fa 2>/dev/null)" ]; then
            local bin_count=$(ls -1 "$sample_bins"/*.fa 2>/dev/null | wc -l)
            log "  Found $bin_count sample-level refined bins at: $sample_bins"
            echo "$sample_bins"
            return 0
        fi

        log "  ERROR: No refined bins found at either location:"
        log "    Treatment-level: $treatment_bins"
        log "    Sample-level: $sample_bins"
    else
        log "  ERROR: No refined bins found at: $treatment_bins"
    fi

    return 1
}

# Function to get reads for quantification
get_reads_for_quant() {
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        # Sample-level: use individual sample reads
        log "Using sample-level reads for quantification..."

        # Check validated reads first, fall back to quality filtered
        local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"

        if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
            local r1="${validated_dir}/validated_1.fastq.gz"
            local r2="${validated_dir}/validated_2.fastq.gz"
            log "  Using validated reads: $validated_dir"
        elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
            local r1="${quality_dir}/filtered_1.fastq.gz"
            local r2="${quality_dir}/filtered_2.fastq.gz"
            log "  Using quality-filtered reads: $quality_dir"
        else
            log "  ERROR: No reads found for sample $sample_name"
            return 1
        fi

        echo "${r1}|${r2}"
        return 0
    else
        # Treatment-level: use merged treatment reads (from coassembly)
        log "Using treatment-level merged reads for quantification..."

        local merged_reads_dir="${OUTPUT_DIR}/coassembly/${treatment}/merged_reads"
        local r1="${merged_reads_dir}/merged_R1.fastq.gz"
        local r2="${merged_reads_dir}/merged_R2.fastq.gz"

        if [ -f "$r1" ] && [ -f "$r2" ]; then
            log "  Using merged reads: $merged_reads_dir"
            echo "${r1}|${r2}"
            return 0
        else
            log "  ERROR: Merged reads not found at: $merged_reads_dir"
            log "  This suggests coassembly was not run or merged reads were not saved"
            return 1
        fi
    fi
}

# Main processing function
stage_metawrap_quant() {
    local treatment="$1"
    local sample_name="${2:-}"  # Optional for treatment-level mode

    if [ -n "$sample_name" ]; then
        # Sample-level mode
        log "Running MetaWRAP quant_bins for $sample_name ($treatment)"
        local output_dir="${OUTPUT_DIR}/metawrap_quant/${treatment}/${sample_name}"
        local entity_desc="sample $sample_name"
    else
        # Treatment-level mode
        log "Running MetaWRAP quant_bins for treatment $treatment"
        local output_dir="${OUTPUT_DIR}/metawrap_quant/${treatment}"
        local entity_desc="treatment $treatment"
    fi

    mkdir -p "$output_dir"

    # Check if already processed
    if [ -f "${output_dir}/bin_abundance_table.tab" ]; then
        log "MetaWRAP quant_bins already completed for $entity_desc, skipping..."
        return 0
    fi

    # Get refined bins directory
    local bins_dir=$(get_refined_bins_dir "$treatment" "$sample_name")
    if [ $? -ne 0 ] || [ -z "$bins_dir" ]; then
        log "ERROR: No refined bins found for $entity_desc"
        return 1
    fi

    # Get reads for quantification
    local reads_info=$(get_reads_for_quant "$treatment" "$sample_name")
    if [ $? -ne 0 ] || [ -z "$reads_info" ]; then
        log "ERROR: Could not locate reads for quantification"
        return 1
    fi

    IFS='|' read -r r1 r2 <<< "$reads_info"

    # Validate read files
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        log "ERROR: Read files not found:"
        log "  R1: $r1"
        log "  R2: $r2"
        return 1
    fi

    log "Read files for quantification:"
    log "  R1: $r1"
    log "  R2: $r2"

    # Count bins
    local bin_count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
    log "Quantifying $bin_count refined bins..."

    if [ $bin_count -eq 0 ]; then
        log "ERROR: No bins found in $bins_dir"
        return 1
    fi

    # Activate MetaWRAP environment
    activate_env metawrap

    # Check if MetaWRAP is available
    if ! command -v metawrap &> /dev/null; then
        log "ERROR: MetaWRAP not available in environment"
        conda deactivate
        return 1
    fi

    # Run MetaWRAP quant_bins
    log "Running MetaWRAP quant_bins..."
    log "Command: metawrap quant_bins -b $bins_dir -o $output_dir -a ${ASSEMBLY_FILE} -t $SLURM_CPUS_PER_TASK $r1 $r2"

    # For quant_bins, we need the assembly file
    # Determine which assembly to use
    local assembly_file=""
    if [ -n "$sample_name" ]; then
        # Sample-level: use individual sample assembly
        assembly_file="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}/contigs.fasta"
    else
        # Treatment-level: use coassembly
        assembly_file="${OUTPUT_DIR}/coassembly/${treatment}/contigs.fasta"
    fi

    if [ ! -f "$assembly_file" ]; then
        log "ERROR: Assembly file not found: $assembly_file"
        conda deactivate
        return 1
    fi

    log "Using assembly: $assembly_file"

    # Run MetaWRAP quant_bins
    metawrap quant_bins \
        -b "$bins_dir" \
        -o "$output_dir" \
        -a "$assembly_file" \
        -t $SLURM_CPUS_PER_TASK \
        "$r1" "$r2" \
        2>&1 | tee "${LOG_DIR}/${treatment}/$([ -n "$sample_name" ] && echo "${sample_name}_" || echo "")metawrap_quant.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    # Check if output was generated
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/bin_abundance_table.tab" ]; then
        log "MetaWRAP quant_bins completed successfully for $entity_desc"

        # Log some statistics
        if [ -f "${output_dir}/bin_abundance_table.tab" ]; then
            local bins_quantified=$(tail -n +2 "${output_dir}/bin_abundance_table.tab" | wc -l)
            log "Successfully quantified $bins_quantified bins"
        fi

        return 0
    else
        log "ERROR: MetaWRAP quant_bins failed for $entity_desc (exit code: $exit_code)"
        return 1
    fi
}

# Validation function
validate_metawrap_quant() {
    local treatment="$1"
    local sample_name="${2:-}"

    if [ -n "$sample_name" ]; then
        local output_dir="${OUTPUT_DIR}/metawrap_quant/${treatment}/${sample_name}"
    else
        local output_dir="${OUTPUT_DIR}/metawrap_quant/${treatment}"
    fi

    # Check if abundance table exists and is not empty
    if [ -f "${output_dir}/bin_abundance_table.tab" ] && [ -s "${output_dir}/bin_abundance_table.tab" ]; then
        local bins_quantified=$(tail -n +2 "${output_dir}/bin_abundance_table.tab" | wc -l)
        if [ $bins_quantified -gt 0 ]; then
            log "Validation successful: Quantified $bins_quantified bins"
            return 0
        fi
    fi

    log "Validation failed: No valid quant_bins results found"
    return 1
}

# Run the MetaWRAP quant_bins stage based on processing mode
if [ "$PROCESSING_MODE" = "treatment-level" ]; then
    # Treatment-level mode: run once for all bins in treatment
    if stage_metawrap_quant "$TREATMENT"; then
        # Validate results
        if validate_metawrap_quant "$TREATMENT"; then
            create_treatment_checkpoint "$TREATMENT" "metawrap_quant"
            log "====== MetaWRAP quant_bins completed successfully for treatment $TREATMENT ======"
        else
            log "ERROR: MetaWRAP quant_bins validation failed for treatment $TREATMENT"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: MetaWRAP quant_bins stage failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

else
    # Sample-level mode: run for individual sample
    if stage_metawrap_quant "$TREATMENT" "$SAMPLE_NAME"; then
        # Validate results
        if validate_metawrap_quant "$TREATMENT" "$SAMPLE_NAME"; then
            create_sample_checkpoint "$SAMPLE_NAME" "metawrap_quant"
            log "====== MetaWRAP quant_bins completed successfully for $SAMPLE_NAME ======"
        else
            log "ERROR: MetaWRAP quant_bins validation failed for $SAMPLE_NAME"
            cleanup_temp_dir "$TEMP_DIR"
            exit 1
        fi
    else
        log "ERROR: MetaWRAP quant_bins stage failed for $SAMPLE_NAME"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
