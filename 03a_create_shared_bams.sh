#!/usr/bin/env bash
#SBATCH --job-name=create_bams
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --account=ofinkel

# 03a_create_shared_bams.sh - Create shared BAM files for all binners
# This script runs BEFORE binning and creates BAM files that will be reused by:
# - MetaWRAP binning (metabat2, maxbin2, concoct)
# - COMEBin
# - SemiBin

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Set up temporary directory
TEMP_DIR=$(setup_temp_dir "create_bams")

# ===== FUNCTION DEFINITIONS =====

# Create BAM file for a single sample
create_bam_for_sample() {
    local assembly="$1"
    local read1="$2"
    local read2="$3"
    local output_dir="$4"
    local sample_name="$5"
    local index_base="$6"

    log "Creating BAM for sample $sample_name..."

    # Check if BAM already exists (skip if --force flag used)
    if [ "${FORCE_RUN:-false}" != "true" ] && \
       [ -f "${output_dir}/${sample_name}.sorted.bam" ] && \
       [ -f "${output_dir}/${sample_name}.sorted.bam.bai" ]; then
        log "BAM already exists for $sample_name, skipping..."
        return 0
    fi

    if [ "${FORCE_RUN:-false}" = "true" ] && [ -f "${output_dir}/${sample_name}.sorted.bam" ]; then
        log "Force mode: Deleting existing BAM for $sample_name"
        rm -f "${output_dir}/${sample_name}.sorted.bam" "${output_dir}/${sample_name}.sorted.bam.bai"
    fi

    # Activate environment with bowtie2 and samtools
    activate_env metawrap-env

    log "Mapping reads with bowtie2..."
    local sam_file="${TEMP_DIR}/${sample_name}.sam"
    local bam_file="${TEMP_DIR}/${sample_name}.bam"
    local sorted_bam="${output_dir}/${sample_name}.sorted.bam"

    bowtie2 \
        -x "$index_base" \
        -1 "$read1" \
        -2 "$read2" \
        -S "$sam_file" \
        -p ${SLURM_CPUS_PER_TASK:-64} \
        --sensitive \
        2>&1 | tee -a "${LOG_DIR}/${TREATMENT}_bam_mapping_${sample_name}.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Bowtie2 mapping failed for $sample_name"
        conda deactivate
        return 1
    fi

    log "Converting SAM to BAM and sorting..."
    samtools view -@ ${SLURM_CPUS_PER_TASK:-64} -bS "$sam_file" > "$bam_file"
    samtools sort -@ ${SLURM_CPUS_PER_TASK:-64} "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"

    # Clean up intermediate files
    rm -f "$sam_file" "$bam_file"

    conda deactivate

    log "✓ BAM file created: $sorted_bam"
    return 0
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

    log "====== Creating Shared BAM Files for Treatment: $TREATMENT ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    mkdir -p "$binning_dir"
    mkdir -p "$shared_bam_dir"

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

    # Build bowtie2 index
    log "Building bowtie2 index..."
    activate_env metawrap-env

    index_base="${TEMP_DIR}/assembly_index"
    bowtie2-build \
        --threads ${SLURM_CPUS_PER_TASK:-64} \
        "$assembly_fasta" \
        "$index_base" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_bowtie2_index.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Failed to build bowtie2 index"
        conda deactivate
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    conda deactivate

    log "✓ Bowtie2 index created"

    # Get all samples in this treatment and create BAM files
    log "Finding all samples in treatment $TREATMENT..."
    total_samples=$(get_total_samples)
    sample_count=0
    failed_count=0

    for i in $(seq 0 $((total_samples - 1))); do
        sample_info=$(get_sample_info_by_index $i 2>/dev/null)
        if [ -n "$sample_info" ]; then
            IFS='|' read -r sample_name treatment r1_path r2_path <<< "$sample_info"

            if [ "$treatment" = "$TREATMENT" ]; then
                log "Processing sample $sample_name ($((sample_count + 1)))..."

                if create_bam_for_sample "$assembly_fasta" "$r1_path" "$r2_path" "$shared_bam_dir" "$sample_name" "$index_base"; then
                    ((sample_count++))
                else
                    log "WARNING: Failed to create BAM for $sample_name"
                    ((failed_count++))
                fi
            fi
        fi
    done

    log "====== BAM Creation Summary ======"
    log "Successful: $sample_count samples"
    log "Failed: $failed_count samples"
    log "Total: $((sample_count + failed_count)) samples"

    if [ $sample_count -eq 0 ]; then
        log "ERROR: No BAM files were created"
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

    log "====== Creating Shared BAM Files for $SAMPLE_NAME ($TREATMENT) ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    mkdir -p "$binning_dir"
    mkdir -p "$shared_bam_dir"

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

    # Build bowtie2 index
    log "Building bowtie2 index..."
    activate_env metawrap-env

    index_base="${TEMP_DIR}/assembly_index"
    bowtie2-build \
        --threads ${SLURM_CPUS_PER_TASK:-64} \
        "$assembly_fasta" \
        "$index_base" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_bowtie2_index.log"

    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        log "ERROR: Failed to build bowtie2 index"
        conda deactivate
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
    conda deactivate

    log "✓ Bowtie2 index created"

    # Create BAM file for this sample
    if create_bam_for_sample "$assembly_fasta" "$R1_PATH" "$R2_PATH" "$shared_bam_dir" "$SAMPLE_NAME" "$index_base"; then
        log "✓ BAM file created successfully"
    else
        log "ERROR: Failed to create BAM file"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

fi

cleanup_temp_dir "$TEMP_DIR"
log "====== Shared BAM creation completed successfully ======"
