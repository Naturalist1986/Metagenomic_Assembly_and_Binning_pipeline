#!/bin/bash
#SBATCH --job-name=binning_treat
#SBATCH --array=0-9%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=48:00:00

# 03_binning_treatment_level.sh - Treatment-level metagenomic binning using MetaWRAP
# Bins co-assembled contigs using reads from all samples in the treatment

# CRITICAL: Ensure OUTPUT_DIR is set and normalize it (remove trailing slash)
if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR environment variable not set!" >&2
    exit 1
fi

# Remove trailing slashes
OUTPUT_DIR="${OUTPUT_DIR%/}"
export OUTPUT_DIR

# Set dependent paths
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"

# Source configuration and utilities FIRST (before setting file paths)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# NOW set these variables AFTER sourcing (so they don't get overwritten by empty initializations)
export TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
export SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"

# Debug: Print all paths
echo "========================================" >&2
echo "Environment Check:" >&2
echo "  OUTPUT_DIR: $OUTPUT_DIR" >&2
echo "  WORK_DIR: $WORK_DIR" >&2
echo "  TREATMENTS_FILE: $TREATMENTS_FILE" >&2
echo "  SAMPLE_INFO_FILE: $SAMPLE_INFO_FILE" >&2
echo "========================================" >&2

# Verify critical files exist
if [ ! -f "$TREATMENTS_FILE" ]; then
    echo "ERROR: Treatments file not found at: $TREATMENTS_FILE" >&2
    echo "Listing WORK_DIR contents:" >&2
    ls -la "$WORK_DIR" 2>&1 >&2 || echo "WORK_DIR doesn't exist!" >&2
    exit 1
fi

if [ ! -f "$SAMPLE_INFO_FILE" ]; then
    echo "ERROR: Sample info file not found at: $SAMPLE_INFO_FILE" >&2
    exit 1
fi

# Get treatment from array task ID
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

# Read treatments into array directly from file
mapfile -t TREATMENTS_ARRAY < "$TREATMENTS_FILE"

# Debug treatment array
echo "Treatment array info:" >&2
echo "  Total treatments: ${#TREATMENTS_ARRAY[@]}" >&2
echo "  Array index: $ARRAY_INDEX" >&2
for i in "${!TREATMENTS_ARRAY[@]}"; do
    echo "    [$i]: '${TREATMENTS_ARRAY[$i]}'" >&2
done

if [ $ARRAY_INDEX -ge ${#TREATMENTS_ARRAY[@]} ]; then
    echo "No treatment found for array index $ARRAY_INDEX" >&2
    exit 0
fi

TREATMENT="${TREATMENTS_ARRAY[$ARRAY_INDEX]}"
echo "Selected treatment: '$TREATMENT'" >&2
export TREATMENT

# Initialize
init_conda
create_treatment_dirs "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Treatment-Level Binning for Treatment: $TREATMENT ======"

# Check if treatment-level binning already completed
if check_treatment_checkpoint "$TREATMENT" "binning"; then
    log "Treatment-level binning already completed for $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main treatment-level binning function
stage_binning_treatment() {
    local treatment="$1"

    log "Running treatment-level binning for: $treatment"

    local binning_dir="${OUTPUT_DIR}/binning/${treatment}"
    local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"

    mkdir -p "$binning_dir"

    # Check if already processed
    if ls "${binning_dir}"/*_bins/*.fa &>/dev/null 2>&1; then
        log "Treatment $treatment already binned, skipping..."
        return 0
    fi

    # Check for co-assembly file
    local assembly_file=""
    for possible_file in \
        "${coassembly_dir}/contigs.fasta" \
        "${coassembly_dir}/scaffolds.fasta" \
        "${coassembly_dir}/final_assembly.fasta" \
        "${coassembly_dir}/final_contigs.fasta" \
        "${coassembly_dir}/assembly.fasta" \
        "${coassembly_dir}/${treatment}_contigs.fasta" \
        "${coassembly_dir}/${treatment}_scaffolds.fasta"; do

        if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
            assembly_file="$possible_file"
            log "Found co-assembly file: $assembly_file"
            break
        fi
    done

    if [ -z "$assembly_file" ]; then
        log "ERROR: No co-assembly file found for treatment $treatment"
        log "Expected location: $coassembly_dir"
        log "Checked files: contigs.fasta, scaffolds.fasta, final_assembly.fasta, final_contigs.fasta"
        return 1
    fi

    # Get all samples for this treatment
    local samples=($(get_samples_for_treatment "$treatment"))
    if [ ${#samples[@]} -eq 0 ]; then
        log "ERROR: No samples found for treatment $treatment"
        return 1
    fi

    log "Found ${#samples[@]} samples for treatment $treatment: ${samples[*]}"

    # Filter contigs by minimum length (1500 bp) for binning
    local filtered_assembly="${TEMP_DIR}/filtered_contigs.fasta"
    log "Filtering contigs to minimum length 1500 bp..."

    awk '/^>/ {if(seq) {if(length(seq) >= 1500) print header "\n" seq}; header=$0; seq=""; next} {seq=seq $0} END {if(seq && length(seq) >= 1500) print header "\n" seq}' "$assembly_file" > "$filtered_assembly"

    local num_contigs=$(grep -c "^>" "$assembly_file")
    local filtered_contigs=$(grep -c "^>" "$filtered_assembly" 2>/dev/null || echo "0")
    log "Contigs: $num_contigs total, $filtered_contigs after filtering (>= 1500 bp)"

    if [ $filtered_contigs -lt 10 ]; then
        log "ERROR: Too few contigs >= 1500 bp ($filtered_contigs) for meaningful binning"
        return 1
    fi

    # Prepare reads from ALL samples for MetaWRAP
    local metawrap_reads_dir="${TEMP_DIR}/metawrap_reads"
    mkdir -p "$metawrap_reads_dir"

    log "Preparing reads from all samples in treatment..."
    local read_pairs=()

    for sample in "${samples[@]}"; do
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample}"
        local read1="${quality_dir}/filtered_1.fastq.gz"
        local read2="${quality_dir}/filtered_2.fastq.gz"

        if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
            log "WARNING: Missing quality-filtered reads for sample $sample"
            log "  Expected: $read1 and $read2"
            continue
        fi

        # MetaWRAP expects files ending in _1.fastq and _2.fastq (uncompressed)
        local metawrap_r1="${metawrap_reads_dir}/${sample}_1.fastq"
        local metawrap_r2="${metawrap_reads_dir}/${sample}_2.fastq"

        log "  Decompressing reads for sample: $sample"
        gunzip -c "$read1" > "$metawrap_r1"
        gunzip -c "$read2" > "$metawrap_r2"

        read_pairs+=("$metawrap_r1" "$metawrap_r2")
    done

    if [ ${#read_pairs[@]} -eq 0 ]; then
        log "ERROR: No valid read pairs found for treatment $treatment"
        return 1
    fi

    log "Prepared ${#read_pairs[@]} read files (${#samples[@]} samples) for binning"

    # Run MetaWRAP binning with all samples
    if run_metawrap_binning_treatment "$filtered_assembly" "$binning_dir" "${read_pairs[@]}"; then
        log "MetaWRAP treatment-level binning completed successfully"
        create_binning_summary_treatment "$treatment" "$binning_dir"
        return 0
    else
        log "ERROR: MetaWRAP binning failed for treatment $treatment"
        return 1
    fi
}

# Run MetaWRAP binning module for treatment-level binning
run_metawrap_binning_treatment() {
    local assembly_file="$1"
    local output_dir="$2"
    shift 2
    local read_files=("$@")

    log "Running MetaWRAP binning module (treatment-level)..."

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
    log "metawrap binning -o $metawrap_output -t $SLURM_CPUS_PER_TASK -a $assembly_file --metabat2 --maxbin2 --concoct ${read_files[*]}"

    metawrap binning \
        -o "$metawrap_output" \
        -t $SLURM_CPUS_PER_TASK \
        -m $((${SLURM_MEM_PER_NODE:-256000} / 1000)) \
        -l 1500 \
        -a "$assembly_file" \
        --metabat2 --maxbin2 --concoct \
        "${read_files[@]}" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/metawrap_binning.log"

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

# Create binning summary for treatment-level binning
create_binning_summary_treatment() {
    local treatment="$1"
    local binning_dir="$2"
    local summary_file="${binning_dir}/binning_summary.txt"

    cat > "$summary_file" << EOF
Treatment-Level Binning Summary
================================

Date: $(date)
Treatment: $treatment

Binning Results:
EOF

    for binner in metabat2 maxbin2 concoct; do
        local binner_dir="${binning_dir}/${binner}_bins"
        if [ -d "$binner_dir" ]; then
            local num_bins=$(ls -1 "$binner_dir"/*.fa 2>/dev/null | wc -l)
            echo "  $binner: $num_bins bins" >> "$summary_file"
        else
            echo "  $binner: 0 bins" >> "$summary_file"
        fi
    done

    log "Created binning summary: $summary_file"
}

# Count bins across all binners
count_bins() {
    local binning_dir="$1"
    local total_bins=0

    for binner in metabat2 maxbin2 concoct; do
        local binner_dir="${binning_dir}/${binner}_bins"
        if [ -d "$binner_dir" ]; then
            local num_bins=$(ls -1 "$binner_dir"/*.fa 2>/dev/null | wc -l)
            total_bins=$((total_bins + num_bins))
        fi
    done

    echo "$total_bins"
}

# Main execution
stage_binning_treatment "$TREATMENT"

if [ $? -eq 0 ]; then
    create_treatment_checkpoint "$TREATMENT" "binning"
    log "====== Treatment-level binning completed for $TREATMENT ======"
    log "Total bins: $(count_bins "${OUTPUT_DIR}/binning/${TREATMENT}")"
else
    log "====== ERROR: Treatment-level binning failed for $TREATMENT ======"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
