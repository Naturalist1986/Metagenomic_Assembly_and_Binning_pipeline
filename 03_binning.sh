#!/bin/bash
#SBATCH --job-name=binning
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# 03_binning.sh - Metagenomic binning using MetaWRAP binning module

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

log "====== Starting Binning for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "binning"; then
    log "Binning already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_binning() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Running binning for $sample_name ($treatment)"
    
    local binning_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}"
    local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    
    mkdir -p "$binning_dir"
    
    # Check if already processed
    if ls "${binning_dir}"/*_bins/*.fa &>/dev/null 2>&1; then
        log "Sample $sample_name already binned, skipping..."
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
    log "Binning $num_contigs contigs for $sample_name"
    
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
    
# MetaWRAP expects files ending in 1.fastq and 2.fastq (no underscore, uncompressed)
local metawrap_r1="${metawrap_reads_dir}/${sample_name}_1.fastq"
local metawrap_r2="${metawrap_reads_dir}/${sample_name}_2.fastq"

log "Creating uncompressed reads for MetaWRAP..."
gunzip -c "$read1" > "$metawrap_r1"
gunzip -c "$read2" > "$metawrap_r2"

log "Created MetaWRAP-compatible read files:"
log "  R1: $metawrap_r1"
log "  R2: $metawrap_r2"
    
    log "Created MetaWRAP-compatible read files:"
    log "  R1: $metawrap_r1"
    log "  R2: $metawrap_r2"
    
    # Run MetaWRAP binning
    if run_metawrap_binning "$filtered_assembly" "$metawrap_r1" "$metawrap_r2" "$binning_dir"; then
        log "MetaWRAP binning completed successfully"
        create_binning_summary "$sample_name" "$treatment" "$binning_dir"
        return 0
    else
        log "ERROR: MetaWRAP binning failed for $sample_name"
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
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_metawrap_binning.log"
    
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
    local sample_name="$1"
    local treatment="$2"
    local binning_dir="$3"
    local summary_file="${binning_dir}/binning_summary.txt"
    
    cat > "$summary_file" << EOF
Binning Summary for $sample_name
================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment
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
            echo "  $(echo $binner | sed 's/_bins$/') bins:" >> "$summary_file"
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
    local sample_name="$1"
    local treatment="$2"
    local binning_dir="${OUTPUT_DIR}/binning/${treatment}/${sample_name}"
    
    # Check if at least one binner produced bins
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${binning_dir}/${binner}" ] && [ "$(ls -A "${binning_dir}/${binner}/"*.fa 2>/dev/null)" ]; then
            return 0
        fi
    done
    
    return 1
}

# Run the binning stage
if stage_binning "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_binning "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "binning"
        log "====== Binning completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Binning validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: Binning stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"