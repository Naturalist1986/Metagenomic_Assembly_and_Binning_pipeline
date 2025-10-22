#!/bin/bash
#SBATCH --job-name=binning_treatment
#SBATCH --array=0-9%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# 03_binning_treatment_level.sh - Treatment-level binning for both coassembly and individual assemblies

# CRITICAL: Ensure OUTPUT_DIR is set
if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR environment variable not set!" >&2
    exit 1
fi

OUTPUT_DIR="${OUTPUT_DIR%/}"
export OUTPUT_DIR
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"
export TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
export SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# Get treatment from array task ID (TREATMENT-LEVEL, not sample-level)
ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

# Read treatment directly from file
TREATMENT=$(sed -n "$((ARRAY_INDEX + 1))p" "$TREATMENTS_FILE" 2>/dev/null | tr -d '\r\n' | xargs)

if [ -z "$TREATMENT" ]; then
    echo "ERROR: No treatment found at index $ARRAY_INDEX" >&2
    exit 1
fi

echo "Processing treatment: '$TREATMENT' (index $ARRAY_INDEX)" >&2
export TREATMENT

# Initialize
init_conda
create_treatment_dirs "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Treatment-Level Binning for $TREATMENT ======"

# Check if treatment-level binning already completed
if check_treatment_checkpoint "$TREATMENT" "binning"; then
    log "Treatment-level binning already completed for $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Get all samples for this treatment
get_treatment_samples() {
    local treatment="$1"
    local samples=()
    
    # Read from sample info file
    while IFS='|' read -r sample_name sample_treatment r1 r2; do
        if [ "$sample_treatment" = "$treatment" ]; then
            samples+=("$sample_name")
        fi
    done < "$SAMPLE_INFO_FILE"
    
    echo "${samples[@]}"
}

# Collect all read files from treatment
collect_treatment_reads() {
    local treatment="$1"
    local temp_dir="$2"
    
    log "Collecting all read files from treatment: $treatment"
    
    local samples=($(get_treatment_samples "$treatment"))
    local read_pairs=()
    
    for sample in "${samples[@]}"; do
        # Check for validated reads first (from stage 0.5), fall back to quality-filtered
        local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample}"
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample}"
        
        local read1=""
        local read2=""
        
        if [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
            read1="${validated_dir}/validated_1.fastq.gz"
            read2="${validated_dir}/validated_2.fastq.gz"
            log "  Found validated reads for sample: $sample"
        elif [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
            read1="${quality_dir}/filtered_1.fastq.gz"
            read2="${quality_dir}/filtered_2.fastq.gz"
            log "  Found quality-filtered reads for sample: $sample"
        fi
        
        if [ -f "$read1" ] && [ -f "$read2" ]; then
            read_pairs+=("$sample|$read1|$read2")
        else
            log "  WARNING: Missing reads for sample: $sample"
        fi
    done
    
    if [ ${#read_pairs[@]} -eq 0 ]; then
        log "ERROR: No valid read pairs found for treatment $treatment"
        return 1
    fi
    
    log "Total samples with reads: ${#read_pairs[@]}"
    
    # Prepare uncompressed reads for MetaWRAP
    local metawrap_reads_dir="${temp_dir}/metawrap_reads"
    mkdir -p "$metawrap_reads_dir"
    
    local r1_files=()
    local r2_files=()
    
    for pair in "${read_pairs[@]}"; do
        IFS='|' read -r sample r1_gz r2_gz <<< "$pair"
        
        local r1_out="${metawrap_reads_dir}/${sample}_1.fastq"
        local r2_out="${metawrap_reads_dir}/${sample}_2.fastq"
        
        log "  Decompressing reads for $sample..."
        gunzip -c "$r1_gz" > "$r1_out"
        gunzip -c "$r2_gz" > "$r2_out"
        
        r1_files+=("$r1_out")
        r2_files+=("$r2_out")
    done
    
    # Return space-separated lists
    echo "${r1_files[*]}|${r2_files[*]}"
}

# Get assembly file for treatment (handles both coassembly and individual modes)
get_treatment_assembly() {
    local treatment="$1"
    
    # Check assembly mode
    if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        # Look for co-assembly
        local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
        
        log "Looking for co-assembly in: $coassembly_dir"
        
        for possible_file in \
            "${coassembly_dir}/contigs.fasta" \
            "${coassembly_dir}/scaffolds.fasta" \
            "${coassembly_dir}/final_assembly.fasta"; do
            
            if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                log "Found co-assembly: $possible_file"
                echo "$possible_file"
                return 0
            fi
        done
        
        log "ERROR: No co-assembly found for treatment $treatment"
        return 1
    else
        # Individual assembly mode - need to concatenate or use first sample
        log "Individual assembly mode - using first sample's assembly"
        
        local samples=($(get_treatment_samples "$treatment"))
        
        if [ ${#samples[@]} -eq 0 ]; then
            log "ERROR: No samples found for treatment $treatment"
            return 1
        fi
        
        # Use first sample's assembly (or could concatenate all)
        local first_sample="${samples[0]}"
        local sample_assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${first_sample}"
        
        for possible_file in \
            "${sample_assembly_dir}/contigs.fasta" \
            "${sample_assembly_dir}/scaffolds.fasta" \
            "${sample_assembly_dir}/final_assembly.fasta"; do
            
            if [ -f "$possible_file" ] && [ -s "$possible_file" ]; then
                log "Using assembly from sample: $first_sample"
                log "WARNING: Only using one sample's assembly. Co-assembly recommended for treatment-level binning."
                echo "$possible_file"
                return 0
            fi
        done
        
        log "ERROR: No assembly found for sample $first_sample"
        return 1
    fi
}

# Main processing function
stage_binning() {
    local treatment="$1"
    
    log "Running treatment-level binning for $treatment"
    log "Assembly mode: $ASSEMBLY_MODE"
    
    local binning_dir="${OUTPUT_DIR}/binning/${treatment}"
    mkdir -p "$binning_dir"
    
    # Check if already processed
    if ls "${binning_dir}"/*_bins/*.fa &>/dev/null 2>&1; then
        log "Treatment $treatment already binned, skipping..."
        return 0
    fi
    
    # Get assembly file
    local assembly_file=$(get_treatment_assembly "$treatment")
    if [ -z "$assembly_file" ]; then
        log "ERROR: No assembly file found for treatment $treatment"
        return 1
    fi
    log "Using assembly file: $assembly_file"
    
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
    
    # Collect all reads from treatment
    local read_info=$(collect_treatment_reads "$treatment" "$TEMP_DIR")
    if [ -z "$read_info" ]; then
        log "ERROR: Failed to collect treatment reads"
        return 1
    fi
    
    IFS='|' read -r r1_files r2_files <<< "$read_info"
    
    log "Read files collected for MetaWRAP:"
    log "  Forward reads: $r1_files"
    log "  Reverse reads: $r2_files"
    
    # Run MetaWRAP binning with all reads
    if run_metawrap_binning "$filtered_assembly" "$r1_files" "$r2_files" "$binning_dir"; then
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
    local r1_files="$2"
    local r2_files="$3"
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
    
    # Run MetaWRAP binning with all read files
    log "Running MetaWRAP binning command..."
    
    metawrap binning \
        -o "$metawrap_output" \
        -t $SLURM_CPUS_PER_TASK \
        -m $((${SLURM_MEM_PER_NODE:-128000} / 1000)) \
        -l 1500 \
        -a "$assembly_file" \
        --metabat2 --maxbin2 --concoct \
        $r1_files $r2_files \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_metawrap_binning.log"
    
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
    
    local samples=($(get_treatment_samples "$treatment"))
    
    cat > "$summary_file" << 'EOF'
Treatment-Level Binning Summary
================================

EOF
    
    # Add dynamic content
    cat >> "$summary_file" << EOF
Date: $(date)
Treatment: $treatment
Assembly mode: $ASSEMBLY_MODE
Samples included: ${#samples[@]}
Sample list: ${samples[@]}
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
            echo "  $(echo $binner | sed 's/_bins$//') bins:" >> "$summary_file"
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
        if [ -d "${binning_dir}/${binner}" ] && [ -n "$(ls -A "${binning_dir}/${binner}/"*.fa 2>/dev/null)" ]; then
            return 0
        fi
    done
    
    return 1
}

# Run the binning stage (treatment-level)
if stage_binning "$TREATMENT"; then
    # Validate results
    if validate_binning "$TREATMENT"; then
        create_treatment_checkpoint "$TREATMENT" "binning"
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