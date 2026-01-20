#!/usr/bin/env bash
#SBATCH --job-name=binning
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --account=ofinkel

# 03_binning.sh - Metagenomic binning using MetaBAT2, MaxBin2, and CONCOCT
# This script runs the 3 binners that were previously wrapped by MetaWRAP binning module
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
TEMP_DIR=$(setup_temp_dir "binning")

# ===== FUNCTION DEFINITIONS =====

# Run MetaBAT2
run_metabat2() {
    local assembly="$1"
    local depth_file="$2"
    local output_dir="$3"

    log "Running MetaBAT2..."

    activate_env metawrap-env

    if ! command -v metabat2 &> /dev/null; then
        log "ERROR: metabat2 not available"
        conda deactivate
        return 1
    fi

    local metabat_dir="${output_dir}/metabat2_bins"
    mkdir -p "$metabat_dir"

    # MetaBAT2 outputs bins with a prefix
    metabat2 \
        -i "$assembly" \
        -a "$depth_file" \
        -o "${metabat_dir}/bin" \
        -t ${SLURM_CPUS_PER_TASK:-64} \
        -m 1500 \
        --seed 1 \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_metabat2.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    if [ $exit_code -eq 0 ]; then
        local bin_count=$(ls -1 "${metabat_dir}"/bin.*.fa 2>/dev/null | wc -l)
        log "✓ MetaBAT2 completed: $bin_count bins"
        return 0
    else
        log "✗ MetaBAT2 failed"
        return 1
    fi
}

# Run MaxBin2
run_maxbin2() {
    local assembly="$1"
    local abundance_file="$2"
    local output_dir="$3"

    log "Running MaxBin2..."

    activate_env metawrap-env

    if ! command -v run_MaxBin.pl &> /dev/null; then
        log "ERROR: MaxBin2 not available"
        conda deactivate
        return 1
    fi

    local maxbin_dir="${output_dir}/maxbin2_bins"
    mkdir -p "$maxbin_dir"

    # MaxBin2 needs abundance file in specific format
    run_MaxBin.pl \
        -contig "$assembly" \
        -abund "$abundance_file" \
        -out "${maxbin_dir}/bin" \
        -thread ${SLURM_CPUS_PER_TASK:-64} \
        -min_contig_length 1500 \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_maxbin2.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    if [ $exit_code -eq 0 ]; then
        local bin_count=$(ls -1 "${maxbin_dir}"/bin.*.fasta 2>/dev/null | wc -l)
        # Rename .fasta to .fa for consistency
        for f in "${maxbin_dir}"/bin.*.fasta; do
            [ -f "$f" ] && mv "$f" "${f%.fasta}.fa"
        done
        log "✓ MaxBin2 completed: $bin_count bins"
        return 0
    else
        log "✗ MaxBin2 failed"
        return 1
    fi
}

# Run CONCOCT
run_concoct() {
    local assembly="$1"
    local bam_files="$2"  # Space-separated list of BAM files
    local output_dir="$3"

    log "Running CONCOCT..."

    activate_env metawrap-env

    if ! command -v concoct &> /dev/null; then
        log "ERROR: CONCOCT not available"
        conda deactivate
        return 1
    fi

    local concoct_dir="${output_dir}/concoct_bins"
    local concoct_work="${TEMP_DIR}/concoct_work"
    mkdir -p "$concoct_dir"
    mkdir -p "$concoct_work"

    # Step 1: Cut contigs into smaller pieces
    log "Cutting contigs for CONCOCT..."
    cut_up_fasta.py "$assembly" \
        -c 10000 \
        -o 0 \
        --merge_last \
        -b "${concoct_work}/contigs_10K.bed" \
        > "${concoct_work}/contigs_10K.fa"

    # Step 2: Generate coverage table from BAM files
    log "Generating coverage table..."
    concoct_coverage_table.py \
        "${concoct_work}/contigs_10K.bed" \
        $bam_files \
        > "${concoct_work}/coverage_table.tsv"

    # Step 3: Run CONCOCT
    log "Running CONCOCT clustering..."
    concoct \
        --composition_file "${concoct_work}/contigs_10K.fa" \
        --coverage_file "${concoct_work}/coverage_table.tsv" \
        -b "${concoct_work}/concoct_output" \
        -t ${SLURM_CPUS_PER_TASK:-64} \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_concoct.log"

    # Step 4: Merge subcontig clustering into original contigs
    log "Merging CONCOCT results..."
    merge_cutup_clustering.py \
        "${concoct_work}/concoct_output_clustering_gt1000.csv" \
        > "${concoct_work}/concoct_clustering_merged.csv"

    # Step 5: Extract bins as FASTA files
    log "Extracting CONCOCT bins..."
    extract_fasta_bins.py \
        "$assembly" \
        "${concoct_work}/concoct_clustering_merged.csv" \
        --output_path "$concoct_dir" 2>&1 | tail -20

    local exit_code=$?
    conda deactivate

    if [ $exit_code -eq 0 ]; then
        # Rename bins to match naming convention (bin.X.fa)
        local counter=1
        for f in "$concoct_dir"/*.fa; do
            if [ -f "$f" ] && [ "$(basename "$f")" != "bin.${counter}.fa" ]; then
                mv "$f" "${concoct_dir}/bin.${counter}.fa"
                ((counter++))
            fi
        done
        local bin_count=$(ls -1 "${concoct_dir}"/bin.*.fa 2>/dev/null | wc -l)
        log "✓ CONCOCT completed: $bin_count bins"
        return 0
    else
        log "✗ CONCOCT failed: extract_fasta_bins.py returned exit code $exit_code"
        log "  Check the output above for error details"
        return 1
    fi
}

# Generate depth file for MetaBAT2 from BAM files
generate_depth_file() {
    local assembly="$1"
    local bam_dir="$2"
    local output_file="$3"

    log "Generating depth file for MetaBAT2..."

    activate_env metawrap-env

    # Find all BAM files
    local bam_files=($(ls "${bam_dir}"/*.sorted.bam 2>/dev/null))

    if [ ${#bam_files[@]} -eq 0 ]; then
        log "ERROR: No BAM files found in $bam_dir"
        conda deactivate
        return 1
    fi

    log "Found ${#bam_files[@]} BAM files"

    # Use jgi_summarize_bam_contig_depths from MetaBAT2
    jgi_summarize_bam_contig_depths \
        --outputDepth "$output_file" \
        "${bam_files[@]}" \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}_depth_calculation.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    if [ $exit_code -eq 0 ] && [ -f "$output_file" ]; then
        log "✓ Depth file created: $output_file"
        return 0
    else
        log "✗ Failed to generate depth file"
        return 1
    fi
}

# Generate abundance file for MaxBin2 from depth file
generate_abundance_file() {
    local depth_file="$1"
    local output_file="$2"

    log "Generating abundance file for MaxBin2..."

    # MaxBin2 needs simple tab-separated format: contig_name <tab> abundance
    # We'll use the mean coverage across all samples
    awk 'NR>1 {sum=0; for(i=4;i<=NF;i++) sum+=$i; print $1"\t"sum/(NF-3)}' "$depth_file" > "$output_file"

    if [ -f "$output_file" ] && [ -s "$output_file" ]; then
        log "✓ Abundance file created: $output_file"
        return 0
    else
        log "✗ Failed to generate abundance file"
        return 1
    fi
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

    log "====== Starting Binning for Treatment: $TREATMENT ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    # Check if shared BAM directory exists
    if [ ! -d "$shared_bam_dir" ]; then
        log "ERROR: Shared BAM directory not found: $shared_bam_dir"
        log "Please run 03a_create_shared_bams.sh first"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
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

    # Check for existing bins
    if [ -d "${binning_dir}/metabat2_bins" ] && \
       [ -d "${binning_dir}/maxbin2_bins" ] && \
       [ -d "${binning_dir}/concoct_bins" ]; then
        log "All three binners already completed for $TREATMENT, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Generate depth file from shared BAMs
    depth_file="${TEMP_DIR}/depth.txt"
    if ! generate_depth_file "$assembly_fasta" "$shared_bam_dir" "$depth_file"; then
        log "ERROR: Failed to generate depth file"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Generate abundance file from depth file
    abundance_file="${TEMP_DIR}/abundance.txt"
    if ! generate_abundance_file "$depth_file" "$abundance_file"; then
        log "ERROR: Failed to generate abundance file"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Get list of BAM files for CONCOCT
    bam_files=$(ls "${shared_bam_dir}"/*.sorted.bam 2>/dev/null | tr '\n' ' ')

    # Run the three binners
    success_count=0

    if [ ! -d "${binning_dir}/metabat2_bins" ]; then
        if run_metabat2 "$assembly_fasta" "$depth_file" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "MetaBAT2 bins already exist, skipping..."
        ((success_count++))
    fi

    if [ ! -d "${binning_dir}/maxbin2_bins" ]; then
        if run_maxbin2 "$assembly_fasta" "$abundance_file" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "MaxBin2 bins already exist, skipping..."
        ((success_count++))
    fi

    if [ ! -d "${binning_dir}/concoct_bins" ]; then
        if run_concoct "$assembly_fasta" "$bam_files" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "CONCOCT bins already exist, skipping..."
        ((success_count++))
    fi

    log "====== Binning Summary ======"
    log "Successful binners: $success_count/3"

    if [ $success_count -eq 0 ]; then
        log "ERROR: All binners failed"
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

    log "====== Starting Binning for $SAMPLE_NAME ($TREATMENT) ======"

    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"
    shared_bam_dir="${binning_dir}/shared_bam_files"

    # Check if shared BAM directory exists
    if [ ! -d "$shared_bam_dir" ]; then
        log "ERROR: Shared BAM directory not found: $shared_bam_dir"
        log "Please run 03a_create_shared_bams.sh first"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
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

    # Check for existing bins
    if [ -d "${binning_dir}/metabat2_bins" ] && \
       [ -d "${binning_dir}/maxbin2_bins" ] && \
       [ -d "${binning_dir}/concoct_bins" ]; then
        log "All three binners already completed for $SAMPLE_NAME, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Generate depth file from shared BAMs
    depth_file="${TEMP_DIR}/depth.txt"
    if ! generate_depth_file "$assembly_fasta" "$shared_bam_dir" "$depth_file"; then
        log "ERROR: Failed to generate depth file"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Generate abundance file from depth file
    abundance_file="${TEMP_DIR}/abundance.txt"
    if ! generate_abundance_file "$depth_file" "$abundance_file"; then
        log "ERROR: Failed to generate abundance file"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Get list of BAM files for CONCOCT
    bam_files=$(ls "${shared_bam_dir}"/*.sorted.bam 2>/dev/null | tr '\n' ' ')

    # Run the three binners
    success_count=0

    if [ ! -d "${binning_dir}/metabat2_bins" ]; then
        if run_metabat2 "$assembly_fasta" "$depth_file" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "MetaBAT2 bins already exist, skipping..."
        ((success_count++))
    fi

    if [ ! -d "${binning_dir}/maxbin2_bins" ]; then
        if run_maxbin2 "$assembly_fasta" "$abundance_file" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "MaxBin2 bins already exist, skipping..."
        ((success_count++))
    fi

    if [ ! -d "${binning_dir}/concoct_bins" ]; then
        if run_concoct "$assembly_fasta" "$bam_files" "$binning_dir"; then
            ((success_count++))
        fi
    else
        log "CONCOCT bins already exist, skipping..."
        ((success_count++))
    fi

    log "====== Binning Summary ======"
    log "Successful binners: $success_count/3"

    if [ $success_count -eq 0 ]; then
        log "ERROR: All binners failed"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
fi

cleanup_temp_dir "$TEMP_DIR"
log "====== Binning completed successfully ======"
