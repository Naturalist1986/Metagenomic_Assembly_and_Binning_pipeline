#!/bin/bash
#SBATCH --job-name=plasmid_detection
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# 02_plasmid_detection.sh - Plasmid detection using PlasClass and MOB-suite

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

log "====== Starting Plasmid Detection for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "plasmid_detection"; then
    log "Plasmid detection already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Main processing function
stage_plasmid_detection() {
    local sample_name="$1"
    local treatment="$2"

    log "Running plasmid detection for $sample_name ($treatment)"
    log "Assembly mode: ${ASSEMBLY_MODE:-individual}"

    # Determine assembly directory based on assembly mode
    local assembly_dir
    if [ "${ASSEMBLY_MODE}" = "coassembly" ]; then
        # In coassembly mode, assembly is per treatment
        assembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
        log "Using coassembly directory: $assembly_dir"
    else
        # In individual mode, assembly is per sample
        assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
        log "Using individual assembly directory: $assembly_dir"
    fi

    local output_dir="${OUTPUT_DIR}/plasmids/${treatment}/${sample_name}"

    mkdir -p "$output_dir"

    # Check if already processed
    if [ -f "${output_dir}/plasmids_final.fasta" ] && [ -f "${output_dir}/non_plasmids.fasta" ]; then
        log "Sample $sample_name already processed, skipping..."
        return 0
    fi

    # Check for input assembly - try multiple possible names
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
        log "ERROR: No assembly found for $sample_name"
        log "  Checked locations in: $assembly_dir"
        log "    contigs.fasta"
        log "    scaffolds.fasta"
        log "    final_contigs.fasta"
        log "    assembly.fasta"
        return 1
    fi
    
    local num_contigs=$(grep -c "^>" "$assembly_file")
    log "Processing $num_contigs contigs for plasmid detection"
    
    # Run PlasClass
    local plasclass_success=false
    if run_plasclass "$assembly_file" "$output_dir"; then
        plasclass_success=true
    fi
    
    # Run MOB-suite
    local mobsuite_success=false
    if run_mobsuite "$assembly_file" "$output_dir"; then
        mobsuite_success=true
    fi
    
    # Combine results
    if [ "$plasclass_success" = true ] || [ "$mobsuite_success" = true ]; then
        combine_plasmid_results "$output_dir" "$assembly_file"
        log "Plasmid detection completed for $sample_name"
        return 0
    else
        log "ERROR: Both PlasClass and MOB-suite failed for $sample_name"
        # Create fallback output (all contigs as non-plasmids)
        create_fallback_output "$assembly_file" "$output_dir"
        log "WARNING: Created fallback output (no plasmids detected)"
        return 0
    fi
}

# Run PlasClass
run_plasclass() {
    local assembly_file="$1"
    local output_dir="$2"
    
    log "Running PlasClass for plasmid classification..."
    
    # Activate PlasClass environment
    activate_env scapp
    
    # Check if PlasClass is available
    if ! command -v classify_fasta.py &> /dev/null; then
        log "WARNING: PlasClass not available, skipping"
        conda deactivate
        return 1
    fi
    
    # Run PlasClass - uses default threshold
    classify_fasta.py \
        -f "$assembly_file" \
        -o "${output_dir}/plasclass_results.tsv" \
        -p ${SLURM_CPUS_PER_TASK:-4} \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_plasclass.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    # Extract plasmid sequences if results exist
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/plasclass_results.tsv" ]; then
        # Get plasmid contig IDs
        awk -F'\t' 'NR>1 && $2 == "plasmid" {print $1}' "${output_dir}/plasclass_results.tsv" > "${output_dir}/plasclass_plasmid_ids.txt"
        
        # Extract plasmid sequences
        if [ -s "${output_dir}/plasclass_plasmid_ids.txt" ]; then
            extract_sequences "$assembly_file" "${output_dir}/plasclass_plasmid_ids.txt" "${output_dir}/plasclass_plasmids.fasta"
        else
            touch "${output_dir}/plasclass_plasmids.fasta"
        fi
    fi
    
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/plasclass_results.tsv" ]; then
        log "PlasClass completed successfully"
        return 0
    else
        log "WARNING: PlasClass failed (exit code: $exit_code)"
        return 1
    fi
}

# Run MOB-suite
run_mobsuite() {
    local assembly_file="$1"
    local output_dir="$2"
    
    log "Running MOB-suite for plasmid detection..."
    
    # Activate MOB-suite environment
    activate_env mobsuite
    
    # Check if MOB-suite is available
    if ! command -v mob_recon &> /dev/null; then
        log "WARNING: MOB-suite not available, skipping"
        conda deactivate
        return 1
    fi
    
    # Create MOB-suite output directory
    local mob_dir="${output_dir}/mobsuite"
    
    # Remove existing directory if it exists
    if [ -d "$mob_dir" ]; then
        log "Removing existing MOB-suite output directory"
        rm -rf "$mob_dir"
    fi
    
    mkdir -p "$mob_dir"
    
    # Run MOB-recon - FIXED: added --force flag
    mob_recon \
        -i "$assembly_file" \
        -o "$mob_dir" \
        -n ${SLURM_CPUS_PER_TASK:-4} \
        --force \
        2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_mobsuite.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -d "$mob_dir" ]; then
        log "MOB-suite completed successfully"
        return 0
    else
        log "WARNING: MOB-suite failed (exit code: $exit_code)"
        return 1
    fi
}

# Create fallback output when both tools fail
create_fallback_output() {
    local assembly_file="$1"
    local output_dir="$2"
    
    log "Creating fallback output (treating all contigs as non-plasmids)"
    
    # Create empty plasmids file
    touch "${output_dir}/plasmids_final.fasta"
    
    # Copy all contigs to non-plasmids
    cp "$assembly_file" "${output_dir}/non_plasmids.fasta"
    
    # Create summary
    local total_contigs=$(grep -c "^>" "$assembly_file")
    create_plasmid_summary "$output_dir" "$total_contigs" "0" "$total_contigs"
}

# Combine plasmid detection results
combine_plasmid_results() {
    local output_dir="$1"
    local assembly_file="$2"
    
    log "Combining plasmid detection results..."
    
    # Create sets of plasmid and chromosome contigs
    local plasmid_contigs_file="${TEMP_DIR}/plasmid_contigs.txt"
    local chromosome_contigs_file="${TEMP_DIR}/chromosome_contigs.txt"
    
    # Initialize files
    > "$plasmid_contigs_file"
    > "$chromosome_contigs_file"
    
    # Get all contig IDs
    grep "^>" "$assembly_file" | sed 's/^>//' | cut -d' ' -f1 > "${TEMP_DIR}/all_contigs.txt"
    
    # Collect plasmid contigs from PlasClass
    if [ -f "${output_dir}/plasclass_results.tsv" ]; then
        log "Processing PlasClass results..."
        awk -F'\t' 'NR>1 && $2 == "plasmid" {print $1}' "${output_dir}/plasclass_results.tsv" >> "$plasmid_contigs_file"
    fi
    
    # Collect plasmid contigs from MOB-suite
    if [ -d "${output_dir}/mobsuite" ]; then
        log "Processing MOB-suite results..."
        # MOB-suite creates separate plasmid files
        find "${output_dir}/mobsuite" -name "plasmid_*.fasta" -exec grep "^>" {} \; 2>/dev/null | \
            sed 's/^>//' | cut -d' ' -f1 >> "$plasmid_contigs_file"
        
        # Also check for other MOB-suite output patterns
        if [ -f "${output_dir}/mobsuite/contig_report.txt" ]; then
            awk -F'\t' 'NR>1 && $3 == "plasmid" {print $1}' "${output_dir}/mobsuite/contig_report.txt" >> "$plasmid_contigs_file" 2>/dev/null || true
        fi
    fi
    
    # Remove duplicates and sort
    sort -u "$plasmid_contigs_file" > "${TEMP_DIR}/plasmid_contigs_unique.txt"
    
    # Get chromosome contigs (all contigs not identified as plasmids)
    comm -23 \
        <(sort "${TEMP_DIR}/all_contigs.txt") \
        <(sort "${TEMP_DIR}/plasmid_contigs_unique.txt") \
        > "$chromosome_contigs_file"
    
    # Extract sequences
    extract_sequences "$assembly_file" "${TEMP_DIR}/plasmid_contigs_unique.txt" "${output_dir}/plasmids_final.fasta"
    extract_sequences "$assembly_file" "$chromosome_contigs_file" "${output_dir}/non_plasmids.fasta"
    
    # Log statistics
    local total_contigs=$(wc -l < "${TEMP_DIR}/all_contigs.txt")
    local plasmid_contigs=$(wc -l < "${TEMP_DIR}/plasmid_contigs_unique.txt")
    local chromosome_contigs=$(wc -l < "$chromosome_contigs_file")
    
    log "Plasmid detection summary:"
    log "  Total contigs: $total_contigs"
    log "  Plasmid contigs: $plasmid_contigs"
    log "  Chromosome contigs: $chromosome_contigs"
    
    # Create summary file
    create_plasmid_summary "$output_dir" "$total_contigs" "$plasmid_contigs" "$chromosome_contigs"
}

# Extract sequences based on ID list
extract_sequences() {
    local input_fasta="$1"
    local id_list="$2"
    local output_fasta="$3"
    
    if [ ! -s "$id_list" ]; then
        # Create empty file if no IDs
        > "$output_fasta"
        return 0
    fi
    
    # Use seqtk if available, otherwise use awk
    if command -v seqtk &> /dev/null; then
        seqtk subseq "$input_fasta" "$id_list" > "$output_fasta"
    else
        # Fallback to awk
        awk -v ids="$id_list" '
            BEGIN {
                while ((getline id < ids) > 0) {
                    wanted[id] = 1
                }
                close(ids)
            }
            /^>/ {
                header = $0
                id = substr($1, 2)
                if (id in wanted) {
                    print_seq = 1
                    print header
                } else {
                    print_seq = 0
                }
                next
            }
            {
                if (print_seq) print
            }
        ' "$input_fasta" > "$output_fasta"
    fi
}

# Create plasmid detection summary
create_plasmid_summary() {
    local output_dir="$1"
    local total_contigs="$2"
    local plasmid_contigs="$3"
    local chromosome_contigs="$4"
    local summary_file="${output_dir}/plasmid_summary.txt"
    
    # Calculate percentage safely
    local plasmid_percentage="0.00"
    if [ "$total_contigs" -gt 0 ]; then
        plasmid_percentage=$(echo "scale=2; $plasmid_contigs * 100 / $total_contigs" | bc -l 2>/dev/null || echo "0.00")
    fi
    
    cat > "$summary_file" << EOF
Plasmid Detection Summary for $SAMPLE_NAME
==========================================

Date: $(date)
Sample: $SAMPLE_NAME
Treatment: $TREATMENT

Results:
  Total contigs: $total_contigs
  Plasmid contigs: $plasmid_contigs
  Chromosome contigs: $chromosome_contigs
  Plasmid percentage: ${plasmid_percentage}%

Output Files:
  Plasmids: ${output_dir}/plasmids_final.fasta
  Chromosomes: ${output_dir}/non_plasmids.fasta
  
Tools Used:
EOF
    
    if [ -f "${output_dir}/plasclass_results.tsv" ]; then
        echo "  - PlasClass: SUCCESS" >> "$summary_file"
    else
        echo "  - PlasClass: FAILED/SKIPPED" >> "$summary_file"
    fi
    
    if [ -d "${output_dir}/mobsuite" ]; then
        echo "  - MOB-suite: SUCCESS" >> "$summary_file"
    else
        echo "  - MOB-suite: FAILED/SKIPPED" >> "$summary_file"
    fi
    
    log "Plasmid detection summary created: $summary_file"
}

# Validation function
validate_plasmid_detection() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/plasmids/${treatment}/${sample_name}"
    
    # Check required output files exist
    if [ -f "${output_dir}/plasmids_final.fasta" ] && [ -f "${output_dir}/non_plasmids.fasta" ]; then
        # Check that non_plasmids.fasta is not empty (should have chromosome contigs)
        if [ -s "${output_dir}/non_plasmids.fasta" ]; then
            return 0
        fi
    fi
    
    return 1
}

# Run the plasmid detection stage
if stage_plasmid_detection "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_plasmid_detection "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "plasmid_detection"
        log "====== Plasmid detection completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Plasmid detection validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: Plasmid detection stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"