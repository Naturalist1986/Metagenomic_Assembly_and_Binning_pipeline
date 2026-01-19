#!/bin/bash
#SBATCH --job-name=bin_collection
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# 09_bin_collection.sh - Collect selected bins, consolidate abundance, and run GTDB-Tk

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# This stage ALWAYS runs at treatment level (consolidates across samples)
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

# Initialize
init_conda
create_treatment_dirs "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Bin Collection for treatment $TREATMENT ======"

# Check if stage already completed (skip if --force flag used)
if [ "${FORCE_RUN:-false}" != "true" ] && check_treatment_checkpoint "$TREATMENT" "bin_collection"; then
    log "Bin collection already completed for treatment $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

if [ "${FORCE_RUN:-false}" = "true" ]; then
    log "Force mode enabled: Re-running bin collection even if previously completed"
fi

# Output directory for this treatment
OUTPUT_COLLECTION_DIR="${OUTPUT_DIR}/bin_collection/${TREATMENT}"
mkdir -p "$OUTPUT_COLLECTION_DIR"

# Function to collect selected bins from stage 7.5
# Function to collect selected bins (Updated to handle Individual and Coassembly modes)
collect_selected_bins() {
    local treatment="$1"
    local output_dir="$2"

    log "Collecting selected bins for treatment $treatment..."

    local selected_bins_base="${OUTPUT_DIR}/selected_bins/${treatment}"
    local collected_bins_dir="${output_dir}/selected_bins"
    mkdir -p "$collected_bins_dir"

    if [ ! -d "$selected_bins_base" ]; then
        log "ERROR: Selected bins directory not found: $selected_bins_base"
        return 1
    fi

    # Check if there are bins directly in the treatment folder (Coassembly Mode)
    local direct_bins=( "$selected_bins_base"/*.fa )
    
    # Check if there are sample sub-directories (Individual Mode)
    local sample_dirs=( "$selected_bins_base"/*/ )

    if [ -f "${direct_bins[0]}" ]; then
        log "  Detected coassembly structure. Copying bins..."
        cp "$selected_bins_base"/*.fa "$collected_bins_dir/"
    elif [ -d "${sample_dirs[0]}" ]; then
        log "  Detected individual assembly structure. Searching sample sub-directories..."
        # Iterate through each sample folder and copy bins
        for s_dir in "${sample_dirs[@]}"; do
            local s_name=$(basename "$s_dir")
            local s_bins=( "$s_dir"/*.fa )
            if [ -f "${s_bins[0]}" ]; then
                log "    Copying bins from sample: $s_name"
                # Optional: rename bins to include sample name to avoid collisions
                for bin_file in "${s_bins[@]}"; do
                    local b_name=$(basename "$bin_file")
                    cp "$bin_file" "${collected_bins_dir}/${s_name}_${b_name}"
                done
            fi
        done
    else
        log "ERROR: No .fa bins found in $selected_bins_base or its sub-directories"
        return 1
    fi

    local copied_count=$(ls -1 "$collected_bins_dir"/*.fa 2>/dev/null | wc -l)
    log "  Successfully collected $copied_count total bins"

    echo "$collected_bins_dir|$copied_count"
    return 0
}
# Function to run CoverM for each sample
# Function to run CoverM for each sample using existing BAM files
run_coverm_for_samples() {
    local treatment="$1"
    local samples_str="$2"
    local bins_dir="$3"

    log "Running CoverM abundance calculation using PRE-EXISTING BAM files for treatment $treatment..."

    # Define the custom BAM source directory
    # Note: Using the path provided, assuming 'treatment' in the path should be the actual treatment name
    local bam_source_dir="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/binning/${treatment}/shared_bam_files"
    
    log "  Source BAM directory: $bam_source_dir"

    # Convert samples string to array
    local samples=($samples_str)

    # Validate bins directory
    if [ ! -d "$bins_dir" ]; then
        log "ERROR: Bins directory not found: $bins_dir"
        return 1
    fi

    # Activate CoverM environment
    activate_env checkm

    # Process each sample
    local successful=0
    local failed=0

    # Get list of bin files
    local bin_files=("$bins_dir"/*.fa)

    for sample in "${samples[@]}"; do
        log "  Processing sample: $sample"

        # Output directory for this sample's results
        local sample_output_dir="${OUTPUT_DIR}/coverm/${treatment}/${sample}"
        mkdir -p "$sample_output_dir"

        # Locate the specific BAM file for this sample
        # We look for common patterns like sample.bam or sample_mapped.bam
        local bam_file=""
        if [ -f "${bam_source_dir}/${sample}.bam" ]; then
            bam_file="${bam_source_dir}/${sample}.bam"
        elif [ -f "${bam_source_dir}/${sample}_mapped.bam" ]; then
            bam_file="${bam_source_dir}/${sample}_mapped.bam"
        else
            # Try to find any bam file containing the sample name
            bam_file=$(find "$bam_source_dir" -maxdepth 1 -name "*${sample}*.bam" | head -n 1)
        fi

        if [ -z "$bam_file" ] || [ ! -f "$bam_file" ]; then
            log "    ERROR: BAM file not found for sample $sample in $bam_source_dir"
            ((failed++))
            continue
        fi

        log "    Using BAM file: $bam_file"

        # Run CoverM in genome mode using the BAM file (-b flag)
        coverm genome \
            --genome-fasta-files "${bin_files[@]}" \
            -b "$bam_file" \
            --output-file "${sample_output_dir}/abundance.tsv" \
            --output-format dense \
            --min-read-aligned-percent 0.75 \
            --min-read-percent-identity 0.95 \
            --min-covered-fraction 0 \
            --methods relative_abundance mean trimmed_mean covered_fraction reads_per_base rpkm tpm \
            --threads $SLURM_CPUS_PER_TASK \
            2>&1 | tee "${LOG_DIR}/${treatment}/${sample}_coverm.log"

        local exit_code=${PIPESTATUS[0]}

        if [ $exit_code -eq 0 ] && [ -f "${sample_output_dir}/abundance.tsv" ]; then
            log "    CoverM completed successfully for $sample"
            ((successful++))
        else
            log "    ERROR: CoverM failed for $sample (exit code: $exit_code)"
            ((failed++))
        fi
    done

    conda deactivate
    return 0
}

# Function to consolidate CoverM abundance files
consolidate_coverm_abundance() {
    local treatment="$1"
    local samples_str="$2"
    local output_dir="$3"

    log "Consolidating CoverM abundance data for treatment $treatment..."

    # Convert samples string to array
    local samples=($samples_str)

    # Output file (matches where Python saves it)
    local consolidated_file="${output_dir}/bin_collection/${treatment}/coverm_abundance_consolidated.tsv"

    # Create Python script for consolidation
    local python_script="${TEMP_DIR}/consolidate_coverm.py"

    cat > "$python_script" << 'PYTHON_EOF'
import sys
import os
import pandas as pd
from collections import defaultdict

treatment = sys.argv[1]
output_dir = sys.argv[2]
samples = sys.argv[3:]

print(f"Consolidating CoverM data for treatment: {treatment}")
print(f"Samples: {', '.join(samples)}")

# Dictionary to store abundance data: bin_name -> {sample_name: abundance}
bin_abundances = defaultdict(dict)
all_bins = set()

# Process each sample's CoverM file
for sample in samples:
    coverm_file = f"{output_dir}/coverm/{treatment}/{sample}/abundance.tsv"

    if not os.path.exists(coverm_file):
        print(f"WARNING: CoverM file not found for sample {sample}: {coverm_file}")
        continue

    print(f"Reading: {coverm_file}")

    try:
        # Read CoverM abundance file
        df = pd.read_csv(coverm_file, sep='\t')

        # Find the relative abundance column
        rel_abund_col = None
        for col in df.columns:
            if 'Relative Abundance' in col or 'relative_abundance' in col.lower():
                rel_abund_col = col
                break

        if rel_abund_col is None:
            print(f"WARNING: No relative abundance column found in {coverm_file}")
            # Try to use second column as fallback
            if len(df.columns) >= 2:
                rel_abund_col = df.columns[1]
                print(f"  Using column '{rel_abund_col}' as abundance")
            else:
                continue

        # Extract bin names and abundances
        for _, row in df.iterrows():
            bin_name = str(row.iloc[0]).strip()

            # Skip unmapped reads
            if bin_name.lower() == 'unmapped' or bin_name == '':
                continue

            # Remove .fa extension if present
            bin_name = bin_name.replace('.fa', '')

            # Get abundance value
            abundance = row[rel_abund_col]

            # Store abundance
            bin_abundances[bin_name][sample] = abundance
            all_bins.add(bin_name)

        print(f"  Processed {len(df)} bins from {sample}")

    except Exception as e:
        print(f"ERROR processing {coverm_file}: {e}")
        continue

# Create consolidated dataframe
if not all_bins:
    print("ERROR: No bins found in any CoverM files")
    sys.exit(1)

print(f"\nConsolidating {len(all_bins)} bins across {len(samples)} samples...")

# Create rows for dataframe
rows = []
for bin_name in sorted(all_bins):
    row = {'Bin': bin_name}
    for sample in samples:
        row[sample] = bin_abundances[bin_name].get(sample, 0.0)
    rows.append(row)

# Create dataframe
consolidated_df = pd.DataFrame(rows)

# Save to TSV
output_file = f"{output_dir}/bin_collection/{treatment}/coverm_abundance_consolidated.tsv"
consolidated_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
print(f"\nConsolidated abundance table saved to: {output_file}")
print(f"  Bins (rows): {len(consolidated_df)}")
print(f"  Samples (columns): {len(samples)}")

# Also save to Excel if openpyxl is available
try:
    excel_file = f"{output_dir}/bin_collection/{treatment}/coverm_abundance_consolidated.xlsx"
    consolidated_df.to_excel(excel_file, index=False, float_format='%.6f')
    print(f"Excel file saved to: {excel_file}")
except ImportError:
    print("Note: openpyxl not available, skipping Excel export")
except Exception as e:
    print(f"Warning: Could not save Excel file: {e}")

print("\nFirst 10 rows of consolidated data:")
print(consolidated_df.head(10).to_string())
PYTHON_EOF

    # Activate base/default conda environment for Python with pandas
    log "  Running Python consolidation script..."

    # Try to use pandas - check if available in base environment
    if python3 -c "import pandas" 2>/dev/null; then
        python3 "$python_script" "$treatment" "$OUTPUT_DIR" "${samples[@]}" 2>&1 | tee "${LOG_DIR}/${treatment}/coverm_consolidation.log"
        local exit_code=${PIPESTATUS[0]}
    else
        log "ERROR: pandas not available in current Python environment"
        log "  Attempting to install pandas..."
        pip3 install pandas openpyxl 2>&1 | tail -20

        python3 "$python_script" "$treatment" "$OUTPUT_DIR" "${samples[@]}" 2>&1 | tee "${LOG_DIR}/${treatment}/coverm_consolidation.log"
        local exit_code=${PIPESTATUS[0]}
    fi

    if [ $exit_code -eq 0 ] && [ -f "$consolidated_file" ]; then
        log "  Successfully consolidated CoverM abundance data"

        # Show summary
        local bin_count=$(tail -n +2 "$consolidated_file" | wc -l)
        local sample_count=$(($(head -1 "$consolidated_file" | tr '\t' '\n' | wc -l) - 1))
        log "  Consolidated table: $bin_count bins × $sample_count samples"

        return 0
    else
        log "ERROR: Failed to consolidate CoverM abundance data"
        return 1
    fi
}

# Function to run GTDB-Tk on selected bins
run_gtdbtk() {
    local treatment="$1"
    local bins_dir="$2"
    local output_dir="$3"

    log "Running GTDB-Tk for treatment $treatment..."
    log "  Bins directory: $bins_dir"

    # GTDB-Tk output directory
    local gtdbtk_dir="${output_dir}/gtdbtk"
    mkdir -p "$gtdbtk_dir"

    # Check if already completed
    if [ -f "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv" ] || [ -f "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv" ]; then
        log "  GTDB-Tk already completed for $treatment"
        return 0
    fi

    # Check that bins directory exists
    if [ ! -d "$bins_dir" ]; then
        log "ERROR: Bins directory not found: $bins_dir"
        return 1
    fi

    # Count bins - use find for more reliable counting
    local bin_count=$(find "$bins_dir" -maxdepth 1 -name "*.fa" -type f | wc -l)

    log "  Found $bin_count .fa files in bins directory"

    if [ $bin_count -eq 0 ]; then
        log "ERROR: No .fa bins found for GTDB-Tk analysis in $bins_dir"
        log "  Directory contents:"
        ls -la "$bins_dir" 2>&1 | head -20 | while read line; do log "    $line"; done
        return 1
    fi

    log "  Running GTDB-Tk classify_wf on $bin_count bins..."

    # Activate GTDB-Tk environment
    activate_env gtdbtk

    # Check if GTDB-Tk is available
    if ! command -v gtdbtk &> /dev/null; then
        log "ERROR: GTDB-Tk not available in environment"
        conda deactivate
        return 1
    fi

    # Set GTDB-Tk database path if not already set
    if [ -z "${GTDBTK_DATA_PATH:-}" ]; then
        export GTDBTK_DATA_PATH="/sci/backup/aerez/aerez/moshea/release226"
        log "  Setting GTDBTK_DATA_PATH to: $GTDBTK_DATA_PATH"
    fi

    if [ ! -d "$GTDBTK_DATA_PATH" ]; then
        log "ERROR: GTDB-Tk database not found at: $GTDBTK_DATA_PATH"
        conda deactivate
        return 1
    fi

    log "  Using GTDB-Tk database: $GTDBTK_DATA_PATH"

    # Mash database location (GTDB-Tk will create if doesn't exist)
    local mash_db="${GTDBTK_DATA_PATH}/mash_db.msh"
    log "  Using mash database: $mash_db (will be created if doesn't exist)"

    # Run GTDB-Tk classify workflow
    log "  Command: gtdbtk classify_wf --genome_dir $bins_dir --out_dir $gtdbtk_dir --extension fa --mash_db $mash_db --cpus $SLURM_CPUS_PER_TASK"

    gtdbtk classify_wf \
        --genome_dir "$bins_dir" \
        --out_dir "$gtdbtk_dir" \
        --extension fa \
        --mash_db "$mash_db" \
        --cpus $SLURM_CPUS_PER_TASK \
        2>&1 | tee "${LOG_DIR}/${treatment}/gtdbtk.log"

    local exit_code=${PIPESTATUS[0]}

    conda deactivate

    # Check if output was generated
    if [ $exit_code -eq 0 ]; then
        # Check for bacterial or archaeal results
        if [ -f "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv" ] || [ -f "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv" ]; then
            log "  GTDB-Tk completed successfully"

            # Log summary
            if [ -f "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv" ]; then
                local bac_count=$(tail -n +2 "${gtdbtk_dir}/gtdbtk.bac120.summary.tsv" | wc -l)
                log "    Bacterial genomes classified: $bac_count"
            fi

            if [ -f "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv" ]; then
                local ar_count=$(tail -n +2 "${gtdbtk_dir}/gtdbtk.ar53.summary.tsv" | wc -l)
                log "    Archaeal genomes classified: $ar_count"
            fi

            return 0
        else
            log "WARNING: GTDB-Tk completed but no summary files found"
            return 1
        fi
    else
        log "ERROR: GTDB-Tk failed with exit code: $exit_code"
        return 1
    fi
}

# Function to create consolidated summary report
create_summary_report() {
    local treatment="$1"
    local output_dir="$2"
    local bin_count="$3"
    local sample_count="$4"

    log "Creating summary report for treatment $treatment..."

    local summary_file="${output_dir}/collection_summary.txt"

    cat > "$summary_file" << EOF
Bin Collection Summary for Treatment: $treatment
================================================

Date: $(date)
Treatment: $treatment
Number of samples: $sample_count
Number of selected bins: $bin_count

Output Files:
-------------
1. Selected bins directory:
   ${output_dir}/selected_bins/
   - Contains all high-quality bins selected in stage 7.5

2. Consolidated CoverM abundance table:
   ${output_dir}/coverm_abundance_consolidated.tsv
   ${output_dir}/coverm_abundance_consolidated.xlsx
   - Bins as rows, samples as columns
   - Values represent relative abundance (%)

3. GTDB-Tk taxonomic classification:
   ${output_dir}/gtdbtk/
   - gtdbtk.bac120.summary.tsv (bacterial genomes)
   - gtdbtk.ar53.summary.tsv (archaeal genomes)
   - Full classification results and phylogenetic trees

Processing Details:
-------------------
- Selected bins source: ${OUTPUT_DIR}/selected_bins/${treatment}/
- CoverM data sources: ${OUTPUT_DIR}/coverm/${treatment}/[samples]/
- GTDB-Tk database: ${GTDBTK_DATA_PATH:-Not set}

Quality Standards (from stage 7.5):
-----------------------------------
Bins were selected based on the best version (orig/strict/permissive)
according to quality score: Completeness - (5 × Contamination)

EOF

    # Add GTDB-Tk summary if available
    if [ -f "${output_dir}/gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
        local bac_count=$(tail -n +2 "${output_dir}/gtdbtk/gtdbtk.bac120.summary.tsv" | wc -l)
        echo "Bacterial genomes classified: $bac_count" >> "$summary_file"
    fi

    if [ -f "${output_dir}/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
        local ar_count=$(tail -n +2 "${output_dir}/gtdbtk/gtdbtk.ar53.summary.tsv" | wc -l)
        echo "Archaeal genomes classified: $ar_count" >> "$summary_file"
    fi

    log "Summary report created: $summary_file"
}

# Main processing function
stage_bin_collection() {
    local treatment="$1"

    log "Running bin collection for treatment $treatment"

    # Step 1: Collect selected bins from stage 7.5
    local bins_info=$(collect_selected_bins "$treatment" "$OUTPUT_COLLECTION_DIR")
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to collect selected bins"
        return 1
    fi

    IFS='|' read -r collected_bins_dir bin_count <<< "$bins_info"

    # Step 2: Get all samples for this treatment
    log "Finding samples for treatment $treatment..."
    local samples=$(get_samples_for_treatment "$treatment")

    if [ -z "$samples" ]; then
        log "ERROR: No samples found for treatment $treatment"
        return 1
    fi

    local sample_count=$(echo "$samples" | wc -w)
    log "  Found $sample_count samples: $samples"

    # Step 3: Run CoverM for each sample using the collected bins
    if ! run_coverm_for_samples "$treatment" "$samples" "$collected_bins_dir"; then
        log "ERROR: Failed to run CoverM for samples"
        return 1
    fi

    # Step 4: Consolidate CoverM abundance data
    if ! consolidate_coverm_abundance "$treatment" "$samples" "$OUTPUT_DIR"; then
        log "ERROR: Failed to consolidate CoverM abundance data"
        return 1
    fi

    # Step 5: Run GTDB-Tk on selected bins
    if ! run_gtdbtk "$treatment" "$collected_bins_dir" "$OUTPUT_COLLECTION_DIR"; then
        log "ERROR: Failed to run GTDB-Tk"
        return 1
    fi

    # Step 6: Create summary report
    create_summary_report "$treatment" "$OUTPUT_COLLECTION_DIR" "$bin_count" "$sample_count"

    log "Bin collection completed for treatment $treatment"
    return 0
}

# Validation function
validate_bin_collection() {
    local treatment="$1"
    local output_dir="${OUTPUT_DIR}/bin_collection/${treatment}"

    log "Validating bin collection for treatment $treatment..."

    # Check that selected bins directory exists and has bins
    if [ ! -d "${output_dir}/selected_bins" ] || [ ! "$(ls -A "${output_dir}/selected_bins"/*.fa 2>/dev/null)" ]; then
        log "Validation failed: No selected bins found"
        return 1
    fi

    # Check that consolidated abundance file exists
    if [ ! -f "${output_dir}/coverm_abundance_consolidated.tsv" ]; then
        log "Validation failed: Consolidated abundance file not found"
        return 1
    fi

    # Check that GTDB-Tk results exist (at least one of bacterial or archaeal)
    if [ ! -f "${output_dir}/gtdbtk/gtdbtk.bac120.summary.tsv" ] && [ ! -f "${output_dir}/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
        log "Validation failed: GTDB-Tk results not found"
        return 1
    fi

    log "Validation successful"
    return 0
}

# Run the bin collection stage
if stage_bin_collection "$TREATMENT"; then
    # Validate results
    if validate_bin_collection "$TREATMENT"; then
        create_treatment_checkpoint "$TREATMENT" "bin_collection"
        log "====== Bin collection completed successfully for treatment $TREATMENT ======"
    else
        log "ERROR: Bin collection validation failed for treatment $TREATMENT"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
else
    log "ERROR: Bin collection stage failed for treatment $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
