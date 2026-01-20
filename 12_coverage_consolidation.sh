#!/bin/bash
#SBATCH --job-name=coverage_consolidation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00

# 12_coverage_consolidation.sh - Consolidate coverage data across all treatments
#
# This script creates:
# 1. A consolidated bin folder with all bins from all treatments
# 2. A combined CoverM abundance file in long format with unmapped reads
# 3. A renormalized abundance file with relative abundances scaled to 100%
# 4. Integration of GTDB-Tk taxonomy and Binette quality metrics

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# This stage runs once for the entire pipeline (not per treatment or sample)
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Coverage Consolidation Across All Treatments ======"

# Output directories
CONSOLIDATED_BINS_DIR="${OUTPUT_DIR}/consolidated_bins"
CONSOLIDATED_DATA_DIR="${OUTPUT_DIR}/consolidated_coverage_data"
mkdir -p "$CONSOLIDATED_BINS_DIR"
mkdir -p "$CONSOLIDATED_DATA_DIR"

# Function to collect and copy all bins from all treatments
collect_all_bins() {
    log "Collecting bins from all treatments..."

    # Get list of treatments
    local treatments_list=$(get_treatments)

    if [ -z "$treatments_list" ]; then
        log "ERROR: No treatments found"
        return 1
    fi

    local treatments=($treatments_list)
    log "  Found ${#treatments[@]} treatments: ${treatments[*]}"

    local total_bins=0

    # Process each treatment
    for treatment in "${treatments[@]}"; do
        log "  Processing treatment: $treatment"

        local selected_bins_base="${OUTPUT_DIR}/bin_collection/${treatment}/selected_bins"

        if [ ! -d "$selected_bins_base" ]; then
            log "    WARNING: Selected bins directory not found: $selected_bins_base"
            continue
        fi

        # Check if there are bins directly in the treatment folder (Coassembly Mode)
        local direct_bins=( "$selected_bins_base"/*.fa )

        # Check if there are sample sub-directories (Individual Mode)
        local sample_dirs=( "$selected_bins_base"/*/ )

        if [ -f "${direct_bins[0]}" ]; then
            # Coassembly mode - prefix with treatment name
            log "    Detected coassembly structure. Copying bins with treatment prefix..."
            for bin_file in "${direct_bins[@]}"; do
                local bin_name=$(basename "$bin_file" .fa)
                local new_name="${treatment}_${bin_name}.fa"
                cp "$bin_file" "${CONSOLIDATED_BINS_DIR}/${new_name}"
                ((total_bins++))
            done
            log "    Copied ${#direct_bins[@]} coassembled bins with prefix: ${treatment}_"

        elif [ -d "${sample_dirs[0]}" ]; then
            # Individual assembly mode - prefix with sample name
            log "    Detected individual assembly structure. Processing sample sub-directories..."
            for s_dir in "${sample_dirs[@]}"; do
                local s_name=$(basename "$s_dir")
                local s_bins=( "$s_dir"/*.fa )

                if [ -f "${s_bins[0]}" ]; then
                    log "      Copying bins from sample: $s_name"
                    for bin_file in "${s_bins[@]}"; do
                        local bin_name=$(basename "$bin_file" .fa)
                        local new_name="${s_name}_${bin_name}.fa"
                        cp "$bin_file" "${CONSOLIDATED_BINS_DIR}/${new_name}"
                        ((total_bins++))
                    done
                fi
            done
        else
            log "    WARNING: No .fa bins found in $selected_bins_base or its sub-directories"
        fi
    done

    log "  Total bins collected: $total_bins"
    echo "$total_bins"
    return 0
}

# Function to consolidate coverage data in long format with unmapped reads
consolidate_coverage_data() {
    log "Consolidating coverage data across all treatments..."

    # Create Python script for consolidation
    local python_script="${TEMP_DIR}/consolidate_coverage.py"

    cat > "$python_script" << 'PYTHON_EOF'
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path

output_dir = sys.argv[1]
consolidated_data_dir = sys.argv[2]

print(f"Output directory: {output_dir}")
print(f"Consolidated data directory: {consolidated_data_dir}")

# Get list of treatments from bin_collection directory
bin_collection_dir = f"{output_dir}/bin_collection"
if not os.path.exists(bin_collection_dir):
    print(f"ERROR: Bin collection directory not found: {bin_collection_dir}")
    sys.exit(1)

treatments = [d for d in os.listdir(bin_collection_dir)
              if os.path.isdir(os.path.join(bin_collection_dir, d))]

if not treatments:
    print("ERROR: No treatments found")
    sys.exit(1)

print(f"\nFound {len(treatments)} treatments: {', '.join(treatments)}")

# Collect all abundance data in long format
all_abundance_data = []

for treatment in treatments:
    print(f"\n{'='*60}")
    print(f"Processing treatment: {treatment}")
    print(f"{'='*60}")

    # Read CoverM abundance data
    abundance_file = f"{output_dir}/bin_collection/{treatment}/coverm_abundance_consolidated.tsv"
    if not os.path.exists(abundance_file):
        print(f"  WARNING: Abundance file not found: {abundance_file}")
        continue

    print(f"  Reading abundance data: {abundance_file}")
    df_abund = pd.read_csv(abundance_file, sep='\t')

    print(f"    Bins: {len(df_abund)}")
    print(f"    Columns: {list(df_abund.columns)}")

    # Melt abundance data to long format
    sample_cols = [col for col in df_abund.columns if col != 'Bin']
    df_long = df_abund.melt(id_vars=['Bin'],
                           value_vars=sample_cols,
                           var_name='Sample',
                           value_name='Relative_Abundance')
    df_long['Treatment'] = treatment

    # Add prefix to bin names to match consolidated bins folder
    # Determine if coassembly or individual mode
    selected_bins_base = f"{output_dir}/bin_collection/{treatment}/selected_bins"

    # Check for coassembly structure
    is_coassembly = any(
        f.endswith('.fa') for f in os.listdir(selected_bins_base)
        if os.path.isfile(os.path.join(selected_bins_base, f))
    ) if os.path.exists(selected_bins_base) else False

    if is_coassembly:
        # Coassembly mode - prefix with treatment
        df_long['Bin_Prefixed'] = treatment + '_' + df_long['Bin'].astype(str)
        print(f"    Coassembly mode: Added treatment prefix '{treatment}_'")
    else:
        # Individual mode - prefix with sample name
        # For individual mode, the original bin name should already have sample prefix
        # from stage 09 bin collection (e.g., "ILDRh2_bin.1")
        df_long['Bin_Prefixed'] = df_long['Sample'] + '_' + df_long['Bin'].astype(str)
        print(f"    Individual mode: Added sample prefix")

    all_abundance_data.append(df_long)
    print(f"  Processed {len(df_long)} abundance records")

# Combine all treatments
if not all_abundance_data:
    print("\nERROR: No abundance data collected")
    sys.exit(1)

df_combined = pd.concat(all_abundance_data, ignore_index=True)

print(f"\n{'='*60}")
print(f"Combined abundance data summary:")
print(f"{'='*60}")
print(f"  Total records: {len(df_combined)}")
print(f"  Unique bins (original): {df_combined['Bin'].nunique()}")
print(f"  Unique bins (prefixed): {df_combined['Bin_Prefixed'].nunique()}")
print(f"  Unique samples: {df_combined['Sample'].nunique()}")
print(f"  Treatments: {df_combined['Treatment'].nunique()}")

# Calculate unmapped reads per sample
print(f"\n{'='*60}")
print(f"Calculating unmapped reads per sample...")
print(f"{'='*60}")

# Group by treatment and sample, sum relative abundances
sample_totals = df_combined.groupby(['Treatment', 'Sample'])['Relative_Abundance'].sum().reset_index()
sample_totals.rename(columns={'Relative_Abundance': 'Mapped_Total'}, inplace=True)
sample_totals['Unmapped'] = 100.0 - sample_totals['Mapped_Total']

# Ensure unmapped is not negative (can happen due to rounding)
sample_totals['Unmapped'] = sample_totals['Unmapped'].clip(lower=0)

print("\nSample totals:")
print(sample_totals.to_string())

# Create unmapped rows
unmapped_rows = []
for _, row in sample_totals.iterrows():
    if row['Unmapped'] > 0.01:  # Only add if unmapped > 0.01%
        unmapped_rows.append({
            'Bin': 'unmapped',
            'Bin_Prefixed': 'unmapped',
            'Sample': row['Sample'],
            'Relative_Abundance': row['Unmapped'],
            'Treatment': row['Treatment']
        })

if unmapped_rows:
    df_unmapped = pd.DataFrame(unmapped_rows)
    df_with_unmapped = pd.concat([df_combined, df_unmapped], ignore_index=True)
    print(f"\nAdded {len(unmapped_rows)} unmapped records")
else:
    df_with_unmapped = df_combined.copy()
    print("\nNo unmapped reads to add (all samples have 100% mapped)")

# Now merge taxonomy data from GTDB-Tk
print(f"\n{'='*60}")
print(f"Merging GTDB-Tk taxonomy data...")
print(f"{'='*60}")

all_taxonomy_data = []

for treatment in treatments:
    bac_file = f"{output_dir}/bin_collection/{treatment}/gtdbtk/gtdbtk.bac120.summary.tsv"
    ar_file = f"{output_dir}/bin_collection/{treatment}/gtdbtk/gtdbtk.ar53.summary.tsv"

    tax_dfs = []
    for tax_file, domain in [(bac_file, 'Bacteria'), (ar_file, 'Archaea')]:
        if os.path.exists(tax_file):
            print(f"  Reading {domain} taxonomy from {treatment}: {tax_file}")
            df_tax = pd.read_csv(tax_file, sep='\t')
            if not df_tax.empty:
                df_tax['domain'] = domain
                df_tax['Treatment'] = treatment
                tax_dfs.append(df_tax)

    if tax_dfs:
        df_tax_combined = pd.concat(tax_dfs, ignore_index=True)
        all_taxonomy_data.append(df_tax_combined)
        print(f"  Found taxonomy for {len(df_tax_combined)} bins in {treatment}")

# Combine taxonomy data
if all_taxonomy_data:
    df_taxonomy_all = pd.concat(all_taxonomy_data, ignore_index=True)
    print(f"\nTotal taxonomy records: {len(df_taxonomy_all)}")

    # Clean bin names in taxonomy data (remove .fa extension)
    if 'user_genome' in df_taxonomy_all.columns:
        df_taxonomy_all['Bin_Original'] = df_taxonomy_all['user_genome'].str.replace('.fa', '', regex=False)

        # Add prefixed bin names to match abundance data
        # Check assembly mode for each treatment
        df_taxonomy_all['Bin_Prefixed'] = df_taxonomy_all.apply(
            lambda row: f"{row['Treatment']}_{row['Bin_Original']}", axis=1
        )

    # Merge taxonomy with abundance data on Bin_Prefixed
    df_merged = df_with_unmapped.merge(
        df_taxonomy_all[['Bin_Prefixed', 'classification', 'domain']],
        on='Bin_Prefixed',
        how='left'
    )

    print(f"\nMerged taxonomy: {len(df_merged)} records")
    print(f"  Bins with taxonomy: {df_merged['classification'].notna().sum()}")
else:
    print("\nWARNING: No taxonomy data found")
    df_merged = df_with_unmapped.copy()
    df_merged['classification'] = None
    df_merged['domain'] = None

# Merge Binette quality metrics
print(f"\n{'='*60}")
print(f"Merging Binette quality metrics...")
print(f"{'='*60}")

all_quality_data = []

for treatment in treatments:
    quality_file = f"{output_dir}/bin_refinement/{treatment}/binette/final_bins_quality_reports.tsv"

    if os.path.exists(quality_file):
        print(f"  Reading quality data from {treatment}: {quality_file}")
        df_quality = pd.read_csv(quality_file, sep='\t')
        if not df_quality.empty:
            df_quality['Treatment'] = treatment
            all_quality_data.append(df_quality)
            print(f"    Found quality metrics for {len(df_quality)} bins")
    else:
        print(f"  WARNING: Quality file not found: {quality_file}")

# Combine quality data
if all_quality_data:
    df_quality_all = pd.concat(all_quality_data, ignore_index=True)
    print(f"\nTotal quality records: {len(df_quality_all)}")

    # Clean bin names in quality data
    if 'Bin Id' in df_quality_all.columns:
        df_quality_all['Bin_Original'] = df_quality_all['Bin Id'].str.replace('.fa', '', regex=False)
    elif 'bin_id' in df_quality_all.columns:
        df_quality_all['Bin_Original'] = df_quality_all['bin_id'].str.replace('.fa', '', regex=False)
    else:
        # Try first column
        df_quality_all['Bin_Original'] = df_quality_all.iloc[:, 0].astype(str).str.replace('.fa', '', regex=False)

    # Add prefixed bin names
    df_quality_all['Bin_Prefixed'] = df_quality_all.apply(
        lambda row: f"{row['Treatment']}_{row['Bin_Original']}", axis=1
    )

    # Select relevant quality columns
    quality_cols = ['Bin_Prefixed']
    for col in ['completeness', 'Completeness', 'contamination', 'Contamination',
                'N50', 'size', 'Size', 'contigs', 'Contigs', 'N50_contigs']:
        if col in df_quality_all.columns:
            quality_cols.append(col)

    df_quality_subset = df_quality_all[quality_cols].copy()

    # Standardize column names
    df_quality_subset.columns = [col.lower() if col != 'Bin_Prefixed' else col
                                  for col in df_quality_subset.columns]

    # Merge quality with main dataframe
    df_final = df_merged.merge(df_quality_subset, on='Bin_Prefixed', how='left')

    print(f"\nMerged quality metrics: {len(df_final)} records")
    print(f"  Bins with quality data: {df_final['completeness'].notna().sum() if 'completeness' in df_final.columns else 0}")
else:
    print("\nWARNING: No quality data found")
    df_final = df_merged.copy()

# Reorder columns for better readability
cols_order = ['Treatment', 'Sample', 'Bin_Prefixed', 'Bin', 'Relative_Abundance']
other_cols = [col for col in df_final.columns if col not in cols_order]
df_final = df_final[cols_order + other_cols]

# Save long format file WITH unmapped reads
output_file_with_unmapped = f"{consolidated_data_dir}/combined_abundance_long_with_unmapped.tsv"
df_final.to_csv(output_file_with_unmapped, sep='\t', index=False, float_format='%.6f')
print(f"\n{'='*60}")
print(f"OUTPUT FILE 1: Long format with unmapped reads")
print(f"{'='*60}")
print(f"Saved to: {output_file_with_unmapped}")
print(f"  Total records: {len(df_final)}")
print(f"  Unique bins: {df_final['Bin_Prefixed'].nunique()}")
print(f"  Unique samples: {df_final['Sample'].nunique()}")

# Create renormalized version WITHOUT unmapped reads
print(f"\n{'='*60}")
print(f"Creating renormalized abundance file without unmapped...")
print(f"{'='*60}")

# Remove unmapped rows
df_no_unmapped = df_final[df_final['Bin_Prefixed'] != 'unmapped'].copy()

# Calculate normalization factors per sample
sample_sums = df_no_unmapped.groupby(['Treatment', 'Sample'])['Relative_Abundance'].sum().reset_index()
sample_sums.rename(columns={'Relative_Abundance': 'Sum_Abundance'}, inplace=True)

# Merge normalization factors
df_no_unmapped = df_no_unmapped.merge(sample_sums, on=['Treatment', 'Sample'], how='left')

# Renormalize to 100%
df_no_unmapped['Relative_Abundance_Renormalized'] = (
    df_no_unmapped['Relative_Abundance'] / df_no_unmapped['Sum_Abundance'] * 100.0
)

# Drop temporary column
df_no_unmapped.drop(columns=['Sum_Abundance'], inplace=True)

# Replace original relative abundance with renormalized
df_no_unmapped['Relative_Abundance'] = df_no_unmapped['Relative_Abundance_Renormalized']
df_no_unmapped.drop(columns=['Relative_Abundance_Renormalized'], inplace=True)

# Save renormalized file WITHOUT unmapped reads
output_file_renormalized = f"{consolidated_data_dir}/combined_abundance_long_renormalized.tsv"
df_no_unmapped.to_csv(output_file_renormalized, sep='\t', index=False, float_format='%.6f')

print(f"\n{'='*60}")
print(f"OUTPUT FILE 2: Long format renormalized (no unmapped)")
print(f"{'='*60}")
print(f"Saved to: {output_file_renormalized}")
print(f"  Total records: {len(df_no_unmapped)}")
print(f"  Unique bins: {df_no_unmapped['Bin_Prefixed'].nunique()}")
print(f"  Unique samples: {df_no_unmapped['Sample'].nunique()}")
print(f"  Sum of abundances per sample: 100.0% (renormalized)")

# Create summary report
print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"Treatments processed: {len(treatments)}")
print(f"Total unique bins: {df_final['Bin_Prefixed'].nunique() - 1}")  # -1 for unmapped
print(f"Total samples: {df_final['Sample'].nunique()}")
print(f"Bins with taxonomy: {df_final['classification'].notna().sum()}")
print(f"Bins with quality metrics: {df_final['completeness'].notna().sum() if 'completeness' in df_final.columns else 0}")

# Show first few rows
print(f"\nFirst 10 rows of combined data with unmapped:")
print(df_final.head(10).to_string())

print(f"\nFirst 10 rows of renormalized data:")
print(df_no_unmapped.head(10).to_string())

print(f"\n{'='*60}")
print(f"Coverage consolidation completed successfully!")
print(f"{'='*60}")

PYTHON_EOF

    # Run Python script
    log "  Running Python consolidation script..."

    # Check if pandas is available
    if python3 -c "import pandas" 2>/dev/null; then
        python3 "$python_script" "$OUTPUT_DIR" "$CONSOLIDATED_DATA_DIR" 2>&1 | tee "${LOG_DIR}/coverage_consolidation.log"
        local exit_code=${PIPESTATUS[0]}
    else
        log "ERROR: pandas not available in current Python environment"
        log "  Attempting to install pandas..."
        pip3 install pandas openpyxl 2>&1 | tail -20

        python3 "$python_script" "$OUTPUT_DIR" "$CONSOLIDATED_DATA_DIR" 2>&1 | tee "${LOG_DIR}/coverage_consolidation.log"
        local exit_code=${PIPESTATUS[0]}
    fi

    if [ $exit_code -eq 0 ]; then
        log "  Successfully consolidated coverage data"
        return 0
    else
        log "ERROR: Failed to consolidate coverage data"
        return 1
    fi
}

# Function to create summary report
create_consolidation_summary() {
    log "Creating consolidation summary report..."

    local summary_file="${CONSOLIDATED_DATA_DIR}/consolidation_summary.txt"

    # Count bins
    local bin_count=$(ls -1 "${CONSOLIDATED_BINS_DIR}"/*.fa 2>/dev/null | wc -l)

    # Count data records
    local records_with_unmapped=0
    local records_renormalized=0

    if [ -f "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_with_unmapped.tsv" ]; then
        records_with_unmapped=$(tail -n +2 "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_with_unmapped.tsv" | wc -l)
    fi

    if [ -f "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_renormalized.tsv" ]; then
        records_renormalized=$(tail -n +2 "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_renormalized.tsv" | wc -l)
    fi

    cat > "$summary_file" << EOF
Coverage and Bin Consolidation Summary
========================================

Date: $(date)
Pipeline: Metagenomics Assembly and Binning

Consolidated Outputs:
---------------------

1. Consolidated Bins Directory:
   ${CONSOLIDATED_BINS_DIR}/
   - Total bins collected: $bin_count
   - Bins are prefixed with treatment name (coassembly) or sample name (individual assembly)

2. Combined Abundance File (with unmapped reads):
   ${CONSOLIDATED_DATA_DIR}/combined_abundance_long_with_unmapped.tsv
   - Format: Long format (one row per bin-sample combination)
   - Records: $records_with_unmapped
   - Includes unmapped reads calculated as: 100% - sum(all bins per sample)
   - Contains GTDB-Tk taxonomy classification
   - Contains Binette quality metrics (completeness, contamination, etc.)

3. Renormalized Abundance File (without unmapped):
   ${CONSOLIDATED_DATA_DIR}/combined_abundance_long_renormalized.tsv
   - Format: Long format (one row per bin-sample combination)
   - Records: $records_renormalized
   - Relative abundances renormalized to 100% (excludes unmapped reads)
   - Contains same taxonomy and quality information as file #2

Data Columns:
-------------
- Treatment: Treatment name
- Sample: Sample name
- Bin_Prefixed: Bin name with treatment/sample prefix
- Bin: Original bin name (without prefix)
- Relative_Abundance: Relative abundance (%)
- classification: GTDB-Tk taxonomy classification string
- domain: Domain (Bacteria/Archaea)
- completeness: Bin completeness (%)
- contamination: Bin contamination (%)
- [Additional quality metrics from Binette as available]

Usage Notes:
------------
- Use the "with_unmapped" file to see the complete picture including unmapped reads
- Use the "renormalized" file for analyses requiring abundances to sum to 100%
- The Bin_Prefixed column ensures unique identification across all treatments/samples
- Unmapped reads represent reads that did not map to any of the recovered bins

Processing Details:
-------------------
- Source data: ${OUTPUT_DIR}/bin_collection/{TREATMENT}/
- Quality reports: ${OUTPUT_DIR}/bin_refinement/{TREATMENT}/binette/
- Coverage data: ${OUTPUT_DIR}/coverm/{TREATMENT}/{SAMPLE}/

EOF

    log "Summary report created: $summary_file"
}

# Validation function
validate_consolidation() {
    log "Validating coverage consolidation outputs..."

    # Check consolidated bins directory
    if [ ! -d "$CONSOLIDATED_BINS_DIR" ]; then
        log "Validation failed: Consolidated bins directory not found"
        return 1
    fi

    local bin_count=$(ls -1 "${CONSOLIDATED_BINS_DIR}"/*.fa 2>/dev/null | wc -l)
    if [ $bin_count -eq 0 ]; then
        log "Validation failed: No bins found in consolidated directory"
        return 1
    fi

    # Check output files
    if [ ! -f "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_with_unmapped.tsv" ]; then
        log "Validation failed: Combined abundance file (with unmapped) not found"
        return 1
    fi

    if [ ! -f "${CONSOLIDATED_DATA_DIR}/combined_abundance_long_renormalized.tsv" ]; then
        log "Validation failed: Renormalized abundance file not found"
        return 1
    fi

    log "Validation successful"
    log "  Consolidated bins: $bin_count"
    return 0
}

# Main processing function
main() {
    log "Starting coverage consolidation..."

    # Step 1: Collect all bins from all treatments
    local total_bins=$(collect_all_bins)
    if [ $? -ne 0 ] || [ "$total_bins" -eq 0 ]; then
        log "ERROR: Failed to collect bins or no bins found"
        return 1
    fi

    # Step 2: Consolidate coverage data with taxonomy and quality
    if ! consolidate_coverage_data; then
        log "ERROR: Failed to consolidate coverage data"
        return 1
    fi

    # Step 3: Create summary report
    create_consolidation_summary

    log "Coverage consolidation completed successfully"
    return 0
}

# Run main processing
if main; then
    if validate_consolidation; then
        log "====== Coverage Consolidation Completed Successfully ======"
        log "Output locations:"
        log "  Consolidated bins: $CONSOLIDATED_BINS_DIR"
        log "  Abundance files: $CONSOLIDATED_DATA_DIR"
    else
        log "ERROR: Validation failed"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi
else
    log "ERROR: Coverage consolidation failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
