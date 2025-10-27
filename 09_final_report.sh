#!/bin/bash
#SBATCH --job-name=final_report
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# 11_final_report.sh - Generate comprehensive final report with taxonomy-labeled abundance plots

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# This stage runs once for the entire pipeline (not per treatment or sample)
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Final Report Generation ======"

# Output directory for final report
FINAL_REPORT_DIR="${OUTPUT_DIR}/final_report"
mkdir -p "$FINAL_REPORT_DIR"

# Function to collect all treatment data
collect_all_treatment_data() {
    log "Collecting data from all treatments..."

    # Get list of treatments
    local treatments_list=$(get_treatments)

    if [ -z "$treatments_list" ]; then
        log "ERROR: No treatments found"
        return 1
    fi

    local treatments=($treatments_list)
    log "  Found ${#treatments[@]} treatments: ${treatments[*]}"

    # Check that bin collection stage completed for all treatments
    local missing_treatments=()
    for treatment in "${treatments[@]}"; do
        local collection_dir="${OUTPUT_DIR}/bin_collection/${treatment}"
        if [ ! -d "$collection_dir" ] || \
           [ ! -f "${collection_dir}/coverm_abundance_consolidated.tsv" ]; then
            missing_treatments+=("$treatment")
        fi
    done

    if [ ${#missing_treatments[@]} -gt 0 ]; then
        log "ERROR: Bin collection not complete for treatments: ${missing_treatments[*]}"
        log "  Please run stage 10 (bin collection) for all treatments first"
        return 1
    fi

    log "  All treatments have completed bin collection"
    echo "${treatments[*]}"
    return 0
}

# Function to create final report plots
create_final_report_plots() {
    local treatments_str="$1"

    log "Creating final report plots..."

    # Create Python plotting script
    local python_script="${TEMP_DIR}/create_final_report.py"

    cat > "$python_script" << 'PYTHON_EOF'
import sys
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
sns.set_palette("husl")

output_dir = sys.argv[1]
treatments = sys.argv[2:]

print(f"Processing {len(treatments)} treatments: {', '.join(treatments)}")
print(f"Output directory: {output_dir}")

# Collect data from all treatments
all_abundance_data = []
all_taxonomy_data = []

for treatment in treatments:
    print(f"\nProcessing treatment: {treatment}")

    # Read CoverM abundance data
    abundance_file = f"{output_dir}/bin_collection/{treatment}/coverm_abundance_consolidated.tsv"
    if not os.path.exists(abundance_file):
        print(f"  WARNING: Abundance file not found: {abundance_file}")
        continue

    print(f"  Reading abundance data: {abundance_file}")
    df_abund = pd.read_csv(abundance_file, sep='\t')

    # Melt abundance data to long format
    sample_cols = [col for col in df_abund.columns if col != 'Bin']
    df_long = df_abund.melt(id_vars=['Bin'], value_vars=sample_cols,
                           var_name='Sample', value_name='Relative_Abundance')
    df_long['Treatment'] = treatment
    all_abundance_data.append(df_long)

    # Read GTDB-Tk taxonomy data
    bac_file = f"{output_dir}/bin_collection/{treatment}/gtdbtk/gtdbtk.bac120.summary.tsv"
    ar_file = f"{output_dir}/bin_collection/{treatment}/gtdbtk/gtdbtk.ar53.summary.tsv"

    tax_dfs = []
    for tax_file, domain in [(bac_file, 'Bacteria'), (ar_file, 'Archaea')]:
        if os.path.exists(tax_file):
            print(f"  Reading {domain} taxonomy: {tax_file}")
            df_tax = pd.read_csv(tax_file, sep='\t')
            if not df_tax.empty:
                df_tax['domain'] = domain
                tax_dfs.append(df_tax)

    if tax_dfs:
        df_tax_combined = pd.concat(tax_dfs, ignore_index=True)
        df_tax_combined['Treatment'] = treatment
        all_taxonomy_data.append(df_tax_combined)
        print(f"  Found taxonomy for {len(df_tax_combined)} bins")
    else:
        print(f"  WARNING: No taxonomy data found for {treatment}")

# Combine all data
if not all_abundance_data:
    print("ERROR: No abundance data collected")
    sys.exit(1)

df_abundance_all = pd.concat(all_abundance_data, ignore_index=True)
print(f"\nTotal abundance records: {len(df_abundance_all)}")
print(f"Unique bins: {df_abundance_all['Bin'].nunique()}")
print(f"Unique samples: {df_abundance_all['Sample'].nunique()}")
print(f"Treatments: {df_abundance_all['Treatment'].nunique()}")

if all_taxonomy_data:
    df_taxonomy_all = pd.concat(all_taxonomy_data, ignore_index=True)
    print(f"Total taxonomy records: {len(df_taxonomy_all)}")
else:
    print("WARNING: No taxonomy data available")
    df_taxonomy_all = pd.DataFrame()

# Parse taxonomy to extract genus and species
def parse_taxonomy(classification_str):
    """Parse GTDB taxonomy string to extract genus and species"""
    if pd.isna(classification_str) or classification_str == '':
        return 'Unclassified', 'Unclassified'

    # GTDB format: d__domain;p__phylum;c__class;o__order;f__family;g__genus;s__species
    parts = classification_str.split(';')

    genus = 'Unclassified'
    species = 'Unclassified'

    for part in parts:
        if part.startswith('g__'):
            genus = part.replace('g__', '').strip()
            if genus == '':
                genus = 'Unclassified'
        elif part.startswith('s__'):
            species = part.replace('s__', '').strip()
            if species == '':
                species = 'Unclassified'

    return genus, species

if not df_taxonomy_all.empty and 'classification' in df_taxonomy_all.columns:
    df_taxonomy_all[['Genus', 'Species']] = df_taxonomy_all['classification'].apply(
        lambda x: pd.Series(parse_taxonomy(x))
    )

    # Create Bin column matching abundance data (remove .fa extension if present)
    if 'user_genome' in df_taxonomy_all.columns:
        df_taxonomy_all['Bin'] = df_taxonomy_all['user_genome'].str.replace('.fa', '', regex=False)

    # Merge taxonomy with abundance data
    df_merged = df_abundance_all.merge(
        df_taxonomy_all[['Bin', 'Treatment', 'Genus', 'Species', 'domain']],
        on=['Bin', 'Treatment'],
        how='left'
    )

    # Fill missing taxonomy
    df_merged['Genus'].fillna('Unclassified', inplace=True)
    df_merged['Species'].fillna('Unclassified', inplace=True)
    df_merged['domain'].fillna('Unknown', inplace=True)

    print(f"\nMerged data: {len(df_merged)} records")
    print(f"Bins with taxonomy: {df_merged['Genus'].ne('Unclassified').sum()}")
else:
    df_merged = df_abundance_all.copy()
    df_merged['Genus'] = 'Unclassified'
    df_merged['Species'] = 'Unclassified'
    df_merged['domain'] = 'Unknown'
    print("\nNo taxonomy data to merge")

# Create shortened labels for plotting
df_merged['Genus_Label'] = df_merged['Genus'].apply(
    lambda x: x if x == 'Unclassified' else x.split('_')[-1] if '_' in x else x
)

df_merged['Species_Label'] = df_merged['Species'].apply(
    lambda x: x if x == 'Unclassified' else x.split('_')[-1] if '_' in x else x
)

# Save merged data
merged_file = f"{output_dir}/final_report/merged_abundance_taxonomy.tsv"
df_merged.to_csv(merged_file, sep='\t', index=False)
print(f"\nSaved merged data to: {merged_file}")

# Create plots directory
plots_dir = f"{output_dir}/final_report/plots"
os.makedirs(plots_dir, exist_ok=True)

print("\n" + "="*80)
print("Creating plots...")
print("="*80)

# Function to create stacked barplot
def create_stacked_barplot(data, group_by, label_col, title, filename,
                          include_unmapped=False, figsize=(14, 8)):
    """Create stacked barplot with taxonomy labels"""

    print(f"\nCreating plot: {title}")
    print(f"  Grouping by: {group_by}")
    print(f"  Labels: {label_col}")
    print(f"  Include unmapped: {include_unmapped}")

    # Prepare data
    plot_data = data.copy()

    # Group data for stacking
    grouped = plot_data.groupby([group_by, label_col])['Relative_Abundance'].sum().reset_index()

    # Pivot for stacking
    pivot_data = grouped.pivot(index=group_by, columns=label_col, values='Relative_Abundance')
    pivot_data = pivot_data.fillna(0)

    # Add unmapped reads if requested
    if include_unmapped:
        # Calculate unmapped reads as 100% minus the sum of all bin abundances in each sample
        # This represents reads that didn't map to any of the assembled bins
        total_mapped = pivot_data.sum(axis=1)  # Sum across all bins per sample
        pivot_data['Unmapped'] = 100 - total_mapped
        pivot_data['Unmapped'] = pivot_data['Unmapped'].clip(lower=0)  # Ensure non-negative

    # Sort columns by abundance (most abundant first)
    col_sums = pivot_data.sum(axis=0).sort_values(ascending=False)
    pivot_data = pivot_data[col_sums.index]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Create color palette
    n_colors = len(pivot_data.columns)
    if 'Unmapped' in pivot_data.columns:
        # Use grey for unmapped, colors for rest
        colors = sns.color_palette("husl", n_colors - 1)
        colors.append((0.7, 0.7, 0.7))  # Grey for unmapped
    else:
        colors = sns.color_palette("husl", n_colors)

    # Create stacked bar plot
    pivot_data.plot(kind='bar', stacked=True, ax=ax, color=colors, width=0.8)

    # Customize plot
    ax.set_xlabel(group_by, fontsize=12, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_ylim(0, 100)

    # Rotate x-axis labels
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Customize legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
             frameon=True, fontsize=9, title=label_col)

    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Tight layout
    plt.tight_layout()

    # Save
    output_path = f"{plots_dir}/{filename}"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()

# Plot 1: Per-treatment stacked barplots with genus labels (no unmapped)
for treatment in treatments:
    treatment_data = df_merged[df_merged['Treatment'] == treatment]
    if not treatment_data.empty:
        create_stacked_barplot(
            treatment_data,
            group_by='Sample',
            label_col='Genus_Label',
            title=f'Bin Abundance by Genus - {treatment}',
            filename=f'abundance_genus_{treatment}.png',
            include_unmapped=False
        )

# Plot 2: Per-treatment stacked barplots with genus labels (with unmapped)
for treatment in treatments:
    treatment_data = df_merged[df_merged['Treatment'] == treatment]
    if not treatment_data.empty:
        create_stacked_barplot(
            treatment_data,
            group_by='Sample',
            label_col='Genus_Label',
            title=f'Bin Abundance by Genus (with unmapped) - {treatment}',
            filename=f'abundance_genus_unmapped_{treatment}.png',
            include_unmapped=True
        )

# Plot 3: Per-treatment stacked barplots with species labels (no unmapped)
for treatment in treatments:
    treatment_data = df_merged[df_merged['Treatment'] == treatment]
    if not treatment_data.empty:
        create_stacked_barplot(
            treatment_data,
            group_by='Sample',
            label_col='Species_Label',
            title=f'Bin Abundance by Species - {treatment}',
            filename=f'abundance_species_{treatment}.png',
            include_unmapped=False
        )

# Plot 4: Per-treatment stacked barplots with species labels (with unmapped)
for treatment in treatments:
    treatment_data = df_merged[df_merged['Treatment'] == treatment]
    if not treatment_data.empty:
        create_stacked_barplot(
            treatment_data,
            group_by='Sample',
            label_col='Species_Label',
            title=f'Bin Abundance by Species (with unmapped) - {treatment}',
            filename=f'abundance_species_unmapped_{treatment}.png',
            include_unmapped=True
        )

# Plot 5: Faceted plot with all treatments - Genus
print("\nCreating faceted plot (Genus)...")
fig, axes = plt.subplots(1, len(treatments), figsize=(7*len(treatments), 6), sharey=True)

if len(treatments) == 1:
    axes = [axes]

for idx, treatment in enumerate(treatments):
    treatment_data = df_merged[df_merged['Treatment'] == treatment]

    # Group and pivot
    grouped = treatment_data.groupby(['Sample', 'Genus_Label'])['Relative_Abundance'].sum().reset_index()
    pivot_data = grouped.pivot(index='Sample', columns='Genus_Label', values='Relative_Abundance')
    pivot_data = pivot_data.fillna(0)

    # Sort by abundance
    col_sums = pivot_data.sum(axis=0).sort_values(ascending=False)
    pivot_data = pivot_data[col_sums.index]

    # Plot
    colors = sns.color_palette("husl", len(pivot_data.columns))
    pivot_data.plot(kind='bar', stacked=True, ax=axes[idx], color=colors, width=0.8, legend=False)

    axes[idx].set_title(treatment, fontsize=12, fontweight='bold')
    axes[idx].set_xlabel('Sample', fontsize=10)
    axes[idx].set_ylim(0, 100)
    axes[idx].grid(axis='y', alpha=0.3, linestyle='--')
    plt.setp(axes[idx].xaxis.get_majorticklabels(), rotation=45, ha='right')

axes[0].set_ylabel('Relative Abundance (%)', fontsize=10, fontweight='bold')

# Create shared legend
handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5),
          frameon=True, fontsize=9, title='Genus')

plt.suptitle('Bin Abundance by Genus - All Treatments', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{plots_dir}/abundance_genus_all_treatments_faceted.png", dpi=300, bbox_inches='tight')
print(f"  Saved: {plots_dir}/abundance_genus_all_treatments_faceted.png")
plt.close()

# Plot 6: Faceted plot with all treatments - Species
print("\nCreating faceted plot (Species)...")
fig, axes = plt.subplots(1, len(treatments), figsize=(7*len(treatments), 6), sharey=True)

if len(treatments) == 1:
    axes = [axes]

for idx, treatment in enumerate(treatments):
    treatment_data = df_merged[df_merged['Treatment'] == treatment]

    # Group and pivot
    grouped = treatment_data.groupby(['Sample', 'Species_Label'])['Relative_Abundance'].sum().reset_index()
    pivot_data = grouped.pivot(index='Sample', columns='Species_Label', values='Relative_Abundance')
    pivot_data = pivot_data.fillna(0)

    # Sort by abundance
    col_sums = pivot_data.sum(axis=0).sort_values(ascending=False)
    pivot_data = pivot_data[col_sums.index]

    # Plot
    colors = sns.color_palette("husl", len(pivot_data.columns))
    pivot_data.plot(kind='bar', stacked=True, ax=axes[idx], color=colors, width=0.8, legend=False)

    axes[idx].set_title(treatment, fontsize=12, fontweight='bold')
    axes[idx].set_xlabel('Sample', fontsize=10)
    axes[idx].set_ylim(0, 100)
    axes[idx].grid(axis='y', alpha=0.3, linestyle='--')
    plt.setp(axes[idx].xaxis.get_majorticklabels(), rotation=45, ha='right')

axes[0].set_ylabel('Relative Abundance (%)', fontsize=10, fontweight='bold')

# Create shared legend
handles, labels = axes[-1].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5),
          frameon=True, fontsize=9, title='Species')

plt.suptitle('Bin Abundance by Species - All Treatments', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{plots_dir}/abundance_species_all_treatments_faceted.png", dpi=300, bbox_inches='tight')
print(f"  Saved: {plots_dir}/abundance_species_all_treatments_faceted.png")
plt.close()

print("\n" + "="*80)
print("Plot generation complete!")
print("="*80)
print(f"\nAll plots saved to: {plots_dir}")
print("\nGenerated plots:")
print(f"  - Per-treatment genus plots (no unmapped): {len(treatments)} files")
print(f"  - Per-treatment genus plots (with unmapped): {len(treatments)} files")
print(f"  - Per-treatment species plots (no unmapped): {len(treatments)} files")
print(f"  - Per-treatment species plots (with unmapped): {len(treatments)} files")
print(f"  - Faceted genus plot (all treatments): 1 file")
print(f"  - Faceted species plot (all treatments): 1 file")
print(f"\nTotal plots: {len(treatments) * 4 + 2}")
PYTHON_EOF

    # Run Python script
    log "  Running Python plotting script..."

    # Check for required Python packages
    if ! python3 -c "import pandas, matplotlib, seaborn" 2>/dev/null; then
        log "  Installing required Python packages..."
        pip3 install pandas matplotlib seaborn 2>&1 | tail -20
    fi

    local treatments=($treatments_str)
    python3 "$python_script" "$OUTPUT_DIR" "${treatments[@]}" 2>&1 | tee "${FINAL_REPORT_DIR}/plotting.log"

    local exit_code=${PIPESTATUS[0]}

    if [ $exit_code -eq 0 ]; then
        log "  Successfully created all plots"
        return 0
    else
        log "ERROR: Failed to create plots"
        return 1
    fi
}

# Function to create summary report
create_summary_report() {
    local treatments_str="$1"

    log "Creating summary report..."

    local summary_file="${FINAL_REPORT_DIR}/report_summary.txt"
    local treatments=($treatments_str)

    cat > "$summary_file" << EOF
Final Report Summary
====================

Date: $(date)
Pipeline output directory: $OUTPUT_DIR
Number of treatments: ${#treatments[@]}
Treatments: ${treatments[*]}

Generated Outputs:
------------------

1. Merged Data:
   ${FINAL_REPORT_DIR}/merged_abundance_taxonomy.tsv
   - Combined abundance and taxonomy data from all treatments
   - Includes relative abundance per sample and taxonomic classification

2. Plots Directory:
   ${FINAL_REPORT_DIR}/plots/

   Per-treatment plots (Genus level):
   $(for t in "${treatments[@]}"; do echo "   - abundance_genus_${t}.png (mapped bins only)"; done)
   $(for t in "${treatments[@]}"; do echo "   - abundance_genus_unmapped_${t}.png (with unmapped reads)"; done)

   Per-treatment plots (Species level):
   $(for t in "${treatments[@]}"; do echo "   - abundance_species_${t}.png (mapped bins only)"; done)
   $(for t in "${treatments[@]}"; do echo "   - abundance_species_unmapped_${t}.png (with unmapped reads)"; done)

   Comparative plots (all treatments):
   - abundance_genus_all_treatments_faceted.png
   - abundance_species_all_treatments_faceted.png

3. Processing Logs:
   ${FINAL_REPORT_DIR}/plotting.log

Description of Plots:
---------------------

1. Per-treatment plots (genus/species, with/without unmapped):
   - Stacked barplots showing bin abundance across samples
   - Each bar represents a sample
   - Colors represent different genera/species
   - Y-axis: Relative abundance (%)
   - Versions with "unmapped" show unclassified reads in grey

2. Faceted plots (all treatments):
   - Side-by-side comparison of all treatments
   - Shared legend for easy comparison
   - Each facet shows one treatment
   - Useful for cross-treatment comparisons

Data Sources:
-------------
$(for t in "${treatments[@]}"; do
    echo "
Treatment: $t
  - Abundance: ${OUTPUT_DIR}/bin_collection/${t}/coverm_abundance_consolidated.tsv"

    if [ -f "${OUTPUT_DIR}/bin_collection/${t}/gtdbtk/gtdbtk.bac120.summary.tsv" ]; then
        local bac_count=$(tail -n +2 "${OUTPUT_DIR}/bin_collection/${t}/gtdbtk/gtdbtk.bac120.summary.tsv" | wc -l)
        echo "  - Taxonomy (Bacteria): ${bac_count} genomes"
    fi

    if [ -f "${OUTPUT_DIR}/bin_collection/${t}/gtdbtk/gtdbtk.ar53.summary.tsv" ]; then
        local ar_count=$(tail -n +2 "${OUTPUT_DIR}/bin_collection/${t}/gtdbtk/gtdbtk.ar53.summary.tsv" | wc -l)
        echo "  - Taxonomy (Archaea): ${ar_count} genomes"
    fi
done)

EOF

    log "Summary report created: $summary_file"
}

# Main execution
main() {
    log "Starting final report generation..."

    # Step 1: Collect all treatment data
    local treatments=$(collect_all_treatment_data)
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to collect treatment data"
        return 1
    fi

    # Step 2: Create plots
    if ! create_final_report_plots "$treatments"; then
        log "ERROR: Failed to create plots"
        return 1
    fi

    # Step 3: Create summary report
    create_summary_report "$treatments"

    log "Final report generation completed successfully"
    log "  Output directory: $FINAL_REPORT_DIR"
    log "  Plots directory: ${FINAL_REPORT_DIR}/plots"
    log "  Summary: ${FINAL_REPORT_DIR}/report_summary.txt"

    return 0
}

# Run main function
if main; then
    log "====== Final Report Generation Complete ======"
else
    log "====== Final Report Generation Failed ======"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
