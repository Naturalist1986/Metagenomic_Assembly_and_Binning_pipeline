#!/usr/bin/env python3
"""
Create stacked barplots of bin abundance across samples for each treatment.
Excludes soil bins and unmapped reads.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Configuration
BASE_DIR = Path("/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap")
MATRIX_DIR = BASE_DIR / "abundance_matrices"
OUTPUT_DIR = BASE_DIR / "abundance_plots"

# Create output directory
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

def load_and_filter_matrix(matrix_file):
    """Load matrix and filter out Soil bins and unmapped"""
    df = pd.read_csv(matrix_file, sep='\t')
    
    # Filter out Soil bins and unmapped
    df_filtered = df[~df['Bin'].str.contains('Soil', case=False, na=False)]
    df_filtered = df_filtered[df_filtered['Bin'] != 'unmapped']
    
    return df_filtered

def create_stacked_barplot(df, treatment, output_dir):
    """Create stacked barplot for a treatment"""
    
    # Set Bin as index
    df_plot = df.set_index('Bin')
    
    # Transpose so samples are on x-axis
    df_plot = df_plot.T
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create stacked bar plot
    df_plot.plot(kind='bar', stacked=True, ax=ax, 
                 colormap='tab20', width=0.8, edgecolor='white', linewidth=0.5)
    
    # Customize plot
    ax.set_xlabel('Sample', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=14, fontweight='bold')
    ax.set_title(f'{treatment} - Bin Abundance Across Samples\n(Soil bins excluded)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Format y-axis as percentage
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
    
    # Customize legend
    ax.legend(title='Bins', bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=8, title_fontsize=10, frameon=True, 
              fancybox=True, shadow=True)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / f"{treatment}_abundance_barplot.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created barplot: {output_file}")
    
    return output_file

def create_simplified_barplot(df, treatment, output_dir):
    """Create simplified barplot with better bin labels"""
    
    # Simplify bin names (remove sample prefix for cleaner labels)
    df_simple = df.copy()
    df_simple['Bin_Short'] = df_simple['Bin'].apply(
        lambda x: x.split('.', 1)[1] if '.' in x else x
    )
    
    # Group by simplified bin name (sum duplicates)
    sample_cols = [col for col in df_simple.columns if col not in ['Bin', 'Bin_Short']]
    df_grouped = df_simple.groupby('Bin_Short')[sample_cols].sum()
    
    # Transpose for plotting
    df_plot = df_grouped.T
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create stacked bar plot
    df_plot.plot(kind='bar', stacked=True, ax=ax,
                 colormap='tab20', width=0.8, edgecolor='white', linewidth=0.5)
    
    # Customize plot
    ax.set_xlabel('Sample', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=14, fontweight='bold')
    ax.set_title(f'{treatment} - Bin Abundance Across Samples (Simplified)\n(Soil bins excluded)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Format y-axis
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
    
    # Customize legend
    ax.legend(title='Bins', bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=8, title_fontsize=10, frameon=True,
              fancybox=True, shadow=True, ncol=1)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / f"{treatment}_abundance_barplot_simplified.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created simplified barplot: {output_file}")
    
    return output_file

def create_summary_table(df, treatment, output_dir):
    """Create summary table of bin presence across samples"""
    
    # Set Bin as index
    df_summary = df.set_index('Bin')
    
    # Create presence/absence (>1% abundance = present)
    df_presence = (df_summary > 1.0).astype(int)
    
    # Add summary columns
    df_presence['Total_Samples'] = df_presence.sum(axis=1)
    df_presence['Mean_Abundance'] = df_summary.mean(axis=1).round(2)
    df_presence['Max_Abundance'] = df_summary.max(axis=1).round(2)
    
    # Sort by total samples present
    df_presence = df_presence.sort_values('Total_Samples', ascending=False)
    
    # Save table
    output_file = output_dir / f"{treatment}_bin_summary.tsv"
    df_presence.to_csv(output_file, sep='\t')
    
    print(f"  ✓ Created summary table: {output_file}")
    
    return output_file

def main():
    print("=" * 60)
    print("Creating Abundance Barplots")
    print("=" * 60)
    print()
    
    # Find all relative abundance matrices
    matrix_files = sorted(MATRIX_DIR.glob("*_relative_abundance_matrix.tsv"))
    
    if not matrix_files:
        print(f"ERROR: No abundance matrices found in {MATRIX_DIR}")
        sys.exit(1)
    
    print(f"Found {len(matrix_files)} treatments to plot")
    print()
    
    # Process each treatment
    for matrix_file in matrix_files:
        treatment = matrix_file.stem.replace("_relative_abundance_matrix", "")
        print(f"Processing {treatment}...")
        
        # Load and filter data
        df = load_and_filter_matrix(matrix_file)
        
        if len(df) == 0:
            print(f"  WARNING: No bins remaining after filtering for {treatment}")
            continue
        
        n_bins = len(df)
        n_samples = len(df.columns) - 1  # Subtract Bin column
        
        print(f"  Bins: {n_bins} (after excluding Soil bins)")
        print(f"  Samples: {n_samples}")
        
        # Create plots
        create_stacked_barplot(df, treatment, OUTPUT_DIR)
        create_simplified_barplot(df, treatment, OUTPUT_DIR)
        
        # Create summary table
        create_summary_table(df, treatment, OUTPUT_DIR)
        
        print()
    
    print("=" * 60)
    print("All plots created successfully!")
    print("=" * 60)
    print()
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print("Files created:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        print(f"  {f.name}")
    print()
    print("View plots:")
    print(f"  ls {OUTPUT_DIR}/*.png")
    print()

if __name__ == "__main__":
    main()
