#!/usr/bin/env python3
"""
Create stacked barplots of bin abundance across samples for each treatment.
Includes taxonomic annotations from GTDB-Tk.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import re

# Configuration
BASE_DIR = Path("/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap")
MATRIX_DIR = BASE_DIR / "abundance_matrices"
GTDBTK_DIR = BASE_DIR / "gtdbtk_results"
OUTPUT_DIR = BASE_DIR / "abundance_plots"

# Create output directory
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 9

def load_taxonomy():
    """Load and parse GTDB-Tk taxonomy results"""
    taxonomy = {}
    
    # Load bacterial taxonomy
    bac_file = GTDBTK_DIR / "gtdbtk.bac120.summary.tsv"
    if bac_file.exists():
        df_bac = pd.read_csv(bac_file, sep='\t')
        for _, row in df_bac.iterrows():
            bin_name = row['user_genome']
            classification = row['classification']
            taxonomy[bin_name] = parse_taxonomy(classification)
    
    # Load archaeal taxonomy
    arc_file = GTDBTK_DIR / "gtdbtk.ar53.summary.tsv"
    if arc_file.exists():
        df_arc = pd.read_csv(arc_file, sep='\t')
        for _, row in df_arc.iterrows():
            bin_name = row['user_genome']
            classification = row['classification']
            taxonomy[bin_name] = parse_taxonomy(classification)
    
    return taxonomy

def parse_taxonomy(classification):
    """Parse GTDB taxonomy string into components"""
    tax_dict = {
        'domain': 'Unknown',
        'phylum': 'Unknown',
        'class': 'Unknown',
        'order': 'Unknown',
        'family': 'Unknown',
        'genus': 'Unknown',
        'species': 'Unknown'
    }
    
    if pd.isna(classification):
        return tax_dict
    
    # Parse taxonomy string
    parts = classification.split(';')
    for part in parts:
        part = part.strip()
        if part.startswith('d__'):
            tax_dict['domain'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('p__'):
            tax_dict['phylum'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('c__'):
            tax_dict['class'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('o__'):
            tax_dict['order'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('f__'):
            tax_dict['family'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('g__'):
            tax_dict['genus'] = part[3:] if part[3:] else 'Unknown'
        elif part.startswith('s__'):
            tax_dict['species'] = part[3:] if part[3:] else 'Unknown'
    
    return tax_dict

def create_bin_label(bin_name, taxonomy):
    """Create informative bin label with taxonomy"""
    if bin_name not in taxonomy:
        return bin_name
    
    tax = taxonomy[bin_name]
    
    # Try to create label from most specific to least specific
    if tax['species'] != 'Unknown' and not tax['species'].startswith('sp'):
        # Clean species name (remove strain info)
        species = tax['species'].split()[0] if ' ' in tax['species'] else tax['species']
        return f"{tax['genus']} {species}"
    elif tax['genus'] != 'Unknown':
        return tax['genus']
    elif tax['family'] != 'Unknown':
        return f"{tax['family']} (family)"
    elif tax['order'] != 'Unknown':
        return f"{tax['order']} (order)"
    elif tax['class'] != 'Unknown':
        return f"{tax['class']} (class)"
    elif tax['phylum'] != 'Unknown':
        return tax['phylum']
    else:
        return bin_name

def load_and_filter_matrix(matrix_file, taxonomy):
    """Load matrix, filter, and add taxonomy labels"""
    df = pd.read_csv(matrix_file, sep='\t')
    
    # Filter out Soil bins and unmapped
    df_filtered = df[~df['Bin'].str.contains('Soil', case=False, na=False)]
    df_filtered = df_filtered[df_filtered['Bin'] != 'unmapped']
    
    # Add taxonomy labels
    df_filtered['Bin_Label'] = df_filtered['Bin'].apply(
        lambda x: create_bin_label(x, taxonomy)
    )
    
    # Add phylum for coloring
    df_filtered['Phylum'] = df_filtered['Bin'].apply(
        lambda x: taxonomy.get(x, {}).get('phylum', 'Unknown')
    )
    
    return df_filtered

def get_phylum_colors(phyla):
    """Generate consistent colors for phyla"""
    unique_phyla = sorted(set(phyla))
    
    # Use a colormap with enough distinct colors
    if len(unique_phyla) <= 20:
        cmap = plt.cm.tab20
    else:
        cmap = plt.cm.hsv
    
    colors = {}
    for i, phylum in enumerate(unique_phyla):
        colors[phylum] = cmap(i / len(unique_phyla))
    
    return colors

def create_annotated_barplot(df, treatment, output_dir, taxonomy):
    """Create stacked barplot with taxonomy annotations"""
    
    # Prepare data
    sample_cols = [col for col in df.columns if col not in ['Bin', 'Bin_Label', 'Phylum']]
    
    # Create plotting dataframe with bin labels
    df_plot = df.set_index('Bin_Label')[sample_cols].T
    
    # Get phylum colors
    phylum_map = dict(zip(df['Bin_Label'], df['Phylum']))
    phyla = [phylum_map[col] for col in df_plot.columns]
    phylum_colors = get_phylum_colors(phyla)
    bin_colors = [phylum_colors[phylum_map[col]] for col in df_plot.columns]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Create stacked bar plot
    df_plot.plot(kind='bar', stacked=True, ax=ax, 
                 color=bin_colors, width=0.8, edgecolor='white', linewidth=0.5)
    
    # Customize plot
    ax.set_xlabel('Sample', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=14, fontweight='bold')
    ax.set_title(f'{treatment} - Bin Abundance Across Samples\n(with taxonomic annotations, Soil bins excluded)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Format y-axis
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
    
    # Create custom legend grouped by phylum
    handles, labels = ax.get_legend_handles_labels()
    
    # Group by phylum
    phylum_groups = {}
    for handle, label in zip(handles, labels):
        phylum = phylum_map.get(label, 'Unknown')
        if phylum not in phylum_groups:
            phylum_groups[phylum] = []
        phylum_groups[phylum].append((handle, label))
    
    # Create legend with phylum grouping
    legend_handles = []
    legend_labels = []
    
    for phylum in sorted(phylum_groups.keys()):
        # Add phylum header (bold)
        legend_labels.append(f"--- {phylum} ---")
        legend_handles.append(plt.Rectangle((0,0),1,1, fc="white", ec="white", alpha=0))
        
        # Add bins in this phylum
        for handle, label in phylum_groups[phylum]:
            legend_handles.append(handle)
            legend_labels.append(f"  {label}")
    
    ax.legend(legend_handles, legend_labels, 
              bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=8, frameon=True, fancybox=True, shadow=True)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / f"{treatment}_abundance_barplot_annotated.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created annotated barplot: {output_file}")
    
    return output_file

def create_phylum_barplot(df, treatment, output_dir):
    """Create phylum-level stacked barplot"""
    
    # Get sample columns
    sample_cols = [col for col in df.columns if col not in ['Bin', 'Bin_Label', 'Phylum']]
    
    # Group by phylum and sum abundances
    phylum_df = df.groupby('Phylum')[sample_cols].sum().T
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Create stacked bar plot
    phylum_df.plot(kind='bar', stacked=True, ax=ax,
                   colormap='tab20', width=0.8, edgecolor='white', linewidth=0.5)
    
    # Customize plot
    ax.set_xlabel('Sample', fontsize=14, fontweight='bold')
    ax.set_ylabel('Relative Abundance (%)', fontsize=14, fontweight='bold')
    ax.set_title(f'{treatment} - Phylum-Level Abundance\n(Soil bins excluded)', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Format y-axis
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
    
    # Customize legend
    ax.legend(title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=10, title_fontsize=11, frameon=True, 
              fancybox=True, shadow=True)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # Tight layout
    plt.tight_layout()
    
    # Save plot
    output_file = output_dir / f"{treatment}_phylum_barplot.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  ✓ Created phylum barplot: {output_file}")
    
    return output_file

def create_taxonomy_table(df, treatment, output_dir):
    """Create summary table with taxonomy"""
    
    # Create summary
    sample_cols = [col for col in df.columns if col not in ['Bin', 'Bin_Label', 'Phylum']]
    
    summary = pd.DataFrame({
        'Bin': df['Bin'],
        'Taxonomy': df['Bin_Label'],
        'Phylum': df['Phylum'],
        'Mean_Abundance': df[sample_cols].mean(axis=1).round(2),
        'Max_Abundance': df[sample_cols].max(axis=1).round(2),
        'Samples_Present': (df[sample_cols] > 1.0).sum(axis=1)
    })
    
    # Sort by mean abundance
    summary = summary.sort_values('Mean_Abundance', ascending=False)
    
    # Save
    output_file = output_dir / f"{treatment}_taxonomy_summary.tsv"
    summary.to_csv(output_file, sep='\t', index=False)
    
    print(f"  ✓ Created taxonomy summary: {output_file}")
    
    return output_file

def main():
    print("=" * 70)
    print("Creating Taxonomically Annotated Abundance Barplots")
    print("=" * 70)
    print()
    
    # Load taxonomy
    print("Loading GTDB-Tk taxonomy...")
    taxonomy = load_taxonomy()
    
    if not taxonomy:
        print("WARNING: No GTDB-Tk results found. Using bin names only.")
        print(f"Expected files in: {GTDBTK_DIR}")
        print()
    else:
        print(f"  Loaded taxonomy for {len(taxonomy)} bins")
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
        df = load_and_filter_matrix(matrix_file, taxonomy)
        
        if len(df) == 0:
            print(f"  WARNING: No bins remaining after filtering for {treatment}")
            continue
        
        n_bins = len(df)
        n_samples = len(df.columns) - 3  # Subtract Bin, Bin_Label, Phylum columns
        n_phyla = df['Phylum'].nunique()
        
        print(f"  Bins: {n_bins} (after excluding Soil bins)")
        print(f"  Samples: {n_samples}")
        print(f"  Phyla: {n_phyla}")
        
        # Create plots
        create_annotated_barplot(df, treatment, OUTPUT_DIR, taxonomy)
        create_phylum_barplot(df, treatment, OUTPUT_DIR)
        
        # Create taxonomy summary
        create_taxonomy_table(df, treatment, OUTPUT_DIR)
        
        print()
    
    print("=" * 70)
    print("All plots created successfully!")
    print("=" * 70)
    print()
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print("Files created:")
    for f in sorted(OUTPUT_DIR.iterdir()):
        if f.suffix in ['.png', '.tsv']:
            print(f"  {f.name}")
    print()
    print("View plots:")
    print(f"  ls {OUTPUT_DIR}/*_annotated.png")
    print(f"  ls {OUTPUT_DIR}/*_phylum_barplot.png")
    print()

if __name__ == "__main__":
    main()