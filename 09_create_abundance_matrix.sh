#!/bin/bash
# 09_create_abundance_matrix.sh - Create abundance matrix from treatment-level CoverM results

BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
COVERM_DIR="${BASE_DIR}/coverm_treatment_level"
OUTPUT_DIR="${BASE_DIR}/abundance_matrices"

echo "========================================="
echo "Creating Abundance Matrices"
echo "========================================="
echo ""

# Check if CoverM results exist
if [ ! -d "$COVERM_DIR" ]; then
    echo "ERROR: CoverM treatment-level directory not found: $COVERM_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Process each treatment
for treatment_dir in "$COVERM_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    
    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"
    
    # Find all abundance files for this treatment
    abundance_files=("${treatment_dir}"/*_abundance.tsv)
    
    if [ ${#abundance_files[@]} -eq 0 ]; then
        echo "  WARNING: No abundance files found"
        continue
    fi
    
    # Output files
    rel_abund_matrix="${OUTPUT_DIR}/${treatment}_relative_abundance_matrix.tsv"
    mean_cov_matrix="${OUTPUT_DIR}/${treatment}_mean_coverage_matrix.tsv"
    tpm_matrix="${OUTPUT_DIR}/${treatment}_TPM_matrix.tsv"
    
    echo "  Creating matrices..."
    
    # Use Python to create matrices
    python3 << 'PYTHON_SCRIPT'
import sys
import os
import pandas as pd
from pathlib import Path

treatment = sys.argv[1]
treatment_dir = sys.argv[2]
output_dir = sys.argv[3]

# Find all abundance files
abundance_files = sorted(Path(treatment_dir).glob("*_abundance.tsv"))

if not abundance_files:
    print(f"  No abundance files found for {treatment}")
    sys.exit(1)

# Dictionary to store dataframes
rel_abund_dfs = []
mean_cov_dfs = []
tpm_dfs = []

for abund_file in abundance_files:
    sample = abund_file.stem.replace("_abundance", "")
    
    try:
        df = pd.read_csv(abund_file, sep='\t')
        
        # Find the relevant columns
        genome_col = df.columns[0]
        rel_col = [c for c in df.columns if 'Relative Abundance' in c][0]
        mean_col = [c for c in df.columns if 'Mean' in c and 'Trimmed' not in c][0]
        tpm_col = [c for c in df.columns if 'TPM' in c][0]
        
        # Extract data
        rel_df = df[[genome_col, rel_col]].rename(columns={rel_col: sample})
        mean_df = df[[genome_col, mean_col]].rename(columns={mean_col: sample})
        tpm_df = df[[genome_col, tpm_col]].rename(columns={tpm_col: sample})
        
        rel_abund_dfs.append(rel_df)
        mean_cov_dfs.append(mean_df)
        tpm_dfs.append(tpm_df)
        
    except Exception as e:
        print(f"  WARNING: Error processing {abund_file}: {e}")
        continue

# Merge all dataframes
if rel_abund_dfs:
    rel_matrix = rel_abund_dfs[0]
    for df in rel_abund_dfs[1:]:
        rel_matrix = pd.merge(rel_matrix, df, on=genome_col, how='outer')
    
    mean_matrix = mean_cov_dfs[0]
    for df in mean_cov_dfs[1:]:
        mean_matrix = pd.merge(mean_matrix, df, on=genome_col, how='outer')
    
    tpm_matrix = tpm_dfs[0]
    for df in tpm_dfs[1:]:
        tpm_matrix = pd.merge(tpm_matrix, df, on=genome_col, how='outer')
    
    # Fill NaN with 0
    rel_matrix = rel_matrix.fillna(0)
    mean_matrix = mean_matrix.fillna(0)
    tpm_matrix = tpm_matrix.fillna(0)
    
    # Rename first column to "Bin"
    rel_matrix = rel_matrix.rename(columns={genome_col: 'Bin'})
    mean_matrix = mean_matrix.rename(columns={genome_col: 'Bin'})
    tpm_matrix = tpm_matrix.rename(columns={genome_col: 'Bin'})
    
    # Save matrices
    rel_file = f"{output_dir}/{treatment}_relative_abundance_matrix.tsv"
    mean_file = f"{output_dir}/{treatment}_mean_coverage_matrix.tsv"
    tpm_file = f"{output_dir}/{treatment}_TPM_matrix.tsv"
    
    rel_matrix.to_csv(rel_file, sep='\t', index=False)
    mean_matrix.to_csv(mean_file, sep='\t', index=False)
    tpm_matrix.to_csv(tpm_file, sep='\t', index=False)
    
    print(f"  ✓ Created matrices with {len(rel_matrix)} bins x {len(rel_matrix.columns)-1} samples")
    print(f"    - {rel_file}")
    print(f"    - {mean_file}")
    print(f"    - {tpm_file}")
else:
    print(f"  ERROR: No valid data found for {treatment}")
    sys.exit(1)
PYTHON_SCRIPT

    python3 -c "
import sys
sys.argv = ['', '$treatment', '$treatment_dir', '$OUTPUT_DIR']
exec(open('/dev/stdin').read())
" << 'PYTHON_SCRIPT'
import sys
import os
import pandas as pd
from pathlib import Path

treatment = sys.argv[1]
treatment_dir = sys.argv[2]
output_dir = sys.argv[3]

# Find all abundance files
abundance_files = sorted(Path(treatment_dir).glob("*_abundance.tsv"))

if not abundance_files:
    print(f"  No abundance files found for {treatment}")
    sys.exit(0)

# Dictionary to store dataframes
rel_abund_dfs = []
mean_cov_dfs = []
tpm_dfs = []

for abund_file in abundance_files:
    sample = abund_file.stem.replace("_abundance", "")
    
    try:
        df = pd.read_csv(abund_file, sep='\t')
        
        # Find the relevant columns
        genome_col = df.columns[0]
        rel_col = [c for c in df.columns if 'Relative Abundance' in c][0]
        mean_col = [c for c in df.columns if 'Mean' in c and 'Trimmed' not in c][0]
        tpm_col = [c for c in df.columns if 'TPM' in c][0]
        
        # Extract data
        rel_df = df[[genome_col, rel_col]].rename(columns={rel_col: sample})
        mean_df = df[[genome_col, mean_col]].rename(columns={mean_col: sample})
        tpm_df = df[[genome_col, tpm_col]].rename(columns={tpm_col: sample})
        
        rel_abund_dfs.append(rel_df)
        mean_cov_dfs.append(mean_df)
        tpm_dfs.append(tpm_df)
        
    except Exception as e:
        print(f"  WARNING: Error processing {abund_file}: {e}")
        continue

# Merge all dataframes
if rel_abund_dfs:
    rel_matrix = rel_abund_dfs[0]
    for df in rel_abund_dfs[1:]:
        rel_matrix = pd.merge(rel_matrix, df, on=genome_col, how='outer')
    
    mean_matrix = mean_cov_dfs[0]
    for df in mean_cov_dfs[1:]:
        mean_matrix = pd.merge(mean_matrix, df, on=genome_col, how='outer')
    
    tpm_matrix = tpm_dfs[0]
    for df in tpm_dfs[1:]:
        tpm_matrix = pd.merge(tpm_matrix, df, on=genome_col, how='outer')
    
    # Fill NaN with 0
    rel_matrix = rel_matrix.fillna(0)
    mean_matrix = mean_matrix.fillna(0)
    tpm_matrix = tpm_matrix.fillna(0)
    
    # Rename first column to "Bin"
    rel_matrix = rel_matrix.rename(columns={genome_col: 'Bin'})
    mean_matrix = mean_matrix.rename(columns={genome_col: 'Bin'})
    tpm_matrix = tpm_matrix.rename(columns={genome_col: 'Bin'})
    
    # Save matrices
    rel_file = f"{output_dir}/{treatment}_relative_abundance_matrix.tsv"
    mean_file = f"{output_dir}/{treatment}_mean_coverage_matrix.tsv"
    tpm_file = f"{output_dir}/{treatment}_TPM_matrix.tsv"
    
    rel_matrix.to_csv(rel_file, sep='\t', index=False)
    mean_matrix.to_csv(mean_file, sep='\t', index=False)
    tpm_matrix.to_csv(tpm_file, sep='\t', index=False)
    
    print(f"  ✓ Created matrices: {len(rel_matrix)} bins x {len(rel_matrix.columns)-1} samples")
else:
    print(f"  ERROR: No valid data for {treatment}")
PYTHON_SCRIPT
    
    echo ""
done

echo "========================================="
echo "Abundance matrices created in:"
echo "  $OUTPUT_DIR"
echo ""
echo "Files created:"
ls -lh "$OUTPUT_DIR"/*.tsv 2>/dev/null | awk '{print "  " $9}'
echo ""
echo "Example usage:"
echo "  # View relative abundance matrix for RH"
echo "  less -S ${OUTPUT_DIR}/RH_relative_abundance_matrix.tsv"
echo ""
echo "  # Load in R for analysis"
echo "  R"
echo "  > library(tidyverse)"
echo "  > rel_abund <- read_tsv('${OUTPUT_DIR}/RH_relative_abundance_matrix.tsv')"
echo "  > head(rel_abund)"
echo "========================================="
