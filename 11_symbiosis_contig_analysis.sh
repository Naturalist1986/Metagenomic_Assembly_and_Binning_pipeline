#!/usr/bin/env bash
#SBATCH --job-name=symbiosis_analysis
#SBATCH --array=0-99%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --account=ofinkel

# 11_symbiosis_contig_analysis.sh - Nitrogen Fixation and Nodulation Contig Analysis
#
# This script performs comprehensive analysis of symbiosis-related contigs:
# 1. Runs Prodigal to predict proteins from coassembly contigs
# 2. Runs InterProScan to annotate predicted proteins
# 3. Runs CoverM contig to calculate trimmed_mean coverage per sample
# 4. Filters contigs with nitrogen fixation/nodulation annotations
# 5. Calculates length-weighted mean coverage for each bin
# 6. Computes Pearson and Spearman correlations between symbiosis contigs and bins
#
# Currently supports coassembly mode only.

set -euo pipefail

# Source configuration and utilities
if [ -n "${PIPELINE_SCRIPT_DIR:-}" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# InterProScan path
INTERPROSCAN_PATH="${INTERPROSCAN_PATH:-/sci/backup/ofinkel/moshea/Tools/interproscan-5.76-107.0/interproscan.sh}"

# Set up temporary directory
TEMP_DIR=$(setup_temp_dir "symbiosis")

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################

# Run Prodigal to predict proteins from contigs
run_prodigal() {
    local contigs_file="$1"
    local output_dir="$2"
    local treatment="$3"

    log "Running Prodigal protein prediction for treatment $treatment..."

    local proteins_file="${output_dir}/predicted_proteins.faa"
    local genes_file="${output_dir}/predicted_genes.fna"
    local gff_file="${output_dir}/predicted_genes.gff"

    # Check if already completed
    if [ -f "$proteins_file" ] && [ -s "$proteins_file" ]; then
        log "  Prodigal output already exists, skipping..."
        echo "$proteins_file"
        return 0
    fi

    mkdir -p "$output_dir"

    # Activate prodigal environment
    activate_env metawrap-env

    if ! command -v prodigal &> /dev/null; then
        log "ERROR: Prodigal not found in environment"
        conda deactivate
        return 1
    fi

    log "  Input contigs: $contigs_file"
    log "  Output proteins: $proteins_file"

    prodigal \
        -i "$contigs_file" \
        -a "$proteins_file" \
        -d "$genes_file" \
        -f gff \
        -o "$gff_file" \
        -p meta \
        2>&1 | tee "${LOG_DIR}/${treatment}_prodigal.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    if [ $exit_code -eq 0 ] && [ -f "$proteins_file" ] && [ -s "$proteins_file" ]; then
        local protein_count=$(grep -c "^>" "$proteins_file" || echo 0)
        log "  Prodigal completed: $protein_count proteins predicted"
        echo "$proteins_file"
        return 0
    else
        log "ERROR: Prodigal failed with exit code $exit_code"
        return 1
    fi
}

# Run InterProScan on predicted proteins
run_interproscan() {
    local proteins_file="$1"
    local output_dir="$2"
    local treatment="$3"

    log "Running InterProScan for treatment $treatment..."

    local ipr_output_base="${output_dir}/interproscan_results"
    local ipr_tsv="${ipr_output_base}.tsv"

    # Check if already completed
    if [ -f "$ipr_tsv" ] && [ -s "$ipr_tsv" ]; then
        log "  InterProScan output already exists, skipping..."
        echo "$ipr_tsv"
        return 0
    fi

    mkdir -p "$output_dir"

    # Check InterProScan availability
    if [ ! -f "$INTERPROSCAN_PATH" ]; then
        log "ERROR: InterProScan not found at $INTERPROSCAN_PATH"
        return 1
    fi

    # Clean proteins file (remove asterisks from stop codons)
    local clean_proteins="${TEMP_DIR}/proteins_cleaned.faa"
    log "  Cleaning protein sequences (removing stop codon asterisks)..."
    sed 's/\*//g' "$proteins_file" > "$clean_proteins"

    # Create temp directory for InterProScan
    local ipr_tmp="${TEMP_DIR}/ipr_tmp"
    mkdir -p "$ipr_tmp"

    log "  Running InterProScan (this may take several hours)..."
    log "  Input: $clean_proteins"
    log "  Output base: $ipr_output_base"

    "$INTERPROSCAN_PATH" \
        -i "$clean_proteins" \
        -b "$ipr_output_base" \
        -T "$ipr_tmp" \
        -cpu ${SLURM_CPUS_PER_TASK:-32} \
        -dp \
        --goterms \
        --pa \
        2>&1 | tee "${LOG_DIR}/${treatment}_interproscan.log"

    local exit_code=${PIPESTATUS[0]}

    # Clean up temp files
    rm -rf "$ipr_tmp" "$clean_proteins"

    if [ $exit_code -eq 0 ] && [ -f "$ipr_tsv" ]; then
        local annotation_count=$(wc -l < "$ipr_tsv")
        log "  InterProScan completed: $annotation_count annotations"
        echo "$ipr_tsv"
        return 0
    else
        log "ERROR: InterProScan failed with exit code $exit_code"
        return 1
    fi
}

# Run CoverM contig to get per-contig coverage
run_coverm_contig() {
    local contigs_file="$1"
    local bam_dir="$2"
    local output_dir="$3"
    local treatment="$4"

    log "Running CoverM contig coverage for treatment $treatment..."

    local coverage_file="${output_dir}/contig_coverage.tsv"

    # Check if already completed
    if [ -f "$coverage_file" ] && [ -s "$coverage_file" ]; then
        log "  CoverM contig output already exists, skipping..."
        echo "$coverage_file"
        return 0
    fi

    mkdir -p "$output_dir"

    # Find all BAM files
    local bam_files=($(find "$bam_dir" -name "*.sorted.bam" -type f 2>/dev/null))

    if [ ${#bam_files[@]} -eq 0 ]; then
        log "ERROR: No BAM files found in $bam_dir"
        return 1
    fi

    log "  Found ${#bam_files[@]} BAM files"
    log "  Output: $coverage_file"

    # Activate CoverM environment
    activate_env checkm

    if ! command -v coverm &> /dev/null; then
        log "ERROR: CoverM not found in environment"
        conda deactivate
        return 1
    fi

    # Run CoverM contig with trimmed_mean method
    coverm contig \
        --bam-files "${bam_files[@]}" \
        --output-file "$coverage_file" \
        --methods trimmed_mean length \
        --min-read-aligned-percent 0.75 \
        --min-read-percent-identity 0.95 \
        --threads ${SLURM_CPUS_PER_TASK:-32} \
        2>&1 | tee "${LOG_DIR}/${treatment}_coverm_contig.log"

    local exit_code=${PIPESTATUS[0]}
    conda deactivate

    if [ $exit_code -eq 0 ] && [ -f "$coverage_file" ]; then
        local contig_count=$(tail -n +2 "$coverage_file" | wc -l)
        log "  CoverM contig completed: $contig_count contigs"
        echo "$coverage_file"
        return 0
    else
        log "ERROR: CoverM contig failed with exit code $exit_code"
        return 1
    fi
}

# Filter contigs with nitrogen fixation/nodulation annotations
filter_symbiosis_contigs() {
    local interproscan_tsv="$1"
    local contigs_file="$2"
    local output_dir="$3"
    local treatment="$4"

    log "Filtering nitrogen fixation/nodulation contigs for treatment $treatment..."

    local symbiosis_contigs_list="${output_dir}/symbiosis_contig_ids.txt"
    local symbiosis_annotations="${output_dir}/symbiosis_annotations.tsv"
    local symbiosis_fasta="${output_dir}/symbiosis_contigs.fasta"

    mkdir -p "$output_dir"

    # Create Python script for filtering
    local filter_script="${TEMP_DIR}/filter_symbiosis.py"

    cat > "$filter_script" << 'PYTHON_EOF'
#!/usr/bin/env python3
"""
Filter InterProScan results for nitrogen fixation and nodulation related annotations.
Search terms (case-insensitive):
- nif (includes nifH, nifD, nifK, nifE, nifN, nifB, etc.)
- nitrogenase
- nod-box
- nodulation
- nitrogen fixation
"""
import sys
import re
import pandas as pd
from collections import defaultdict

interproscan_file = sys.argv[1]
output_contigs = sys.argv[2]
output_annotations = sys.argv[3]

print(f"Reading InterProScan results: {interproscan_file}")

# Search patterns (case-insensitive)
# Note: We use word boundaries for 'nif' to avoid false positives
search_patterns = [
    r'\bnif[A-Z]?\b',           # nif, nifH, nifD, nifK, nifE, nifN, nifB, etc.
    r'\bnitrogenase\b',         # nitrogenase
    r'\bnod[-_]?box\b',         # nod-box, nodbox, nod_box
    r'\bnodulation\b',          # nodulation
    r'\bnitrogen\s*fixation\b', # nitrogen fixation
    r'\bfix[A-Z]?\b',           # fix genes (fixA, fixB, fixC, etc.) - nitrogen fixation related
]

# Combine into single pattern
combined_pattern = re.compile('|'.join(search_patterns), re.IGNORECASE)

# InterProScan TSV columns (standard format)
# 0: Protein ID
# 1: MD5 digest
# 2: Length
# 3: Analysis (e.g., Pfam, TIGRFAM)
# 4: Signature accession
# 5: Signature description
# 6: Start
# 7: Stop
# 8: E-value
# 9: Status
# 10: Date
# 11: InterPro accession
# 12: InterPro description
# 13: GO terms (optional)
# 14: Pathways (optional)

symbiosis_proteins = set()
symbiosis_contigs = set()
annotations = []

try:
    with open(interproscan_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 6:
                continue

            protein_id = fields[0]
            analysis = fields[3] if len(fields) > 3 else ''
            signature_desc = fields[5] if len(fields) > 5 else ''
            interpro_desc = fields[12] if len(fields) > 12 else ''
            go_terms = fields[13] if len(fields) > 13 else ''
            pathways = fields[14] if len(fields) > 14 else ''

            # Combine all annotation text for searching
            all_text = ' '.join([signature_desc, interpro_desc, go_terms, pathways])

            # Check if any symbiosis pattern matches
            match = combined_pattern.search(all_text)
            if match:
                symbiosis_proteins.add(protein_id)

                # Extract contig ID from protein ID
                # Prodigal format: contig_ID_protein_number (e.g., NODE_1_length_1000_cov_10_1)
                # The contig ID is everything except the last underscore and number
                parts = protein_id.rsplit('_', 1)
                if len(parts) == 2 and parts[1].isdigit():
                    contig_id = parts[0]
                else:
                    contig_id = protein_id

                symbiosis_contigs.add(contig_id)

                annotations.append({
                    'protein_id': protein_id,
                    'contig_id': contig_id,
                    'analysis': analysis,
                    'signature_desc': signature_desc,
                    'interpro_desc': interpro_desc,
                    'matched_term': match.group(),
                    'go_terms': go_terms,
                    'pathways': pathways
                })

except Exception as e:
    print(f"ERROR reading file: {e}")
    sys.exit(1)

print(f"Found {len(symbiosis_proteins)} symbiosis-related proteins")
print(f"Found {len(symbiosis_contigs)} unique contigs with symbiosis annotations")

# Write contig IDs
with open(output_contigs, 'w') as f:
    for contig_id in sorted(symbiosis_contigs):
        f.write(f"{contig_id}\n")

print(f"Contig IDs written to: {output_contigs}")

# Write detailed annotations
if annotations:
    df = pd.DataFrame(annotations)
    df.to_csv(output_annotations, sep='\t', index=False)
    print(f"Annotations written to: {output_annotations}")

    # Print summary of matched terms
    term_counts = df['matched_term'].str.lower().value_counts()
    print("\nMatched term summary:")
    for term, count in term_counts.items():
        print(f"  {term}: {count}")
else:
    # Create empty annotations file with header
    with open(output_annotations, 'w') as f:
        f.write("protein_id\tcontig_id\tanalysis\tsignature_desc\tinterpro_desc\tmatched_term\tgo_terms\tpathways\n")
    print("No symbiosis annotations found")

print("\nDone!")
PYTHON_EOF

    # Run the filter script
    python3 "$filter_script" "$interproscan_tsv" "$symbiosis_contigs_list" "$symbiosis_annotations" \
        2>&1 | tee "${LOG_DIR}/${treatment}_filter_symbiosis.log"

    local exit_code=${PIPESTATUS[0]}

    if [ $exit_code -ne 0 ]; then
        log "ERROR: Symbiosis filtering failed"
        return 1
    fi

    # Extract symbiosis contigs to FASTA
    if [ -f "$symbiosis_contigs_list" ] && [ -s "$symbiosis_contigs_list" ]; then
        local contig_count=$(wc -l < "$symbiosis_contigs_list")
        log "  Found $contig_count symbiosis-related contigs"

        # Use seqtk or Python to extract sequences
        if command -v seqtk &> /dev/null; then
            seqtk subseq "$contigs_file" "$symbiosis_contigs_list" > "$symbiosis_fasta"
        else
            # Fallback: use Python
            python3 << EXTRACT_EOF
import sys
contigs_to_extract = set()
with open("$symbiosis_contigs_list", 'r') as f:
    for line in f:
        contigs_to_extract.add(line.strip())

from Bio import SeqIO
with open("$symbiosis_fasta", 'w') as out_f:
    for record in SeqIO.parse("$contigs_file", "fasta"):
        # Match by exact ID or first field
        record_id = record.id.split()[0]
        if record_id in contigs_to_extract:
            SeqIO.write(record, out_f, "fasta")
EXTRACT_EOF
        fi

        if [ -f "$symbiosis_fasta" ]; then
            local extracted=$(grep -c "^>" "$symbiosis_fasta" || echo 0)
            log "  Extracted $extracted contigs to $symbiosis_fasta"
        fi
    else
        log "  No symbiosis-related contigs found"
        touch "$symbiosis_fasta"
    fi

    echo "$symbiosis_contigs_list|$symbiosis_annotations|$symbiosis_fasta"
    return 0
}

# Calculate bin coverage and run correlation analysis
run_correlation_analysis() {
    local coverage_file="$1"
    local symbiosis_contigs_list="$2"
    local bins_dir="$3"
    local output_dir="$4"
    local treatment="$5"

    log "Running correlation analysis for treatment $treatment..."

    local correlation_output="${output_dir}/symbiosis_bin_correlations.tsv"
    local summary_output="${output_dir}/correlation_summary.tsv"

    mkdir -p "$output_dir"

    # Create comprehensive Python analysis script
    local analysis_script="${TEMP_DIR}/correlation_analysis.py"

    cat > "$analysis_script" << 'PYTHON_EOF'
#!/usr/bin/env python3
"""
Correlation analysis between symbiosis contigs and bin coverage profiles.

This script:
1. Loads contig coverage data (trimmed_mean from CoverM)
2. Maps contigs to bins
3. Calculates length-weighted mean coverage for each bin (excluding symbiosis contigs)
4. Computes Pearson and Spearman correlations between each symbiosis contig and each bin
5. Outputs comprehensive correlation results
"""
import sys
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict

# Input arguments
coverage_file = sys.argv[1]
symbiosis_contigs_file = sys.argv[2]
bins_dir = sys.argv[3]
output_dir = sys.argv[4]
treatment = sys.argv[5]

print(f"=== Correlation Analysis for Treatment: {treatment} ===")
print(f"Coverage file: {coverage_file}")
print(f"Symbiosis contigs file: {symbiosis_contigs_file}")
print(f"Bins directory: {bins_dir}")
print(f"Output directory: {output_dir}")

# ============================================================================
# Step 1: Load coverage data
# ============================================================================
print("\n[Step 1] Loading coverage data...")

try:
    coverage_df = pd.read_csv(coverage_file, sep='\t')
    print(f"  Loaded coverage for {len(coverage_df)} contigs")
    print(f"  Columns: {list(coverage_df.columns)}")
except Exception as e:
    print(f"ERROR loading coverage file: {e}")
    sys.exit(1)

# Identify contig ID column and coverage columns
contig_col = coverage_df.columns[0]  # Usually 'Contig' or first column

# Find trimmed_mean columns (sample coverage) and length column
coverage_cols = []
length_col = None

for col in coverage_df.columns:
    if 'Trimmed Mean' in col or 'trimmed_mean' in col.lower():
        coverage_cols.append(col)
    elif col.lower() == 'length' or 'Length' in col:
        length_col = col

if not coverage_cols:
    print("ERROR: No trimmed_mean coverage columns found")
    print(f"  Available columns: {list(coverage_df.columns)}")
    sys.exit(1)

print(f"  Found {len(coverage_cols)} coverage columns (samples)")
print(f"  Length column: {length_col}")

# Extract sample names from column headers
# CoverM format: "sample.sorted Trimmed Mean" -> "sample"
sample_names = []
for col in coverage_cols:
    sample = col.replace(' Trimmed Mean', '').replace('.sorted', '')
    sample_names.append(sample)

print(f"  Samples: {sample_names}")

# Create clean coverage matrix
coverage_matrix = coverage_df[[contig_col] + coverage_cols].copy()
coverage_matrix.columns = ['contig_id'] + sample_names

# Add length information
if length_col:
    coverage_matrix['length'] = coverage_df[length_col]
else:
    print("WARNING: No length column found, using default length=1000 for all contigs")
    coverage_matrix['length'] = 1000

# Clean contig IDs (remove description after first space if present)
coverage_matrix['contig_id'] = coverage_matrix['contig_id'].str.split().str[0]

print(f"  Coverage matrix shape: {coverage_matrix.shape}")

# ============================================================================
# Step 2: Load symbiosis contig IDs
# ============================================================================
print("\n[Step 2] Loading symbiosis contig IDs...")

symbiosis_contigs = set()
if os.path.exists(symbiosis_contigs_file) and os.path.getsize(symbiosis_contigs_file) > 0:
    with open(symbiosis_contigs_file, 'r') as f:
        for line in f:
            contig_id = line.strip()
            if contig_id:
                symbiosis_contigs.add(contig_id)

print(f"  Loaded {len(symbiosis_contigs)} symbiosis contig IDs")

if len(symbiosis_contigs) == 0:
    print("WARNING: No symbiosis contigs found. Creating empty output files.")
    # Create empty output files
    pd.DataFrame(columns=['treatment', 'symbiosis_contig', 'bin', 'pearson_r', 'pearson_p',
                          'spearman_r', 'spearman_p']).to_csv(
        os.path.join(output_dir, 'symbiosis_bin_correlations.tsv'), sep='\t', index=False)
    pd.DataFrame(columns=['treatment', 'symbiosis_contig', 'best_matching_bin', 'best_pearson_r',
                          'best_spearman_r']).to_csv(
        os.path.join(output_dir, 'correlation_summary.tsv'), sep='\t', index=False)
    print("Empty output files created.")
    sys.exit(0)

# ============================================================================
# Step 3: Load bin-to-contig mappings
# ============================================================================
print("\n[Step 3] Loading bin-to-contig mappings...")

bin_to_contigs = defaultdict(set)
contig_to_bin = {}

# Find all bin files (.fa) in the bins directory
bin_files = glob.glob(os.path.join(bins_dir, '*.fa'))

if not bin_files:
    print(f"ERROR: No bin files found in {bins_dir}")
    sys.exit(1)

print(f"  Found {len(bin_files)} bin files")

for bin_file in bin_files:
    bin_name = os.path.basename(bin_file).replace('.fa', '')

    with open(bin_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                contig_id = line[1:].strip().split()[0]
                bin_to_contigs[bin_name].add(contig_id)
                contig_to_bin[contig_id] = bin_name

print(f"  Mapped {len(contig_to_bin)} contigs to {len(bin_to_contigs)} bins")

# ============================================================================
# Step 4: Calculate bin coverage profiles (length-weighted mean)
# ============================================================================
print("\n[Step 4] Calculating bin coverage profiles (length-weighted mean)...")

# Exclude symbiosis contigs from bin profile calculation to avoid circularity
bin_coverage_profiles = {}

for bin_name, bin_contigs in bin_to_contigs.items():
    # Get contigs for this bin, excluding symbiosis contigs
    core_contigs = bin_contigs - symbiosis_contigs

    if len(core_contigs) == 0:
        print(f"  WARNING: Bin {bin_name} has no core contigs after excluding symbiosis contigs")
        continue

    # Get coverage data for core contigs
    bin_coverage_data = coverage_matrix[coverage_matrix['contig_id'].isin(core_contigs)]

    if len(bin_coverage_data) == 0:
        print(f"  WARNING: No coverage data found for bin {bin_name}")
        continue

    # Calculate length-weighted mean coverage for each sample
    total_length = bin_coverage_data['length'].sum()

    if total_length == 0:
        print(f"  WARNING: Bin {bin_name} has zero total length")
        continue

    weighted_coverage = {}
    for sample in sample_names:
        weighted_sum = (bin_coverage_data[sample] * bin_coverage_data['length']).sum()
        weighted_coverage[sample] = weighted_sum / total_length

    bin_coverage_profiles[bin_name] = weighted_coverage

print(f"  Calculated coverage profiles for {len(bin_coverage_profiles)} bins")

if len(bin_coverage_profiles) == 0:
    print("ERROR: No valid bin coverage profiles could be calculated")
    sys.exit(1)

# ============================================================================
# Step 5: Get symbiosis contig coverage profiles
# ============================================================================
print("\n[Step 5] Getting symbiosis contig coverage profiles...")

symbiosis_coverage = coverage_matrix[coverage_matrix['contig_id'].isin(symbiosis_contigs)]

print(f"  Found coverage data for {len(symbiosis_coverage)} symbiosis contigs")

# Check which symbiosis contigs are missing from coverage data
missing_contigs = symbiosis_contigs - set(symbiosis_coverage['contig_id'])
if missing_contigs:
    print(f"  WARNING: {len(missing_contigs)} symbiosis contigs not found in coverage data")

# ============================================================================
# Step 6: Compute correlations
# ============================================================================
print("\n[Step 6] Computing Pearson and Spearman correlations...")

correlation_results = []

for _, contig_row in symbiosis_coverage.iterrows():
    contig_id = contig_row['contig_id']
    contig_profile = np.array([contig_row[s] for s in sample_names])

    # Skip if all values are zero or NaN
    if np.all(contig_profile == 0) or np.all(np.isnan(contig_profile)):
        continue

    for bin_name, bin_profile_dict in bin_coverage_profiles.items():
        bin_profile = np.array([bin_profile_dict[s] for s in sample_names])

        # Skip if bin profile is all zeros
        if np.all(bin_profile == 0):
            continue

        # Calculate Pearson correlation
        try:
            pearson_r, pearson_p = stats.pearsonr(contig_profile, bin_profile)
        except Exception:
            pearson_r, pearson_p = np.nan, np.nan

        # Calculate Spearman correlation
        try:
            spearman_r, spearman_p = stats.spearmanr(contig_profile, bin_profile)
        except Exception:
            spearman_r, spearman_p = np.nan, np.nan

        # Get contig length and assigned bin (if any)
        contig_length = contig_row['length']
        assigned_bin = contig_to_bin.get(contig_id, 'unbinned')

        correlation_results.append({
            'treatment': treatment,
            'symbiosis_contig': contig_id,
            'contig_length': contig_length,
            'assigned_bin': assigned_bin,
            'compared_bin': bin_name,
            'is_same_bin': assigned_bin == bin_name,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p
        })

print(f"  Computed {len(correlation_results)} correlations")

# ============================================================================
# Step 7: Save results
# ============================================================================
print("\n[Step 7] Saving results...")

# Create results dataframe
results_df = pd.DataFrame(correlation_results)

# Save full correlation results
correlation_output = os.path.join(output_dir, 'symbiosis_bin_correlations.tsv')
results_df.to_csv(correlation_output, sep='\t', index=False, float_format='%.6f')
print(f"  Full correlations saved to: {correlation_output}")

# Create summary: best matching bin for each symbiosis contig
if len(results_df) > 0:
    summary_data = []

    for contig_id in results_df['symbiosis_contig'].unique():
        contig_results = results_df[results_df['symbiosis_contig'] == contig_id]

        # Find best matching bin by Pearson correlation
        best_pearson_idx = contig_results['pearson_r'].idxmax()
        best_pearson_row = contig_results.loc[best_pearson_idx]

        # Find best matching bin by Spearman correlation
        best_spearman_idx = contig_results['spearman_r'].idxmax()
        best_spearman_row = contig_results.loc[best_spearman_idx]

        contig_length = contig_results['contig_length'].iloc[0]
        assigned_bin = contig_results['assigned_bin'].iloc[0]

        summary_data.append({
            'treatment': treatment,
            'symbiosis_contig': contig_id,
            'contig_length': contig_length,
            'assigned_bin': assigned_bin,
            'best_pearson_bin': best_pearson_row['compared_bin'],
            'best_pearson_r': best_pearson_row['pearson_r'],
            'best_pearson_p': best_pearson_row['pearson_p'],
            'pearson_matches_assigned': best_pearson_row['compared_bin'] == assigned_bin,
            'best_spearman_bin': best_spearman_row['compared_bin'],
            'best_spearman_r': best_spearman_row['spearman_r'],
            'best_spearman_p': best_spearman_row['spearman_p'],
            'spearman_matches_assigned': best_spearman_row['compared_bin'] == assigned_bin
        })

    summary_df = pd.DataFrame(summary_data)
    summary_output = os.path.join(output_dir, 'correlation_summary.tsv')
    summary_df.to_csv(summary_output, sep='\t', index=False, float_format='%.6f')
    print(f"  Summary saved to: {summary_output}")

    # Print statistics
    print("\n=== Analysis Summary ===")
    print(f"Total symbiosis contigs analyzed: {len(summary_df)}")
    print(f"Total bins compared: {len(bin_coverage_profiles)}")

    binned_contigs = summary_df[summary_df['assigned_bin'] != 'unbinned']
    if len(binned_contigs) > 0:
        pearson_match_rate = binned_contigs['pearson_matches_assigned'].mean() * 100
        spearman_match_rate = binned_contigs['spearman_matches_assigned'].mean() * 100
        print(f"Binned symbiosis contigs: {len(binned_contigs)}")
        print(f"  Pearson best match = assigned bin: {pearson_match_rate:.1f}%")
        print(f"  Spearman best match = assigned bin: {spearman_match_rate:.1f}%")

    unbinned_contigs = summary_df[summary_df['assigned_bin'] == 'unbinned']
    print(f"Unbinned symbiosis contigs: {len(unbinned_contigs)}")

    # High correlation summary
    high_pearson = results_df[results_df['pearson_r'] > 0.8]
    high_spearman = results_df[results_df['spearman_r'] > 0.8]
    print(f"\nHigh correlations (r > 0.8):")
    print(f"  Pearson: {len(high_pearson)} contig-bin pairs")
    print(f"  Spearman: {len(high_spearman)} contig-bin pairs")
else:
    print("  No correlation results to save")

# Also save bin coverage profiles for reference
bin_profiles_df = pd.DataFrame(bin_coverage_profiles).T
bin_profiles_df.index.name = 'bin'
bin_profiles_output = os.path.join(output_dir, 'bin_coverage_profiles.tsv')
bin_profiles_df.to_csv(bin_profiles_output, sep='\t', float_format='%.6f')
print(f"  Bin coverage profiles saved to: {bin_profiles_output}")

# Save symbiosis contig coverage profiles
symbiosis_profiles = symbiosis_coverage[['contig_id', 'length'] + sample_names].copy()
symbiosis_profiles_output = os.path.join(output_dir, 'symbiosis_contig_coverage.tsv')
symbiosis_profiles.to_csv(symbiosis_profiles_output, sep='\t', index=False, float_format='%.6f')
print(f"  Symbiosis contig coverage saved to: {symbiosis_profiles_output}")

print("\n=== Correlation analysis completed ===")
PYTHON_EOF

    # Run the analysis script
    python3 "$analysis_script" \
        "$coverage_file" \
        "$symbiosis_contigs_list" \
        "$bins_dir" \
        "$output_dir" \
        "$treatment" \
        2>&1 | tee "${LOG_DIR}/${treatment}_correlation_analysis.log"

    local exit_code=${PIPESTATUS[0]}

    if [ $exit_code -eq 0 ]; then
        log "  Correlation analysis completed successfully"
        return 0
    else
        log "ERROR: Correlation analysis failed"
        return 1
    fi
}

# Create final summary report
create_summary_report() {
    local output_dir="$1"
    local treatment="$2"

    local report_file="${output_dir}/symbiosis_analysis_report.txt"

    cat > "$report_file" << EOF
Symbiosis Contig Analysis Report
================================

Treatment: $treatment
Date: $(date)

Analysis Pipeline:
------------------
1. Prodigal protein prediction from coassembly contigs
2. InterProScan functional annotation
3. CoverM contig coverage (trimmed_mean) calculation
4. Nitrogen fixation/nodulation contig filtering
5. Bin coverage profile calculation (length-weighted mean)
6. Pearson and Spearman correlation analysis

Search Terms Used:
------------------
- nif (nifH, nifD, nifK, nifE, nifN, nifB, etc.)
- nitrogenase
- nod-box / nodbox
- nodulation
- nitrogen fixation
- fix genes (fixA, fixB, fixC, etc.)

Output Files:
-------------
EOF

    # Add file listings
    echo "" >> "$report_file"
    echo "Generated files:" >> "$report_file"

    for f in "${output_dir}"/*.{tsv,fasta,faa,txt} 2>/dev/null; do
        if [ -f "$f" ]; then
            local size=$(du -h "$f" | cut -f1)
            local lines=$(wc -l < "$f" 2>/dev/null || echo "N/A")
            echo "  $(basename "$f"): $size, $lines lines" >> "$report_file"
        fi
    done

    # Add summary statistics if available
    if [ -f "${output_dir}/correlation_summary.tsv" ]; then
        echo "" >> "$report_file"
        echo "Summary Statistics:" >> "$report_file"
        echo "-------------------" >> "$report_file"

        local total_contigs=$(tail -n +2 "${output_dir}/correlation_summary.tsv" | wc -l)
        echo "  Total symbiosis contigs: $total_contigs" >> "$report_file"
    fi

    log "Summary report created: $report_file"
}

###############################################################################
# MAIN EXECUTION
###############################################################################

# This script only supports coassembly mode
if [ "${ASSEMBLY_MODE:-coassembly}" != "coassembly" ]; then
    log "ERROR: This script currently only supports coassembly mode"
    log "  Set ASSEMBLY_MODE=coassembly in your environment"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Get treatment from array task
TREATMENTS_ARRAY=($(get_treatments))
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ $TASK_ID -ge ${#TREATMENTS_ARRAY[@]} ]; then
    log "Array task $TASK_ID has no treatment to process"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

TREATMENT="${TREATMENTS_ARRAY[$TASK_ID]}"
export TREATMENT

log "====== Starting Symbiosis Contig Analysis for Treatment: $TREATMENT ======"

# Define paths
assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"
binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
binette_dir="${OUTPUT_DIR}/bin_refinement/${TREATMENT}/binette/final_bins"
shared_bam_dir="${binning_dir}/shared_bam_files"
output_dir="${OUTPUT_DIR}/symbiosis_analysis/${TREATMENT}"

mkdir -p "$output_dir"

# Find assembly file
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
log "Binette bins: $binette_dir"
log "BAM files: $shared_bam_dir"

# Check prerequisites
if [ ! -d "$binette_dir" ]; then
    log "ERROR: Binette bins directory not found: $binette_dir"
    log "  Please run stage 4b (Binette) first"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

if [ ! -d "$shared_bam_dir" ]; then
    log "ERROR: Shared BAM directory not found: $shared_bam_dir"
    log "  Please run stage 3a (create_shared_bams) first"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# ============================================================================
# Stage 1: Prodigal protein prediction
# ============================================================================
log "=== Stage 1: Prodigal Protein Prediction ==="

proteins_file=$(run_prodigal "$assembly_fasta" "${output_dir}/prodigal" "$TREATMENT")
if [ $? -ne 0 ] || [ -z "$proteins_file" ]; then
    log "ERROR: Prodigal failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# ============================================================================
# Stage 2: InterProScan annotation
# ============================================================================
log "=== Stage 2: InterProScan Annotation ==="

interproscan_tsv=$(run_interproscan "$proteins_file" "${output_dir}/interproscan" "$TREATMENT")
if [ $? -ne 0 ] || [ -z "$interproscan_tsv" ]; then
    log "ERROR: InterProScan failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# ============================================================================
# Stage 3: CoverM contig coverage
# ============================================================================
log "=== Stage 3: CoverM Contig Coverage ==="

coverage_file=$(run_coverm_contig "$assembly_fasta" "$shared_bam_dir" "${output_dir}/coverage" "$TREATMENT")
if [ $? -ne 0 ] || [ -z "$coverage_file" ]; then
    log "ERROR: CoverM contig failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# ============================================================================
# Stage 4: Filter symbiosis contigs
# ============================================================================
log "=== Stage 4: Filter Symbiosis Contigs ==="

filter_result=$(filter_symbiosis_contigs "$interproscan_tsv" "$assembly_fasta" "${output_dir}/symbiosis" "$TREATMENT")
if [ $? -ne 0 ]; then
    log "ERROR: Symbiosis filtering failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

IFS='|' read -r symbiosis_contigs_list symbiosis_annotations symbiosis_fasta <<< "$filter_result"

# ============================================================================
# Stage 5: Correlation analysis
# ============================================================================
log "=== Stage 5: Correlation Analysis ==="

if ! run_correlation_analysis "$coverage_file" "$symbiosis_contigs_list" "$binette_dir" "${output_dir}/correlations" "$TREATMENT"; then
    log "ERROR: Correlation analysis failed"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# ============================================================================
# Create summary report
# ============================================================================
create_summary_report "$output_dir" "$TREATMENT"

# Cleanup
cleanup_temp_dir "$TEMP_DIR"

log "====== Symbiosis Contig Analysis Completed for Treatment: $TREATMENT ======"
