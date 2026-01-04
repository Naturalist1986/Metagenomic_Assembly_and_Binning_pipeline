# EukFinder Workflow for ALL Bins

This workflow runs EukFinder on **ALL bins** from each binning tool (not just the largest ones), then maps treatment-level coassembly reads to non-prokaryotic sequences.

## Workflow Overview

```
1. Run EukFinder on ALL bins
   ↓
2. Generate consolidated results table
   ↓
3. Map reads to non-prokaryotic sequences
   ↓
4. Get summary of % reads mapped
```

## Quick Start

### Complete Workflow (One Command)

```bash
# Set your output directory
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"

# Step 1: Run EukFinder on ALL bins
./submit_eukfinder_all_bins.sh

# Wait for jobs to complete, then:

# Step 2: Generate EukFinder summary table
./create_eukfinder_summary_table.sh

# Step 3: Map reads to non-prokaryotic sequences and get summary
./submit_nonprokaryotic_mapping.sh

# After mapping jobs complete:

# Step 4: Generate mapping summary
./summarize_nonprokaryotic_mapping.sh
```

## Detailed Steps

### Step 1: Run EukFinder on ALL Bins

```bash
./submit_eukfinder_all_bins.sh
```

This will:
- Identify ALL bins from binning results (not just largest)
- Submit SLURM array jobs with `--account=ofinkel`
- Process each bin with EukFinder

**Expected output structure:**
```
binning/
└── carK/
    ├── concoct_bins/
    │   ├── bin.1.fa
    │   ├── bin.2.fa
    │   └── bin.3.fa   (ALL bins processed)
    ├── metabat2_bins/
    └── maxbin2_bins/
```

**Results saved to:**
```
eukfinder/
└── carK/
    ├── carK_concoct_bin.1/
    │   └── Eukfinder_results/
    ├── carK_concoct_bin.2/
    │   └── Eukfinder_results/
    └── carK_concoct_bin.3/
        └── Eukfinder_results/
```

### Step 2: Generate EukFinder Summary Table

After EukFinder jobs complete:

```bash
./create_eukfinder_summary_table.sh
```

**Output:** `eukfinder_consolidated_results.tsv`

Columns include:
- Treatment, Sample, Binner, Bin
- Total sequences and size
- Euk, Bact, Arch, Unk, EUnk counts and percentages

### Step 3: Map Reads to Non-Prokaryotic Sequences

```bash
./submit_nonprokaryotic_mapping.sh
```

This will:
1. **Collect** all non-prokaryotic sequences (Euk, Unk, EUnk) from ALL bins per treatment
2. **Map** treatment coassembly reads (`coassembly/<treatment>/merged_reads/merged_R1.fastq.gz`)
3. **Calculate** what percentage of reads mapped

### Step 4: Generate Mapping Summary

After mapping jobs complete:

```bash
./summarize_nonprokaryotic_mapping.sh
```

**Output:** `nonprokaryotic_mapping_summary.tsv`

Shows for each treatment:
- Total read pairs
- Reads mapped to non-prokaryotic sequences
- **Percentage mapped** ← Key metric!

## Differences from Largest Bins Workflow

| Feature | Largest Bins | ALL Bins |
|---------|-------------|----------|
| Script | `submit_eukfinder.sh` | `submit_eukfinder_all_bins.sh` |
| Bins processed | 1-2 largest per binner | ALL bins |
| Number of jobs | ~6-12 per treatment | Could be 100s |
| Time | Faster | Longer |
| Coverage | Partial | Complete |

## Example: carK Treatment

### Largest Bins Workflow
```
carK/
├── metabat2: 2 largest bins
├── maxbin2: 2 largest bins
└── concoct: 2 largest bins
Total: ~6 bins
```

### ALL Bins Workflow
```
carK/
├── metabat2: ALL bins (e.g., 45 bins)
├── maxbin2: ALL bins (e.g., 38 bins)
└── concoct: ALL bins (e.g., 52 bins)
Total: ~135 bins
```

## Monitoring

```bash
# Check EukFinder jobs
squeue -u $USER | grep eukfinder_all

# Check mapping jobs
squeue -u $USER | grep nonprok_map

# View logs
tail -f ${OUTPUT_DIR}/logs/eukfinder/eukfinder_all_*.out
tail -f ${OUTPUT_DIR}/logs/nonprokaryotic_mapping/map_*.out
```

## Expected Output

### EukFinder Summary (`eukfinder_consolidated_results.tsv`)

```
Treatment  Binner    Bin        Total_Seqs  Euk_Seqs  Euk_%
carK       concoct   bin.1      5234        210       4.01%
carK       concoct   bin.2      8921        45        0.50%
carK       concoct   bin.3      3456        12        0.35%
...
```

### Mapping Summary (`nonprokaryotic_mapping_summary.tsv`)

```
Treatment  Total_Reads  Mapped_Reads  Percent_Mapped
carK       9000000      210000        2.33%
carR       6800000      102000        1.50%
ces        7200000      94000         1.31%
hok        7500000      180000        2.40%
mtz        6000000      62000         1.03%
RH         8500000      145000        1.71%
```

## Interpreting Results

**High mapping percentage (>2%):**
- Significant non-prokaryotic content
- Multiple bins contain eukaryotic sequences
- Treatment may have eukaryotic contamination/symbiosis

**Individual bin analysis:**
- Bins with high Euk_% are likely eukaryotic organisms
- Bins with high Unk_% may be novel or divergent organisms
- EUnk category captures potential eukaryotes with uncertain classification

## Resource Requirements

For ALL bins workflow:

**Per job:**
- CPUs: 48 cores
- Memory: 128 GB
- Time: up to 24 hours
- Account: `ofinkel`

**Estimated totals** (example with 600 bins):
- Jobs: ~600 EukFinder jobs + 6 mapping jobs
- Wall time: Jobs run in parallel (no limit on concurrency)
- Storage: ~1-5 GB per bin for EukFinder results

## Troubleshooting

### Too many bins - jobs taking too long

If you have thousands of bins, consider:
1. Use largest bins workflow instead: `submit_eukfinder.sh`
2. Filter by bin size first:
   ```bash
   # Only process bins > 1 MB
   awk -F'|' '$6 > 1000000' all_bins_list.txt > filtered_bins_list.txt
   export BINS_LIST_FILE="filtered_bins_list.txt"
   sbatch --array=0-N 10_eukfinder.sh
   ```

### Database locking errors

Run taxonomy update first:
```bash
./update_eukfinder_taxonomy.sh
```

### Out of memory

Some very large bins may need more memory. Increase in `submit_eukfinder_all_bins.sh`:
```bash
#SBATCH --mem=256G
```

## Files Created

```
all_bins_list.txt                           # List of ALL bins for EukFinder
eukfinder/<treatment>/<bin>/               # EukFinder results per bin
eukfinder_consolidated_results.tsv         # Summary of all EukFinder results
nonprokaryotic_sequences/<treatment>_nonprokaryotic.fasta  # Combined non-prok seqs
nonprokaryotic_mapping/<treatment>/        # Mapping results per treatment
nonprokaryotic_mapping_summary.tsv         # Final % mapped summary
```

## Best Practices

1. **Run on all bins** when:
   - You need comprehensive analysis
   - You want to find all potential eukaryotes
   - Computational resources are available

2. **Use largest bins** when:
   - Quick screening needed
   - Limited compute time
   - Most eukaryotic content is in large bins

3. **Filter bins** by:
   - Size (>1 MB recommended for meaningful analysis)
   - Completeness (if CheckM2 data available)
   - Contamination levels

## Next Steps After Analysis

1. **Identify high-eukaryotic bins:**
   ```bash
   # Show bins with >1% eukaryotic sequences
   awk -F'\t' '$9 > 1.0' eukfinder_consolidated_results.tsv | sort -t$'\t' -k9 -nr
   ```

2. **Extract eukaryotic sequences:**
   - Look in `Eukfinder_results/Euk.fasta` for each high-scoring bin
   - Perform taxonomic classification (e.g., BLAST, DIAMOND)
   - Assemble into eukaryotic MAGs if desired

3. **Analyze eukaryotic content:**
   - Gene prediction (e.g., Augustus, GeneMark-ES)
   - Functional annotation
   - Phylogenetic analysis

## See Also

- `EUKFINDER_WORKFLOW.md` - Original workflow (largest bins only)
- `NONPROKARYOTIC_MAPPING_WORKFLOW.md` - Detailed mapping documentation
- EukFinder GitHub: https://github.com/RogerLab/Eukfinder
