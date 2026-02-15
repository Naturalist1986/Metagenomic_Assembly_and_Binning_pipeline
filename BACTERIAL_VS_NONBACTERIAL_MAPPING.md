# Bacterial vs Non-Bacterial Mapping Analysis

## Overview

This workflow quantifies what percentage of sequencing reads originate from **bacterial** vs **non-bacterial** (archaeal, eukaryotic, unknown, etc.) sources by mapping reads to sequences classified by EukFinder.

## Key Innovation: Per-Binner Mapping

### The Problem
EukFinder runs on bins from **multiple binning tools** (MetaBAT2, MaxBin2, CONCOCT, DAS Tool, Binette, COMEBin, SemiBin). These different binners often produce bins with **overlapping contigs** - the same genomic sequence may appear in bins from different tools.

If we naively combined all sequences from all binners into a single reference and mapped reads to it, we would:
- **Artificially inflate mapping rates** (reads mapping to duplicated contigs counted multiple times)
- **Overestimate** the total fraction of reads explained by binned sequences
- **Cannot determine** which binner's results are most accurate

### The Solution
We map reads to **each binner's sequences separately** and then **average the results**:

1. For each binner (e.g., MetaBAT2):
   - Collect bacterial sequences from that binner's bins only
   - Collect non-bacterial sequences from that binner's bins only
   - Map reads to these binner-specific references
   - Calculate percentages

2. Average percentages across all binners

3. Report both:
   - **Average** (recommended for interpretation)
   - **Per-binner** (to assess consistency)

This approach:
- ✅ Avoids double-counting reads mapping to overlapping contigs
- ✅ Provides fair comparison across binning tools
- ✅ Shows variability across binners (indicates reliability)
- ✅ Gives accurate representation of read origin

## Workflow

### Prerequisites
- EukFinder must be completed on all bins
- Reads must be available (quality-filtered or assembled)

### Usage

**Single command to run everything:**
```bash
export OUTPUT_DIR="/path/to/output"
./submit_bacterial_vs_nonbacterial_mapping.sh
```

This will:
1. Identify all treatments with EukFinder results
2. Identify all binners used (metabat2, maxbin2, concoct, etc.)
3. Submit SLURM array jobs to map reads for each treatment
4. Generate mapping statistics

**After jobs complete:**
```bash
./summarize_bacterial_vs_nonbacterial.sh
```

This creates:
- `bacterial_vs_nonbacterial_summary.tsv` - Spreadsheet-compatible table
- `bacterial_vs_nonbacterial_summary.txt` - Human-readable report

### Output Table Format

```
Treatment    Total_Reads  Bacterial_Reads  Bacterial_%  Non-Bacterial_Reads  Non-Bacterial_%  Unmapped_Reads  Unmapped_%
-----------  -----------  ---------------  -----------  -------------------  ---------------  --------------  ----------
treatment_A  10000000     7500000          75.00%       1500000              15.00%           1000000         10.00%
treatment_B  12000000     9000000          75.00%       2000000              16.67%           1000000         8.33%
```

### Interpretation

**Bacterial %**: Percentage of reads mapping to bacterial sequences (Bact category)

**Non-Bacterial %**: Percentage of reads mapping to:
- Archaeal (Arch)
- Eukaryotic (Euk)
- Unknown (Unk)
- Eukaryotic + Unknown (EUnk)
- Miscellaneous (Misc)

**Unmapped %**: Percentage of reads that don't map to any binned sequences
- May indicate:
  - Low-abundance organisms not assembled
  - Host contamination (if metagenomic sample)
  - Quality-filtered reads
  - Sequencing errors

## Technical Details

### Mapping Parameters
- **Tool**: BBMap
- **Identity threshold**: 95% (minid=0.95)
- **Ambiguous reads**: Random assignment
- **Resources**: 16 CPUs, 64GB RAM per job

### Directory Structure
```
${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping/
  treatment1/
    metabat2_bacterial.fasta         # Bacterial sequences from MetaBAT2 bins
    metabat2_nonbacterial.fasta      # Non-bacterial sequences from MetaBAT2 bins
    metabat2_bacterial.bam           # Reads mapped to MetaBAT2 bacterial
    metabat2_nonbacterial.bam        # Reads mapped to MetaBAT2 non-bacterial
    metabat2_bacterial_stats.txt     # BBMap stats for bacterial mapping
    maxbin2_bacterial.fasta          # Bacterial sequences from MaxBin2 bins
    maxbin2_nonbacterial.fasta       # Non-bacterial sequences from MaxBin2 bins
    ...                              # Same for all other binners
    mapping_summary.txt              # Final averaged results
  treatment2/
    ...
  bacterial_vs_nonbacterial_summary.tsv   # Cross-treatment table
  bacterial_vs_nonbacterial_summary.txt   # Formatted report
```

### Per-Binner Results
Each treatment's `mapping_summary.txt` includes:
```
Results (Averaged Across Binners):
-----------------------------------
Bacterial (Bact):
  Mapped reads: 7500000
  Percentage: 75.00%

Non-Bacterial (Arch + Euk + Unk + EUnk + Misc):
  Mapped reads: 1500000
  Percentage: 15.00%

Unmapped:
  Unmapped reads: 1000000
  Percentage: 10.00%

Per-Binner Results:
-------------------

metabat2:
  Bacterial: 7600000 (76.00%)
  Non-bacterial: 1400000 (14.00%)
  Unmapped: 1000000 (10.00%)

maxbin2:
  Bacterial: 7400000 (74.00%)
  Non-bacterial: 1600000 (16.00%)
  Unmapped: 1000000 (10.00%)
```

**What to look for:**
- **Low variability** across binners (good) = consistent classification
- **High variability** across binners = uncertainty in binning results
- **One binner very different** = may indicate binner-specific issues

## Example Workflow

```bash
# 1. Set output directory
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"

# 2. Run EukFinder on all bins (if not already done)
./submit_eukfinder_all_bins.sh

# 3. Wait for EukFinder to complete
squeue -u $USER | grep eukfinder

# 4. Run bacterial vs non-bacterial mapping
./submit_bacterial_vs_nonbacterial_mapping.sh

# 5. Monitor mapping jobs
squeue -u $USER | grep bact_map

# 6. After mapping completes, generate summary table
./summarize_bacterial_vs_nonbacterial.sh

# 7. View results
cat bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.txt

# 8. Import TSV into Excel/R for visualization
# File: bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.tsv
```

## Troubleshooting

### No treatments found
- Ensure EukFinder has completed successfully
- Check `${OUTPUT_DIR}/eukfinder_output/` for treatment directories

### No binners identified
- Check EukFinder directory structure
- Verify bin directories follow naming: `{sample}_{binner}_{binname}/`

### Could not find reads
- Script searches multiple locations:
  - `${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads/`
  - `${OUTPUT_DIR}/assembly/${TREATMENT}/`
  - `${OUTPUT_DIR}/quality_filtering/${TREATMENT}/`
  - `${OUTPUT_DIR}/validate_repair/${TREATMENT}/`
  - `${INPUT_DIR}/${TREATMENT}/`
- Ensure reads are in one of these locations

### BBMap memory issues
- Increase `--mem` in submission script if needed
- Default: 64GB (sufficient for most datasets)

## Related Scripts

- `submit_eukfinder_all_bins.sh` - Run EukFinder on all bins (prerequisite)
- `summarize_bacterial_vs_nonbacterial.sh` - Generate comparison table (post-processing)
- `10_eukfinder.sh` - Core EukFinder processing script

## References

- **EukFinder**: Classifies contigs as bacterial, archaeal, eukaryotic, or unknown
- **BBMap**: Fast and accurate read mapping tool
- **Binning tools**: MetaBAT2, MaxBin2, CONCOCT, DAS Tool, Binette, COMEBin, SemiBin

## Citation

If you use this workflow, please cite:
- EukFinder publication (if available)
- BBMap: Bushnell, B. (2014). BBMap. sourceforge.net/projects/bbmap/
- Relevant binning tool papers

## Contact

For issues or questions:
- File an issue in the pipeline repository
- Contact the pipeline maintainer
