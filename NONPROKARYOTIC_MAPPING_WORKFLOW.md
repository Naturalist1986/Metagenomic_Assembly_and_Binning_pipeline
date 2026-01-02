# Non-Prokaryotic Read Mapping Workflow

This workflow analyzes what percentage of metagenomic reads map to non-prokaryotic (potentially eukaryotic) sequences identified by EukFinder.

## Overview

The workflow consists of 4 main scripts:

1. **collect_nonprokaryotic_sequences.sh** - Collects all non-prokaryotic sequences per treatment
2. **map_to_nonprokaryotic.sh** - Maps treatment reads to non-prokaryotic sequences
3. **submit_nonprokaryotic_mapping.sh** - Automated submission wrapper
4. **summarize_nonprokaryotic_mapping.sh** - Creates final summary table

## Prerequisites

- EukFinder analysis must be completed (run `submit_eukfinder.sh` first)
- Treatment-level coassembly reads must be available
- Bowtie2 and SAMtools must be installed (available in metawrap-env)

## What Gets Collected

From each EukFinder result, the workflow collects sequences from:
- **Euk.fasta** - Eukaryotic sequences
- **Unk.fasta** - Unknown sequences
- **EUnk.fasta** - Eukaryotic + Unknown sequences
- **Misc.fasta** - Miscellaneous sequences

**Excluded** (prokaryotic):
- Bact.fasta - Bacterial sequences
- Arch.fasta - Archaeal sequences

## Quick Start

```bash
# Run the complete workflow with one command
./submit_nonprokaryotic_mapping.sh
```

This will:
1. Collect non-prokaryotic sequences from all EukFinder results
2. Submit SLURM jobs to map reads to these sequences
3. Tell you how to generate the summary when jobs complete

After jobs finish:
```bash
# Generate the final summary
./summarize_nonprokaryotic_mapping.sh
```

## Manual Workflow

### Step 1: Collect Non-Prokaryotic Sequences

```bash
./collect_nonprokaryotic_sequences.sh
```

This creates combined FASTA files per treatment:
```
nonprokaryotic_sequences/
├── carK_nonprokaryotic.fasta
├── carR_nonprokaryotic.fasta
├── ces_nonprokaryotic.fasta
├── hok_nonprokaryotic.fasta
├── mtz_nonprokaryotic.fasta
├── RH_nonprokaryotic.fasta
└── collection_summary.txt
```

### Step 2: Map Reads

Submit the mapping jobs:

```bash
# Determine number of treatments
NUM_TREATMENTS=$(ls -1 nonprokaryotic_sequences/*_nonprokaryotic.fasta | wc -l)
MAX_INDEX=$((NUM_TREATMENTS - 1))

# Submit array job
sbatch --array=0-${MAX_INDEX} --account=ofinkel map_to_nonprokaryotic.sh
```

Or use the wrapper:
```bash
./submit_nonprokaryotic_mapping.sh
```

### Step 3: Generate Summary

After mapping jobs complete:

```bash
./summarize_nonprokaryotic_mapping.sh
```

## Output Structure

```
nonprokaryotic_sequences/
├── <treatment>_nonprokaryotic.fasta  # Combined sequences per treatment
└── collection_summary.txt

nonprokaryotic_mapping/
├── <treatment>/
│   ├── bowtie2_index/               # Bowtie2 index files
│   ├── <treatment>.sorted.bam        # Sorted alignment
│   ├── <treatment>.sorted.bam.bai    # BAM index
│   ├── mapping_stats.txt             # Detailed statistics
│   └── mapping_summary.tsv           # Simple TSV summary
├── nonprokaryotic_mapping_summary.tsv  # Combined results (all treatments)
└── nonprokaryotic_mapping_report.txt   # Human-readable report
```

## Output Files

### nonprokaryotic_mapping_summary.tsv

Tab-separated table with columns:
- **Treatment** - Treatment name
- **Total_Read_Pairs** - Total number of read pairs
- **Mapped_Read_Pairs** - Number mapping to non-prokaryotic sequences
- **Percent_Mapped** - Percentage of reads mapped
- **Nonprok_Sequences** - Number of non-prokaryotic sequences
- **Nonprok_Size_bp** - Total size of non-prokaryotic sequences

### nonprokaryotic_mapping_report.txt

Human-readable report including:
- Overall statistics across all treatments
- Per-treatment results table
- Top 5 treatments by mapping percentage
- Detailed statistics for each treatment

## Interpreting Results

**High mapping percentage (>5%)**:
- Indicates substantial non-prokaryotic (potentially eukaryotic) content
- These samples may have significant eukaryotic contamination or symbiosis

**Low mapping percentage (<1%)**:
- Most reads are prokaryotic (bacteria/archaea)
- Typical for many microbial communities

**Moderate mapping percentage (1-5%)**:
- Some eukaryotic presence
- Could be environmental eukaryotes, protists, or fungi

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER | grep nonprok_map

# View output logs
tail -f ${OUTPUT_DIR}/logs/nonprokaryotic_mapping/map_*.out

# Check for errors
tail -f ${OUTPUT_DIR}/logs/nonprokaryotic_mapping/map_*.err

# Check completed treatments
ls ${OUTPUT_DIR}/nonprokaryotic_mapping/*/.mapping_complete
```

## Configuration

SLURM parameters can be modified in `map_to_nonprokaryotic.sh`:

```bash
#SBATCH --cpus-per-task=16    # Threads for mapping
#SBATCH --mem=64G              # Memory allocation
#SBATCH --time=12:00:00        # Time limit
#SBATCH --account=ofinkel      # SLURM account
```

Bowtie2 mapping parameters (in `map_to_nonprokaryotic.sh`):
- `--very-sensitive` - High sensitivity mapping mode
- `--no-unal` - Don't output unmapped reads (saves space)

## Troubleshooting

### No non-prokaryotic sequences found
- Ensure EukFinder has completed successfully
- Check that EukFinder results exist in `${OUTPUT_DIR}/eukfinder/`
- Verify EukFinder identified some non-prokaryotic sequences

### Read files not found
The script searches for reads in:
1. `${OUTPUT_DIR}/assembly/<treatment>/filtered_1.fastq.gz`
2. `${OUTPUT_DIR}/quality_filtering/<treatment>/filtered_1.fastq.gz`
3. `${OUTPUT_DIR}/coassembly_reads/<treatment>_1.fastq.gz`

Ensure your coassembly reads are in one of these locations.

### Bowtie2 index build fails
- Check that FASTA files are not empty
- Ensure sufficient disk space
- Verify bowtie2 is installed: `conda activate metawrap-env; bowtie2 --version`

### Out of memory errors
- Increase `--mem` in the SLURM directives
- Reduce `--cpus-per-task` to use less memory
- Check memory usage: `seff <job_id>`

## Example Output

```
Non-Prokaryotic Read Mapping Summary Report
============================================
Generated: Thu Jan  2 2026 10:30:00 IST

Total treatments analyzed: 6

Overall Statistics:
-------------------
Total read pairs across all treatments: 45000000
Total read pairs mapped to non-prokaryotic sequences: 850000
Overall percentage mapped: 1.8889%

Per-Treatment Results:
----------------------
Treatment       Total Reads    Mapped Reads    % Mapped  Nonprok Seqs   Nonprok Size
---------       -----------    ------------    --------  ------------   ------------
hok                 7500000          180000      2.4000%          89        3223345
carK                9000000          210000      2.3333%         210        4404962
RH                  8500000          145000      1.7059%          61         319327
carR                6800000          102000      1.5000%          76        1430010
ces                 7200000           94000      1.3056%          69         132407
mtz                 6000000           62000      1.0333%          10          17272

Top 5 Treatments by Mapping Percentage:
----------------------------------------
Treatment         % Mapped
---------         --------
hok                 2.4000%
carK                2.3333%
RH                  1.7059%
carR                1.5000%
ces                 1.3056%
```

## Integration with EukFinder Workflow

This workflow is designed to run after EukFinder:

```bash
# 1. Run EukFinder on largest bins
./submit_eukfinder.sh

# 2. Wait for EukFinder to complete, then generate summary
./create_eukfinder_summary_table.sh

# 3. Map reads to non-prokaryotic sequences
./submit_nonprokaryotic_mapping.sh

# 4. After mapping completes, generate summary
./summarize_nonprokaryotic_mapping.sh
```

## Advanced Usage

### Reprocess specific treatment

```bash
# Set the treatment manually
export SLURM_ARRAY_TASK_ID=0  # First treatment (0-indexed)
bash map_to_nonprokaryotic.sh
```

### Custom read locations

If your reads are in a non-standard location, edit the read search paths in `map_to_nonprokaryotic.sh` around lines 70-85.

### Different mapping sensitivity

Edit `map_to_nonprokaryotic.sh` and change bowtie2 parameters:
- `--very-sensitive` → `--sensitive` (faster, less sensitive)
- `--very-sensitive` → `--very-sensitive-local` (local alignment)

## References

- Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/
- SAMtools: http://www.htslib.org/
- EukFinder: https://github.com/RogerLab/Eukfinder
