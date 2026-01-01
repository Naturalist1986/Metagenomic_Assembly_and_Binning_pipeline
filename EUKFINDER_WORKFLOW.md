# EukFinder Workflow for Largest Bins

This workflow runs EukFinder on the largest bins from each binning tool to identify eukaryotic sequences.

## Overview

The EukFinder workflow consists of three main components:

1. **identify_largest_bins.sh** - Identifies the largest bins from each binning tool
2. **10_eukfinder.sh** - Main SLURM script that runs EukFinder on individual bins
3. **submit_eukfinder.sh** - Submission wrapper that automates the entire process

## Prerequisites

- Binning stage must be completed (03_binning.sh)
- EukFinder conda environment must be installed and named `Eukfinder`
- EukFinder databases must be set up

### Installing EukFinder

```bash
# Create conda environment
conda create -n Eukfinder -c bioconda eukfinder

# Activate environment
conda activate Eukfinder

# Set up databases (first time only)
# Follow instructions at: https://github.com/RogerLab/Eukfinder/wiki
```

## Usage

### Quick Start (Recommended)

Run the submission wrapper which handles everything automatically:

```bash
# Submit EukFinder jobs for the 2 largest bins per binner (default)
./submit_eukfinder.sh

# Or specify a different number of largest bins (e.g., 3)
./submit_eukfinder.sh 3
```

This script will:
1. Identify the largest bins from each binning tool
2. Calculate the correct SLURM array size
3. Submit the job array automatically

### Manual Workflow

If you prefer to run steps manually:

#### Step 1: Identify largest bins

```bash
# Find the 2 largest bins from each binner
./identify_largest_bins.sh 2

# Output will be saved to: ${OUTPUT_DIR}/largest_bins_list.txt
```

#### Step 2: Submit EukFinder jobs

```bash
# Check how many bins were found
TOTAL_BINS=$(grep -v "^#" ${OUTPUT_DIR}/largest_bins_list.txt | wc -l)
MAX_INDEX=$((TOTAL_BINS - 1))

# Submit array job with correct indices
sbatch --array=0-${MAX_INDEX}%10 10_eukfinder.sh
```

## Configuration

EukFinder parameters can be customized in `00_config_utilities.sh` or via environment variables:

```bash
# Number of threads (default: 48)
export EUKFINDER_THREADS=48

# Number of chunks for PLAST processing (default: 6)
export EUKFINDER_CHUNKS=6

# Update taxonomy database - set to True on first run (default: False)
export EUKFINDER_TAXONOMY_UPDATE=False

# E-value threshold for PLAST (default: 0.01)
export EUKFINDER_EVALUE=0.01

# Percent identity threshold (default: 60)
export EUKFINDER_PID=60

# Percent coverage threshold (default: 30)
export EUKFINDER_COV=30

# Minimum hit length for Centrifuge (default: 100)
export EUKFINDER_MHLEN=100
```

## Output Structure

Results are organized by treatment and sample:

```
${OUTPUT_DIR}/
└── eukfinder/
    └── <treatment>/
        └── <sample>/
            ├── Eukfinder_results/
            │   ├── Euk.fasta      # Eukaryotic sequences
            │   ├── Bact.fasta     # Bacterial sequences
            │   ├── Arch.fasta     # Archaeal sequences
            │   ├── Unk.fasta      # Unknown sequences
            │   ├── EUnk.fasta     # Eukaryotic unknown
            │   └── Misc.fasta     # Miscellaneous
            ├── Intermediate_data/
            │   ├── centrifuge_classification.txt
            │   └── plast_results/
            └── eukfinder_summary.txt
```

## Output Files

### Eukfinder_results/
Contains FASTA files with sequences classified into categories:
- **Euk.fasta** - Confidently identified eukaryotic sequences
- **Bact.fasta** - Bacterial sequences
- **Arch.fasta** - Archaeal sequences
- **Unk.fasta** - Unknown sequences
- **EUnk.fasta** - Possible eukaryotic sequences with uncertain classification
- **Misc.fasta** - Miscellaneous sequences

### Intermediate_data/
Contains raw classification results from Centrifuge and PLAST

### eukfinder_summary.txt
Summary statistics for each bin analyzed, including:
- Bin information (name, size, contig count)
- Number of sequences in each classification category
- Total sequence sizes for each category

## Monitoring Jobs

```bash
# Check job status
squeue -u $USER | grep eukfinder

# View output logs
tail -f ${LOG_DIR}/eukfinder/eukfinder_*.out

# View error logs
tail -f ${LOG_DIR}/eukfinder/eukfinder_*.err

# Check completed jobs
ls ${OUTPUT_DIR}/eukfinder/*/*/.*.done
```

## Interpreting Results

After EukFinder completes:

1. Check the summary files for overview of classifications
2. Focus on bins with sequences in `Euk.fasta` - these are likely eukaryotic
3. Review `EUnk.fasta` for potential eukaryotic sequences needing further analysis
4. Large numbers of sequences in `Bact.fasta` or `Arch.fasta` indicate prokaryotic contamination

## Troubleshooting

### No bins found
- Ensure binning stage (03_binning.sh) has completed successfully
- Check that bins exist in: `${OUTPUT_DIR}/binning/`
- Verify bin files are not empty

### EukFinder command not found
- Activate the correct conda environment: `conda activate Eukfinder`
- Verify EukFinder is installed: `conda list | grep eukfinder`

### Database errors
- On first run, set `EUKFINDER_TAXONOMY_UPDATE=True`
- Ensure sufficient disk space for databases
- Check EukFinder installation and database setup

### Out of memory errors
- Reduce `--cpus-per-task` or increase `--mem` in 10_eukfinder.sh
- Consider reducing `EUKFINDER_CHUNKS`

## References

- EukFinder GitHub: https://github.com/RogerLab/Eukfinder
- EukFinder Wiki: https://github.com/RogerLab/Eukfinder/wiki
- Long reads/contigs guide: https://github.com/RogerLab/Eukfinder/wiki/Eukfinder_with_Long_read_or_contig_data
