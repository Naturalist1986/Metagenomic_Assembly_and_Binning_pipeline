# EukFinder Workflow for Largest Bins

This workflow runs EukFinder on the largest bins from each binning tool to identify eukaryotic sequences.

## Overview

The EukFinder workflow consists of three main components:

1. **identify_largest_bins.sh** - Identifies the largest bins from each binning tool
2. **10_eukfinder.sh** - Main SLURM script that runs EukFinder on individual bins
3. **submit_eukfinder.sh** - Submission wrapper that automates the entire process

## Prerequisites

- Binning stage must be completed (03_binning.sh)
- EukFinder conda environment must be installed and named `eukfinder`
- EukFinder databases must be set up and configured

### Installing EukFinder

```bash
# Create conda environment
conda create -n eukfinder -c bioconda eukfinder

# Activate environment
conda activate eukfinder

# Set up databases (first time only)
# Follow instructions at: https://github.com/RogerLab/Eukfinder/wiki
```

### Configuring Database Paths

**IMPORTANT**: You must configure THREE database paths in `00_config_utilities.sh` before running EukFinder.

Edit lines 25-27 in `00_config_utilities.sh`:

```bash
# EukFinder database paths
export EUKFINDER_CENTRIFUGE_DB="/path/to/your/centrifuge/database"
export EUKFINDER_PLAST_DB="/path/to/your/plast/database"
export EUKFINDER_PLAST_ID_MAP="/path/to/your/plast/id_map"
```

Replace the placeholder paths with your actual database locations.

**To find these paths**, check your EukFinder database directory:

```bash
ls -la /path/to/eukfinder_databases/
```

Common structures:
- **Centrifuge DB**: Look for files with `.cf` extensions (e.g., `centrifuge_db.1.cf`)
- **PLAST DB**: Look for FASTA files (e.g., `uniprot_sprot.fasta` or similar)
- **PLAST ID Map**: Look for mapping files (e.g., `idmapping.dat` or `uniprot_idmap.txt`)

Example configuration:

```bash
export EUKFINDER_CENTRIFUGE_DB="/sci/backup/ofinkel/moshea/eukfinder_databases/centrifuge/p_compressed"
export EUKFINDER_PLAST_DB="/sci/backup/ofinkel/moshea/eukfinder_databases/plast/uniprot_sprot.fasta"
export EUKFINDER_PLAST_ID_MAP="/sci/backup/ofinkel/moshea/eukfinder_databases/plast/idmapping.dat"
```

Alternatively, you can set these as environment variables before running:

```bash
export EUKFINDER_CENTRIFUGE_DB="/your/centrifuge/path"
export EUKFINDER_PLAST_DB="/your/plast/path"
export EUKFINDER_PLAST_ID_MAP="/your/id_map/path"
./submit_eukfinder.sh
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
- Activate the correct conda environment: `conda activate eukfinder`
- Verify EukFinder is installed: `conda list | grep eukfinder`

### Database errors
- **IMPORTANT**: Verify database paths are correctly set in `00_config_utilities.sh`
  ```bash
  # Check your current settings
  grep EUKFINDER_.*_DB 00_config_utilities.sh
  ```
- On first run, set `EUKFINDER_TAXONOMY_UPDATE=True`
- Ensure sufficient disk space for databases
- Check EukFinder installation and database setup
- Verify database files exist at specified paths:
  ```bash
  ls -lh $EUKFINDER_CENTRIFUGE_DB
  ls -lh $EUKFINDER_PLAST_DB
  ```

### Out of memory errors
- Reduce `--cpus-per-task` or increase `--mem` in 10_eukfinder.sh
- Consider reducing `EUKFINDER_CHUNKS`

## References

- EukFinder GitHub: https://github.com/RogerLab/Eukfinder
- EukFinder Wiki: https://github.com/RogerLab/Eukfinder/wiki
- Long reads/contigs guide: https://github.com/RogerLab/Eukfinder/wiki/Eukfinder_with_Long_read_or_contig_data
