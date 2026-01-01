# EukFinder Setup Guide for Your System

## Your Configuration Paths

- **Database Location**: `/sci/backup/ofinkel/moshea/eukfinder_databases/`
- **Results Location**: `/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/`

## Step 1: Check Your Database Structure

First, verify what's in your EukFinder database directory:

```bash
ls -la /sci/backup/ofinkel/moshea/eukfinder_databases/
```

You need to identify:
- **Centrifuge database files** (usually have extensions like `.cf`, `.1.cf`, `.2.cf`, `.3.cf`)
- **PLAST database files** (usually FASTA files or database indices)

Common structures:
```
# Option A: Separate subdirectories
eukfinder_databases/
├── centrifuge/
│   ├── database.1.cf
│   ├── database.2.cf
│   └── database.3.cf
└── plast/
    └── database.fasta

# Option B: All in one directory
eukfinder_databases/
├── centrifuge_db.1.cf
├── centrifuge_db.2.cf
├── centrifuge_db.3.cf
└── plast_db.fasta

# Option C: EukFinder default structure
eukfinder_databases/
├── Centrifuge_data/
└── PLAST_data/
```

## Step 2: Configure Database Paths in 00_config_utilities.sh

Based on what you find, edit `00_config_utilities.sh` lines 24-26:

### If you have separate subdirectories:
```bash
# EukFinder database paths
export EUKFINDER_CENTRIFUGE_DB="${EUKFINDER_CENTRIFUGE_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/centrifuge}"
export EUKFINDER_PLAST_DB="${EUKFINDER_PLAST_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/plast}"
```

### If using EukFinder default structure:
```bash
# EukFinder database paths
export EUKFINDER_CENTRIFUGE_DB="${EUKFINDER_CENTRIFUGE_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/Centrifuge_data}"
export EUKFINDER_PLAST_DB="${EUKFINDER_PLAST_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/PLAST_data}"
```

### If all databases are in the same directory:
```bash
# EukFinder database paths (point to the directory containing both)
export EUKFINDER_CENTRIFUGE_DB="${EUKFINDER_CENTRIFUGE_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases}"
export EUKFINDER_PLAST_DB="${EUKFINDER_PLAST_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases}"
```

## Step 3: Set OUTPUT_DIR Environment Variable

Before running the workflow, set your OUTPUT_DIR to point to your results:

```bash
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"
```

You can add this to your `~/.bashrc` or set it each time before running.

## Step 4: Verify Your Binning Results Exist

Check that your binning results are in the expected location:

```bash
ls -la ${OUTPUT_DIR}/binning/
```

Expected structure:
```
binning/
├── <treatment1>/
│   ├── <sample1>/
│   │   ├── metabat2_bins/
│   │   │   ├── bin.1.fa
│   │   │   └── bin.2.fa
│   │   ├── maxbin2_bins/
│   │   └── concoct_bins/
│   └── <sample2>/
└── <treatment2>/
```

## Step 5: Run the EukFinder Workflow

Once everything is configured:

```bash
# Navigate to the pipeline directory
cd /path/to/Metagenomic_Assembly_and_Binning_pipeline

# Set the OUTPUT_DIR
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"

# Run the submission script (finds 2 largest bins per binner by default)
./submit_eukfinder.sh

# Or specify a different number of bins (e.g., 3 largest)
./submit_eukfinder.sh 3
```

## Step 6: Monitor the Jobs

```bash
# Check job status
squeue -u $USER | grep eukfinder

# View output logs (replace JOBID with your actual job ID)
tail -f ${OUTPUT_DIR}/logs/eukfinder/eukfinder_JOBID_*.out

# Check for errors
tail -f ${OUTPUT_DIR}/logs/eukfinder/eukfinder_JOBID_*.err
```

## Step 7: Check Results

Results will be in:
```bash
${OUTPUT_DIR}/eukfinder/<treatment>/<sample>/
```

Each bin analyzed will have:
- `Eukfinder_results/` - Classified sequences (Euk.fasta, Bact.fasta, etc.)
- `Intermediate_data/` - Raw classification data
- `eukfinder_summary.txt` - Summary statistics

## Quick Start Commands

Here's everything in one script:

```bash
#!/bin/bash

# Set your paths
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"

# Verify database structure first
echo "Checking database structure..."
ls -la /sci/backup/ofinkel/moshea/eukfinder_databases/

# Verify binning results exist
echo "Checking binning results..."
ls -la ${OUTPUT_DIR}/binning/

# After verifying, update 00_config_utilities.sh with correct database paths
# Then run:
cd /path/to/Metagenomic_Assembly_and_Binning_pipeline
./submit_eukfinder.sh
```

## Troubleshooting

### If databases are not found:
```bash
# Check if database paths are correct
echo $EUKFINDER_CENTRIFUGE_DB
echo $EUKFINDER_PLAST_DB

# Verify files exist
ls -lh $EUKFINDER_CENTRIFUGE_DB
ls -lh $EUKFINDER_PLAST_DB
```

### If no bins are found:
```bash
# Check that OUTPUT_DIR is set correctly
echo $OUTPUT_DIR

# Verify binning directory structure
find ${OUTPUT_DIR}/binning -name "*.fa" | head -20
```

### If you need to update taxonomy on first run:
Add to your environment before running:
```bash
export EUKFINDER_TAXONOMY_UPDATE=True
```

## Next Steps

1. Log into your HPC system
2. Check your database structure (Step 1)
3. Update `00_config_utilities.sh` with correct paths (Step 2)
4. Set `OUTPUT_DIR` environment variable (Step 3)
5. Verify binning results (Step 4)
6. Run `./submit_eukfinder.sh` (Step 5)
