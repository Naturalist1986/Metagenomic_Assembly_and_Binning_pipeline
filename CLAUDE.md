# CLAUDE.md - AI Assistant Guide for Metagenomic Pipeline

This document provides guidance for AI assistants working with the Metagenomic Assembly and Binning Pipeline codebase.

## Project Overview

This is a **production-grade SLURM-based metagenomic pipeline** for processing environmental DNA sequencing data. It takes raw FASTQ reads through assembly, binning, quality assessment, and taxonomic classification to produce high-quality Metagenome-Assembled Genomes (MAGs).

**Key characteristics:**
- 40+ shell scripts totaling ~17,000 lines of code
- Two processing modes: sample-level and treatment-level (coassembly)
- Checkpoint-based job recovery for HPC environments
- Multiple binning algorithms with consensus refinement

## Repository Structure

```
/
├── 00_config_utilities.sh      # Core config & 45+ utility functions (1,191 lines)
├── run_pipeline.sh             # Master orchestration script (1,750 lines)
│
├── Core Pipeline (numbered 00-12):
│   ├── 00_quality_filtering.sh     # Trimmomatic QC
│   ├── 01_validate_repair.sh       # BBTools read repair
│   ├── 02_assembly.sh              # MetaSPAdes (per-sample)
│   ├── 02_coassembly.sh            # MetaSPAdes (per-treatment)
│   ├── 03_binning.sh               # MetaBAT2, MaxBin2, CONCOCT
│   ├── 03a_create_shared_bams.sh   # BAM files for binning
│   ├── 03b_comebin.sh              # COMEBin binner
│   ├── 03c_semibin.sh              # SemiBin binner
│   ├── 04_bin_refinement.sh        # DAS Tool consensus
│   ├── 04b_binette.sh              # Binette alternative
│   ├── 04c_binspreader.sh          # Graph-aware refinement
│   ├── 05_bin_reassembly.sh        # MetaWRAP reassembly
│   ├── 06_magpurify.sh             # Contamination removal
│   ├── 07_checkm2.sh               # Quality assessment
│   ├── 07b_gunc.sh                 # Chimerism detection
│   ├── 08_bin_selection.sh         # Best version selection
│   ├── 09_bin_collection.sh        # Consolidation & GTDB-Tk
│   ├── 10_eukfinder.sh             # Eukaryotic detection
│   ├── 11_symbiosis_contig_analysis.sh  # Nitrogen fixation genes
│   └── 12_coverage_consolidation.sh     # Abundance tables
│
├── Optional/Utility Scripts:
│   ├── -01_merge_lanes.sh          # Multi-lane/run merging
│   ├── plasmid_detection.sh        # PlasClass & MOB-suite
│   ├── annotate_bins.sh            # Functional annotation
│   ├── final_report.sh             # Cross-treatment visualization
│   ├── submit_eukfinder*.sh        # EukFinder job submission
│   └── identify_*_bins.sh          # Bin listing utilities
│
└── Documentation:
    ├── README.md                   # Comprehensive user guide
    ├── SETUP_GUIDE.md              # EukFinder setup
    ├── EUKFINDER_WORKFLOW.md       # EukFinder workflows
    └── ASSEMBLY_SUCCESS_RATE.md    # Metrics documentation
```

## Critical Files to Understand

1. **`00_config_utilities.sh`** - The foundation of the entire pipeline:
   - All configuration variables (paths, thresholds, SLURM settings)
   - 45+ utility functions used by all other scripts
   - Sample storage and management functions
   - Checkpoint system functions
   - Conda environment handling

2. **`run_pipeline.sh`** - Master orchestration:
   - Parses sample sheets (CSV, TSV, Excel)
   - Submits SLURM array jobs in dependency order
   - Handles both processing modes
   - Command-line argument parsing

## Code Conventions

### Shell Script Standards

```bash
# All scripts use strict mode
set -euo pipefail

# Source config as first action
source "$(dirname "$0")/00_config_utilities.sh"

# Use log() function for all output
log "Starting process..."
log "ERROR: Something failed"
```

### Naming Conventions

- **Scripts**: Numbered stages (`00_`, `01_`), variants use letters (`03a_`, `04b_`)
- **Functions**: `snake_case` (e.g., `check_sample_checkpoint`, `activate_env`)
- **Variables**: `UPPER_SNAKE_CASE` for exports, `lower_snake_case` for locals
- **Directories**: Stage outputs follow `{stage_name}/{treatment}/{sample}/` pattern

### Conda Environment Handling

```bash
# CRITICAL: Use the activate_env/deactivate_env functions
# They handle set -u compatibility issues with conda

activate_env "environment_name"
# ... do work ...
deactivate_env

# NEVER do this directly (will fail with set -u):
# conda activate env_name
```

### SLURM Array Job Pattern

All processing scripts follow this pattern:

```bash
#SBATCH --array=0-N%10        # Array size, max 10 parallel
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G

ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

# Sample-level mode:
SAMPLE_INFO=$(get_sample_info_by_index $ARRAY_INDEX)
IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"

# Treatment-level mode:
TREATMENT=$(sed -n "$((ARRAY_INDEX + 1))p" "$TREATMENTS_FILE")
```

### Checkpoint System

```bash
# Check if already completed (skip if done)
if check_sample_checkpoint "$SAMPLE_NAME" "assembly"; then
    log "Skipping - already completed"
    exit 0
fi

# ... do processing ...

# Mark as complete
create_sample_checkpoint "$SAMPLE_NAME" "assembly"
```

### Error Handling Patterns

```bash
# Temp directory with cleanup trap
TEMP_DIR=$(setup_temp_dir)
trap "cleanup_temp_dir '$TEMP_DIR'" EXIT

# Optional variables with defaults (set -u safe)
export VALUE="${VALUE:-default}"

# Check for unset variables before use
if [ -z "${OPTIONAL_VAR:-}" ]; then
    OPTIONAL_VAR="default"
fi
```

## Key Configuration Variables

Located in `00_config_utilities.sh`:

```bash
# Directory paths
INPUT_DIR       # Raw FASTQ input location
OUTPUT_DIR      # Pipeline output root
WORK_DIR        # Processing workspace
CONDA_BASE      # Miniconda installation path

# Database paths
MAGPURIFYDB                # MAGpurify database
EUKFINDER_CENTRIFUGE_DB    # Centrifuge taxonomy DB
EUKFINDER_PLAST_DB         # PLAST protein DB

# Quality thresholds
MIN_COMPLETENESS=90        # CheckM2 completeness threshold
MAX_CONTAMINATION=5        # CheckM2 contamination threshold

# Processing modes
ASSEMBLY_MODE              # "individual" or "coassembly"
TREATMENT_LEVEL_BINNING    # true/false
```

## Processing Modes

### Sample-Level (Default)
- Each sample processed independently
- Array jobs indexed by sample number
- Paths: `{stage}/{treatment}/{sample}/`
- Best for: High-depth samples (>5 GB)

### Treatment-Level (Coassembly)
- Samples merged within treatments, then co-assembled
- Array jobs indexed by treatment number
- Paths: `{stage}/{treatment}/`
- Best for: Low-depth samples, biological replicates
- Enabled with `--treatment-level-binning` flag

## Common Development Tasks

### Adding a New Pipeline Stage

1. Create numbered script (e.g., `XX_new_stage.sh`)
2. Source `00_config_utilities.sh`
3. Add SLURM directives
4. Use `check_sample_checkpoint` / `create_sample_checkpoint`
5. Update `run_pipeline.sh` to submit the new stage

### Adding a New Configuration Variable

1. Add to `00_config_utilities.sh` with default:
   ```bash
   export NEW_VAR="${NEW_VAR:-default_value}"
   ```
2. Document in README.md if user-facing

### Modifying Checkpoint Logic

- Sample checkpoints: `${OUTPUT_DIR}/checkpoints/{treatment}/{sample}/{stage}.done`
- Treatment checkpoints: `${OUTPUT_DIR}/checkpoints/{treatment}/{stage}.done`

## Testing Changes

```bash
# Test single sample manually (bypass SLURM)
SLURM_ARRAY_TASK_ID=0 bash ./XX_script.sh

# Test with specific treatment
TREATMENT="test_treatment" bash ./XX_script.sh

# Force rerun (ignore checkpoints)
FORCE_RUN=true bash ./XX_script.sh
```

## Important Gotchas

1. **Conda + set -u**: Conda activation scripts have unbound variables. Always use `activate_env`/`deactivate_env` functions.

2. **Return values in pipes**: Use `>&2` for logging in functions that return values:
   ```bash
   func() {
       echo "status" >&2  # Goes to stderr (for logging)
       echo "result"      # Goes to stdout (return value)
   }
   ```

3. **Path length limits**: CheckM2 temp directories must use short paths (e.g., `/tmp/checkm2_*`) to avoid hitting filesystem limits.

4. **SLURM array size**: Calculate from sample/treatment counts, not hardcoded.

5. **Multi-lane samples**: Some samples span multiple sequencing lanes/runs. The `-01_merge_lanes.sh` script handles detection and merging.

## Dependencies

### Conda Environments Required
- `metagenome_assembly` - SPAdes, MEGAHIT
- `checkm` - CheckM2, CoverM
- `metabat`, `maxbin`, `concoct` - Binners
- `dastool` - Bin consolidation
- `magpurify` - Contamination removal
- `bbtools` - BBMap, repair.sh
- `metawrap-env` - MetaWRAP
- `gtdbtk` - Taxonomic classification
- Optional: `eukfinder`, `binette`, `binspreader`, `gunc`, `comebin`, `semibin`

### External Databases
- MAGpurify database
- EukFinder databases (Centrifuge + PLAST)
- GTDB-Tk database

## Commit Message Guidelines

Follow conventional commits style:
- `Fix:` for bug fixes
- `Add:` for new features
- `Update:` for changes to existing functionality
- `Docs:` for documentation changes

Example: `Fix: Add deactivate_env function to handle conda deactivate with set -u`

## Quick Reference: Key Functions

From `00_config_utilities.sh`:

| Function | Purpose |
|----------|---------|
| `log "message"` | Log with timestamp |
| `activate_env "name"` | Safely activate conda environment |
| `deactivate_env` | Safely deactivate conda environment |
| `init_sample_storage` | Initialize sample tracking |
| `add_sample_info` | Add sample to tracking |
| `get_sample_info_by_index N` | Get Nth sample info (pipe-delimited) |
| `get_treatments` | List all treatments |
| `get_samples_for_treatment "name"` | List samples in treatment |
| `check_sample_checkpoint "sample" "stage"` | Check if stage completed |
| `create_sample_checkpoint "sample" "stage"` | Mark stage as completed |
| `setup_temp_dir` | Create temp directory |
| `cleanup_temp_dir "/path"` | Remove temp directory |
| `validate_read_counts "r1" "r2"` | Verify paired read counts match |

## Additional Documentation

- `README.md` - Comprehensive user guide with workflows and troubleshooting
- `SETUP_GUIDE.md` - EukFinder database setup
- `EUKFINDER_WORKFLOW.md` - EukFinder analysis workflows
- `ASSEMBLY_SUCCESS_RATE.md` - Assembly metrics calculation
