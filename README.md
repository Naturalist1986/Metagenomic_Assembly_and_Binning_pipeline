# Metagenomic Assembly and Binning Pipeline

A comprehensive SLURM-based pipeline for metagenomic assembly, binning, and quality assessment. Supports both individual sample processing and treatment-level coassembly workflows.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Pipeline Stages](#pipeline-stages)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Sample-Level Binning](#sample-level-binning-mode)
  - [Treatment-Level Binning (Coassembly)](#treatment-level-binning-mode-coassembly)
- [Output Structure](#output-structure)
- [Workflow Diagrams](#workflow-diagrams)
- [Quality Metrics](#quality-metrics)
- [Troubleshooting](#troubleshooting)

---

## Overview

This pipeline processes metagenomic sequencing data through quality control, assembly, binning, refinement, reassembly, purification, and abundance quantification. It supports two distinct operational modes:

1. **Sample-Level Binning**: Each sample is processed independently from raw reads to final bins
2. **Treatment-Level Binning**: Samples within a treatment group are co-assembled and binned together

---

## Features

- **Dual Processing Modes**: Sample-level or treatment-level (coassembly) binning
- **Comprehensive Quality Control**: FastQC, MultiQC, read validation
- **Advanced Binning**: Multiple binners (MetaBAT2, MaxBin2, CONCOCT) with DAS Tool refinement
- **Bin Improvement Pipeline**:
  - Reassembly with MetaWRAP (generates orig/strict/permissive versions)
  - Contamination removal with MAGpurify
  - Quality assessment with CheckM2
  - Intelligent bin version selection
- **Dual Quantification**:
  - MetaWRAP quant_bins for original refined bins
  - CoverM for selected high-quality bins
- **SLURM Integration**: Efficient array job scheduling
- **Checkpoint System**: Resume failed jobs without reprocessing

---

## Prerequisites

### Required Software

The pipeline requires the following conda environments:

- **metagenome_assembly**: megahit, spades, metawrap
- **checkm**: checkm2, coverm
- **metabat**: metabat2
- **maxbin**: maxbin2
- **concoct**: concoct
- **dastool**: DAS Tool
- **magpurify**: MAGpurify

### Required Files

1. **Sample Information File** (`samples.csv`):
   ```csv
   sample_name,treatment,read1,read2
   sample1,treatment1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
   sample2,treatment1,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
   sample3,treatment2,/path/to/sample3_R1.fastq.gz,/path/to/sample3_R2.fastq.gz
   ```

2. **Treatments File** (`treatments.txt`) - for treatment-level binning:
   ```
   treatment1
   treatment2
   treatment3
   ```

3. **Configuration File** (`00_config_utilities.sh`):
   - Set paths to data and conda environments
   - Configure SLURM parameters
   - Set quality thresholds

---

## Pipeline Stages

### Pre-Processing Stages (-1 to 04)

| Stage | Script | Description | Mode |
|-------|--------|-------------|------|
| -1 | `run_fastqc.sh` | Initial quality control with FastQC | Per sample |
| 00 | `00_quality_filtering.sh` | Quality filtering and adapter trimming | Per sample |
| 00b | `00b_validate_reads.sh` | Read validation and error correction | Per sample |
| 01 | `01_assembly.sh` | Individual sample assembly | Per sample |
| 01b | `01b_coassembly.sh` | Treatment-level coassembly | Per treatment |
| 02 | `02_binning.sh` | Multiple binning algorithms | Depends on mode |
| 03 | `03_bin_consolidation.sh` | Consolidate binner outputs | Depends on mode |
| 04 | `04_bin_refinement.sh` | DAS Tool refinement | Depends on mode |

### Post-Refinement Stages (05 to 08)

| Stage | Script | Description | Mode |
|-------|--------|-------------|------|
| 05 | `05_bin_reassembly.sh` | MetaWRAP reassembly | Depends on mode |
| 06 | `06_magpurify.sh` | Contamination removal | Depends on mode |
| 07 | `07_checkm2.sh` | Quality assessment | Depends on mode |
| 07b | `07b_bin_selection.sh` | Select best bin versions | Depends on mode |
| 08a | `08a_metawrap_quant.sh` | Quantify refined bins | Depends on mode |
| 08 | `08_coverm.sh` | Abundance calculation | Per sample |

---

## Configuration

### Key Configuration Variables

Edit `00_config_utilities.sh` to set:

```bash
# Paths
OUTPUT_DIR="/path/to/output"
SAMPLES_FILE="/path/to/samples.csv"
TREATMENTS_FILE="/path/to/treatments.txt"

# Conda environments
CONDA_BASE="/path/to/miniconda3"

# Mode selection
TREATMENT_LEVEL_BINNING=false  # Set to 'true' for coassembly mode

# Quality thresholds
MIN_READ_LENGTH=50
MIN_CONTIG_LENGTH=1000
MIN_BIN_SIZE=500000

# SLURM settings
DEFAULT_PARTITION="standard"
DEFAULT_QOS="normal"
```

---

## Usage

## Sample-Level Binning Mode

**Best for**: Samples with sufficient sequencing depth individually, no biological replication

### Configuration

```bash
# In 00_config_utilities.sh
TREATMENT_LEVEL_BINNING=false
```

### Workflow

```bash
# 1. Quality Control
sbatch -1_run_fastqc.sh
sbatch --array=0-N%10 00_quality_filtering.sh
sbatch --array=0-N%10 00b_validate_reads.sh

# 2. Assembly and Binning (per sample)
sbatch --array=0-N%10 01_assembly.sh
sbatch --array=0-N%10 02_binning.sh
sbatch --array=0-N%10 03_bin_consolidation.sh
sbatch --array=0-N%10 04_bin_refinement.sh

# 3. Bin Improvement (per sample)
sbatch --array=0-N%10 05_bin_reassembly.sh
sbatch --array=0-N%10 06_magpurify.sh
sbatch --array=0-N%10 07_checkm2.sh
sbatch --array=0-N%10 07b_bin_selection.sh

# 4. Quantification (per sample)
sbatch --array=0-N%10 08a_metawrap_quant.sh
sbatch --array=0-N%10 08_coverm.sh
```

Where `N` is the number of samples minus 1.

### Output Structure

```
output/
├── assembly/
│   └── treatment/
│       └── sample/
│           └── contigs.fasta
├── bin_refinement/
│   └── treatment/
│       └── sample/
│           └── dastool_DASTool_bins/*.fa
├── reassembly/
│   └── treatment/
│       └── sample/
│           └── reassembled_bins/
│               ├── bin.orig.fa
│               ├── bin.strict.fa
│               └── bin.permissive.fa
├── magpurify/
│   └── treatment/
│       └── sample/
│           └── purified_bins/*.fa
├── checkm2/
│   └── treatment/
│       └── sample/
│           └── quality_report.tsv
├── selected_bins/
│   └── treatment/
│       └── sample/
│           └── *.fa  (best version of each bin)
└── coverm/
    └── treatment/
        └── sample/
            └── abundance.tsv
```

---

## Treatment-Level Binning Mode (Coassembly)

**Best for**: Multiple samples per treatment, low individual depth, biological replicates

### Configuration

```bash
# In 00_config_utilities.sh
TREATMENT_LEVEL_BINNING=true
```

### Workflow

```bash
# 1. Quality Control (per sample)
sbatch -1_run_fastqc.sh
sbatch --array=0-N%10 00_quality_filtering.sh
sbatch --array=0-N%10 00b_validate_reads.sh

# 2. Coassembly and Binning (per treatment)
sbatch --array=0-M%5 01b_coassembly.sh        # M = num treatments - 1
sbatch --array=0-M%5 02_binning.sh
sbatch --array=0-M%5 03_bin_consolidation.sh
sbatch --array=0-M%5 04_bin_refinement.sh

# 3. Bin Improvement (ONCE per treatment)
sbatch --array=0-M%5 05_bin_reassembly.sh
sbatch --array=0-M%5 06_magpurify.sh
sbatch --array=0-M%5 07_checkm2.sh
sbatch --array=0-M%5 07b_bin_selection.sh

# 4. Quantification
sbatch --array=0-M%5 08a_metawrap_quant.sh    # Once per treatment
sbatch --array=0-N%10 08_coverm.sh             # Per sample (dual mode)
```

### Key Differences in Treatment-Level Mode

1. **Coassembly**: All samples in a treatment are assembled together
   - Merged reads saved in `coassembly/${treatment}/merged_reads/`
   - Single assembly: `coassembly/${treatment}/contigs.fasta`

2. **Binning**: Performed once per treatment
   - Bins represent the collective microbial community
   - Higher depth improves binning quality

3. **Reassembly**: Uses merged treatment reads
   - MetaWRAP maps ALL treatment reads to bins
   - Produces orig/strict/permissive versions

4. **Quality Control**: Run once on all bin versions
   - MAGpurify processes all bins together
   - CheckM2 evaluates all versions simultaneously

5. **Bin Selection**: Chooses best version per bin
   - Compares orig/strict/permissive using CheckM2 scores
   - Quality score = Completeness - (5 × Contamination)

6. **Dual Mode Quantification**:
   - **Stage 08a**: Quantifies original refined bins using merged reads
   - **Stage 08**: Each sample quantifies the same selected bins
   - Enables cross-sample comparison with consistent bin set

### Output Structure (Treatment-Level)

```
output/
├── coassembly/
│   └── treatment/
│       ├── contigs.fasta
│       └── merged_reads/
│           ├── merged_R1.fastq.gz
│           ├── merged_R2.fastq.gz
│           └── merged_singletons.fastq.gz
├── bin_refinement/
│   └── treatment/
│       └── dastool_DASTool_bins/*.fa
├── reassembly/
│   └── treatment/
│       └── reassembled_bins/
│           ├── bin.orig.fa
│           ├── bin.strict.fa
│           └── bin.permissive.fa
├── magpurify/
│   └── treatment/
│       └── purified_bins/*.fa
├── checkm2/
│   └── treatment/
│       └── quality_report.tsv
├── selected_bins/
│   └── treatment/
│       └── *.fa  (best version: orig, strict, or permissive)
├── metawrap_quant/
│   └── treatment/
│       └── bin_abundance_table.tab
└── coverm/
    └── treatment/
        └── sample1/
        │   └── abundance.tsv  (same bins across all samples)
        └── sample2/
            └── abundance.tsv
```

---

## Workflow Diagrams

### Sample-Level Workflow

```
Raw Reads (Sample 1, 2, 3, ...)
    ↓
Quality Filtering (per sample)
    ↓
Assembly (per sample)
    ↓
Binning (per sample)
    ↓
Refinement (per sample)
    ↓
Reassembly (per sample)
    ↓
MAGpurify (per sample)
    ↓
CheckM2 (per sample)
    ↓
Bin Selection (per sample)
    ↓
Quantification (per sample)
```

### Treatment-Level Workflow

```
Raw Reads (Samples in Treatment)
    ↓
Quality Filtering (per sample)
    ↓
Merge Reads (per treatment) ─────────────┐
    ↓                                      │
Coassembly (per treatment)                │
    ↓                                      │
Binning (per treatment)                   │
    ↓                                      │
Refinement (per treatment)                │
    ↓                                      │
Reassembly (ONCE per treatment) ←─────────┘
    ↓ (uses merged reads)
MAGpurify (ONCE per treatment)
    ↓ (all bin versions)
CheckM2 (ONCE per treatment)
    ↓ (quality assessment)
Bin Selection (ONCE per treatment)
    ↓ (selects best versions)
    ├─→ MetaWRAP Quant (merged reads)
    └─→ CoverM (per sample, same bins)
```

---

## Quality Metrics

### Bin Quality Categories (MIMAG Standards)

| Category | Completeness | Contamination |
|----------|--------------|---------------|
| High Quality | ≥90% | ≤5% |
| Medium Quality | ≥50% | ≤10% |
| Low Quality | <50% | >10% |

### Quality Score Formula

```
Quality Score = Completeness - (5 × Contamination)
```

Used in stage 07b to select the best bin version.

### CheckM2 Output Interpretation

- **Completeness**: Percentage of expected genes present
- **Contamination**: Percentage of duplicate/unexpected genes
- **Strain Heterogeneity**: Genetic diversity within bin

---

## Advanced Features

### Checkpoint System

The pipeline automatically tracks completed stages:

```bash
# Check if sample completed a stage
check_sample_checkpoint "$SAMPLE_NAME" "assembly"

# Check if treatment completed a stage
check_treatment_checkpoint "$TREATMENT" "reassembly"
```

Checkpoints are stored in:
- `${OUTPUT_DIR}/checkpoints/${treatment}/${sample}/${stage}.done`
- `${OUTPUT_DIR}/checkpoints/${treatment}/${stage}.done`

### Resume Failed Jobs

If a job fails, simply resubmit:

```bash
# Automatically skips completed samples/treatments
sbatch --array=0-N%10 05_bin_reassembly.sh
```

### Bin Version Selection Logic

Stage 07b evaluates all three reassembly versions:

1. **orig**: Original reassembly (same as no suffix)
2. **strict**: More stringent contig filtering
3. **permissive**: More lenient contig inclusion

Selection criteria:
- Calculate quality score for each version
- Select version with highest score
- Prioritize completeness while penalizing contamination

---

## Output Files

### Key Output Files

| File | Description |
|------|-------------|
| `quality_report.tsv` | CheckM2 quality metrics for all bins |
| `bin_selection_report.txt` | Shows which version was selected per bin |
| `abundance.tsv` | CoverM abundance metrics |
| `bin_abundance_table.tab` | MetaWRAP quantification results |
| `abundance_enhanced.tsv` | Combined abundance + quality data |

### Important Intermediate Files

| File | Description |
|------|-------------|
| `contigs.fasta` | Assembly contigs |
| `dastool_DASTool_bins/*.fa` | Refined bins from DAS Tool |
| `merged_reads/*.fastq.gz` | Merged treatment reads (coassembly) |
| `reassembled_bins/*.fa` | Reassembled bin versions |
| `purified_bins/*.fa` | MAGpurify cleaned bins |
| `selected_bins/*.fa` | Final selected high-quality bins |

---

## Troubleshooting

### Common Issues

#### 1. No bins found for treatment-level processing

**Symptom**: Error message: "No refined bins found for treatment X"

**Solution**: Ensure `TREATMENT_LEVEL_BINNING=true` in config and treatments file exists

#### 2. Merged reads not found

**Symptom**: Stage 05 can't find merged reads for coassembly

**Solution**:
- Rerun stage 01b (coassembly)
- Or stage 05 will create merged reads on-demand

#### 3. Array job index out of range

**Symptom**: "No sample/treatment found for array index"

**Solution**: Adjust array range based on number of samples/treatments:
```bash
# For 10 samples: --array=0-9%10
# For 5 treatments: --array=0-4%5
```

#### 4. CheckM2 quality report empty

**Symptom**: Stage 07b fails because quality_report.tsv is empty

**Solution**:
- Check that stage 06 (MAGpurify) completed successfully
- Ensure bins meet minimum size thresholds (500KB, 5 contigs)
- Review CheckM2 logs for errors

#### 5. CoverM not finding selected_bins

**Symptom**: Stage 08 falls back to older bin sources

**Solution**:
- Verify stage 07b completed successfully
- Check that `selected_bins/${treatment}/` directory exists and contains .fa files

### Memory Requirements

If jobs fail with out-of-memory errors:

| Stage | Recommended Memory |
|-------|-------------------|
| Assembly | 64-128 GB |
| Binning | 32-64 GB |
| Reassembly | 64-128 GB |
| MAGpurify | 32-64 GB |
| CheckM2 | 64-128 GB |
| CoverM | 32-64 GB |

Adjust in script headers:
```bash
#SBATCH --mem=128G
```

### Time Limits

Adjust based on dataset size:

| Stage | Typical Time |
|-------|--------------|
| Assembly | 4-12 hours |
| Coassembly | 12-24 hours |
| Binning | 4-8 hours |
| Reassembly | 8-16 hours |
| CheckM2 | 4-8 hours |

---

## Best Practices

### Sample-Level Mode

- **Use when**: Each sample has >5 GB sequencing data
- **Advantages**:
  - Sample-specific bins
  - Faster processing per sample
  - Better for samples with different communities

### Treatment-Level Mode

- **Use when**:
  - Samples have <5 GB individual depth
  - Biological replicates exist
  - Community composition similar within treatment
- **Advantages**:
  - Improved binning from combined depth
  - Consistent bins across samples for comparison
  - Better for low-biomass samples

### General Recommendations

1. **Quality filtering**: Don't skip stages -1, 00, 00b
2. **Reassembly versions**: Always run stage 07b to select best versions
3. **Dual quantification**: Use both 08a (merged) and 08 (per-sample) for full picture
4. **Checkpoints**: Don't delete checkpoint files unless reprocessing from scratch
5. **Log files**: Check logs in `${LOG_DIR}/${treatment}/` for troubleshooting

---

## Citation

If you use this pipeline, please cite the following tools:

- **MetaWRAP**: Uritskiy et al. (2018) Microbiome
- **DAS Tool**: Sieber et al. (2018) Nature Microbiology
- **CheckM2**: Chklovski et al. (2023) Nature Methods
- **MAGpurify**: Nayfach et al. (2019) Bioinformatics
- **CoverM**: https://github.com/wwood/CoverM
- **MEGAHIT**: Li et al. (2015) Bioinformatics
- **MetaBAT2**: Kang et al. (2019) PeerJ
- **MaxBin2**: Wu et al. (2016) Bioinformatics
- **CONCOCT**: Alneberg et al. (2014) Nature Methods

---

## Support

For issues, questions, or contributions, please open an issue on the GitHub repository.

---

## License

This pipeline is provided as-is for academic and research use.

---

**Last Updated**: 2025-01-XX
**Pipeline Version**: 2.0
