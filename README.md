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
| -1 | `-01_merge_lanes.sh` | Lane/Run detection and merging | Per sample |
| 00 | `00_quality_filtering.sh` | Quality filtering and adapter trimming | Per sample |
| 0.5 | `00b_validate_repair.sh` | Read validation and repair | Per sample |
| 01 | `01_assembly.sh` | Individual sample assembly | Per sample |
| 01b | `01b_coassembly.sh` | Treatment-level coassembly | Per treatment |
| 02 | `02_plasmid_detection.sh` | Plasmid detection (PlasClass, MOB-suite) | Per sample |
| 03 | `03_binning.sh` | Multiple binning algorithms | Depends on mode |
| 04 | `04_bin_refinement.sh` | DAS Tool refinement | Depends on mode |

### Post-Refinement Stages (05 to 10)

| Stage | Script | Description | Mode |
|-------|--------|-------------|------|
| 05 | `05_bin_reassembly.sh` | MetaWRAP reassembly | Depends on mode |
| 06 | `06_magpurify.sh` | Contamination removal | Depends on mode |
| 07 | `07_checkm2.sh` | Quality assessment | Depends on mode |
| 7.5 | `07b_bin_selection.sh` | Select best bin versions | Depends on mode |
| 08 | `08a_metawrap_quant.sh` | Quantify refined bins | Depends on mode |
| 09 | `08b_bin_collection.sh` | Bin collection, consolidate abundance, GTDB-Tk | Per treatment |
| 10 | `09_final_report.sh` | Generate taxonomy-labeled abundance plots | Single execution |

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
# 1. Lane/Run Merge and Quality Control
sbatch --array=0-N%10 -01_merge_lanes.sh
sbatch --array=0-N%10 00_quality_filtering.sh
sbatch --array=0-N%10 00b_validate_repair.sh

# 2. Assembly, Plasmid Detection, and Binning (per sample)
sbatch --array=0-N%10 01_assembly.sh
sbatch --array=0-N%10 02_plasmid_detection.sh
sbatch --array=0-N%10 03_binning.sh
sbatch --array=0-N%10 04_bin_refinement.sh

# 3. Bin Improvement (per sample)
sbatch --array=0-N%10 05_bin_reassembly.sh
sbatch --array=0-N%10 06_magpurify.sh
sbatch --array=0-N%10 07_checkm2.sh
sbatch --array=0-N%10 07b_bin_selection.sh

# 4. Quantification and Reporting
sbatch --array=0-N%10 08a_metawrap_quant.sh
sbatch --array=0-N%10 08b_bin_collection.sh  # Per treatment
sbatch 09_final_report.sh  # Single execution for all treatments
```

Where `N` is the number of samples minus 1.

### Output Structure

```
output/
├── merged_reads/          # Stage -1: Merged lanes/runs
│   └── treatment/
│       └── sample/
│           ├── merged_R1.fastq.gz
│           └── merged_R2.fastq.gz
├── assembly/             # Stage 1: Assembly
│   └── treatment/
│       └── sample/
│           └── contigs.fasta
├── plasmid_detection/    # Stage 2: Plasmid detection
│   └── treatment/
│       └── sample/
│           └── plasmid_results/
├── bin_refinement/       # Stage 4: Refined bins
│   └── treatment/
│       └── sample/
│           └── dastool_DASTool_bins/*.fa
├── reassembly/           # Stage 5: Reassembled bins
│   └── treatment/
│       └── sample/
│           └── reassembled_bins/
│               ├── bin.orig.fa
│               ├── bin.strict.fa
│               └── bin.permissive.fa
├── magpurify/            # Stage 6: Purified bins
│   └── treatment/
│       └── sample/
│           └── purified_bins/*.fa
├── checkm2/              # Stage 7: Quality assessment
│   └── treatment/
│       └── sample/
│           └── quality_report.tsv
├── selected_bins/        # Stage 7.5: Best bin versions
│   └── treatment/
│       └── sample/
│           └── *.fa
├── metawrap_quant/       # Stage 8: Quantification
│   └── treatment/
│       └── sample/
│           └── bin_abundance_table.tab
├── bin_collection/       # Stage 9: Per-treatment collection
│   └── treatment/
│       ├── all_bins/     # All selected bins
│       ├── gtdbtk/       # GTDB-Tk results
│       └── coverm_abundance_consolidated.tsv
└── final_report/         # Stage 10: Final visualizations
    ├── taxonomy_abundance_plots/
    └── summary_report.html
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
# 1. Lane/Run Merge and Quality Control (per sample)
sbatch --array=0-N%10 -01_merge_lanes.sh
sbatch --array=0-N%10 00_quality_filtering.sh
sbatch --array=0-N%10 00b_validate_repair.sh

# 2. Coassembly, Plasmid Detection, and Binning (per treatment)
sbatch --array=0-M%5 01b_coassembly.sh        # M = num treatments - 1
sbatch --array=0-N%10 02_plasmid_detection.sh  # Per sample
sbatch --array=0-M%5 03_binning.sh            # Per treatment
sbatch --array=0-M%5 04_bin_refinement.sh

# 3. Bin Improvement (ONCE per treatment)
sbatch --array=0-M%5 05_bin_reassembly.sh
sbatch --array=0-M%5 06_magpurify.sh
sbatch --array=0-M%5 07_checkm2.sh
sbatch --array=0-M%5 07b_bin_selection.sh

# 4. Quantification and Reporting
sbatch --array=0-M%5 08a_metawrap_quant.sh    # Once per treatment
sbatch --array=0-M%5 08b_bin_collection.sh    # Per treatment: consolidate abundance, GTDB-Tk
sbatch 09_final_report.sh                     # Single execution for all treatments
```

### Key Differences in Treatment-Level Mode

1. **Lane/Run Merge**: Automatically detects and merges multiple sequencing runs
   - Handles both multi-lane and multi-run data
   - Supports complex sample sheet formats

2. **Coassembly**: All samples in a treatment are assembled together
   - Merged reads saved in `coassembly/${treatment}/merged_reads/`
   - Single assembly: `coassembly/${treatment}/contigs.fasta`

3. **Plasmid Detection**: Performed per sample (stage 2)
   - Uses PlasClass and MOB-suite
   - Identifies plasmid contigs and mobile elements

4. **Binning**: Performed once per treatment (stage 3)
   - Bins represent the collective microbial community
   - Higher depth improves binning quality

5. **Reassembly**: Uses merged treatment reads (stage 5)
   - MetaWRAP maps ALL treatment reads to bins
   - Produces orig/strict/permissive versions

6. **Quality Control**: Run once on all bin versions (stages 6-7)
   - MAGpurify processes all bins together
   - CheckM2 evaluates all versions simultaneously

7. **Bin Selection**: Chooses best version per bin (stage 7.5)
   - Compares orig/strict/permissive using CheckM2 scores
   - Quality score = Completeness - (5 × Contamination)

8. **Bin Collection**: Per treatment consolidation (stage 9)
   - Consolidates CoverM abundance across all samples
   - Runs GTDB-Tk for taxonomic classification
   - Creates taxonomy-enhanced abundance tables

9. **Final Report**: Single execution for entire pipeline (stage 10)
   - Generates taxonomy-labeled abundance plots
   - Aggregates data from all treatments
   - Produces comprehensive visualization suite

### Output Structure (Treatment-Level)

```
output/
├── merged_reads/          # Stage -1: Per-sample merged lanes/runs
│   └── treatment/
│       └── sample/
│           ├── merged_R1.fastq.gz
│           └── merged_R2.fastq.gz
├── coassembly/           # Stage 1b: Co-assembly
│   └── treatment/
│       ├── contigs.fasta
│       └── merged_reads/
│           ├── merged_R1.fastq.gz
│           ├── merged_R2.fastq.gz
│           └── merged_singletons.fastq.gz
├── plasmid_detection/    # Stage 2: Per-sample plasmid detection
│   └── treatment/
│       └── sample/
│           └── plasmid_results/
├── bin_refinement/       # Stage 4: Treatment-level refined bins
│   └── treatment/
│       └── dastool_DASTool_bins/*.fa
├── reassembly/           # Stage 5: Treatment-level reassembly
│   └── treatment/
│       └── reassembled_bins/
│           ├── bin.orig.fa
│           ├── bin.strict.fa
│           └── bin.permissive.fa
├── magpurify/            # Stage 6: Treatment-level purification
│   └── treatment/
│       └── purified_bins/*.fa
├── checkm2/              # Stage 7: Treatment-level quality
│   └── treatment/
│       └── quality_report.tsv
├── selected_bins/        # Stage 7.5: Best bin versions
│   └── treatment/
│       └── *.fa  (best version: orig, strict, or permissive)
├── metawrap_quant/       # Stage 8: Treatment-level quantification
│   └── treatment/
│       └── bin_abundance_table.tab
├── bin_collection/       # Stage 9: Per-treatment collection
│   └── treatment/
│       ├── all_bins/     # All selected bins
│       ├── gtdbtk/       # GTDB-Tk taxonomic classification
│       │   ├── classify/
│       │   └── identify/
│       └── coverm_abundance_consolidated.tsv
└── final_report/         # Stage 10: Final visualizations
    ├── taxonomy_abundance_plots/
    │   ├── phylum_level.png
    │   ├── class_level.png
    │   ├── order_level.png
    │   ├── family_level.png
    │   └── genus_level.png
    └── summary_report.html
```

---

## Workflow Diagrams

### Sample-Level Workflow

```
Raw Reads (Sample 1, 2, 3, ...)
    ↓
Stage -1: Lane/Run Merge (per sample)
    ↓
Stage 0: Quality Filtering (per sample)
    ↓
Stage 0.5: Read Validation & Repair (per sample)
    ↓
Stage 1: Assembly (per sample)
    ↓
Stage 2: Plasmid Detection (per sample)
    ↓
Stage 3: Binning (per sample)
    ↓
Stage 4: Refinement (per sample)
    ↓
Stage 5: Reassembly (per sample)
    ↓
Stage 6: MAGpurify (per sample)
    ↓
Stage 7: CheckM2 (per sample)
    ↓
Stage 7.5: Bin Selection (per sample)
    ↓
Stage 8: MetaWRAP Quant (per sample)
    ↓
Stage 9: Bin Collection (per treatment)
    ├─→ Consolidate CoverM abundance
    └─→ GTDB-Tk taxonomic classification
    ↓
Stage 10: Final Report (single execution)
    └─→ Taxonomy-labeled abundance plots
```

### Treatment-Level Workflow

```
Raw Reads (Samples in Treatment)
    ↓
Stage -1: Lane/Run Merge (per sample)
    ↓
Stage 0: Quality Filtering (per sample)
    ↓
Stage 0.5: Read Validation & Repair (per sample)
    ↓
Stage 1b: Merge Reads & Coassembly (per treatment) ─────┐
    ↓                                                      │
Stage 2: Plasmid Detection (per sample)                  │
    ↓                                                      │
Stage 3: Binning (per treatment)                         │
    ↓                                                      │
Stage 4: Refinement (per treatment)                      │
    ↓                                                      │
Stage 5: Reassembly (ONCE per treatment) ←───────────────┘
    ↓ (uses merged treatment reads)
Stage 6: MAGpurify (ONCE per treatment)
    ↓ (all bin versions)
Stage 7: CheckM2 (ONCE per treatment)
    ↓ (quality assessment)
Stage 7.5: Bin Selection (ONCE per treatment)
    ↓ (selects best versions)
Stage 8: MetaWRAP Quant (per treatment)
    ↓
Stage 9: Bin Collection (per treatment)
    ├─→ Consolidate CoverM abundance
    └─→ GTDB-Tk taxonomic classification
    ↓
Stage 10: Final Report (single execution)
    └─→ Taxonomy-labeled abundance plots for all treatments
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
| `quality_report.tsv` | CheckM2 quality metrics for all bins (Stage 7) |
| `bin_selection_report.txt` | Shows which version was selected per bin (Stage 7.5) |
| `bin_abundance_table.tab` | MetaWRAP quantification results (Stage 8) |
| `coverm_abundance_consolidated.tsv` | Consolidated CoverM abundance across samples (Stage 9) |
| `gtdbtk.bac120.summary.tsv` | GTDB-Tk bacterial taxonomic classification (Stage 9) |
| `gtdbtk.ar53.summary.tsv` | GTDB-Tk archaeal taxonomic classification (Stage 9) |
| `taxonomy_abundance_plots/*.png` | Final taxonomy-labeled visualizations (Stage 10) |
| `summary_report.html` | Comprehensive pipeline summary (Stage 10) |

### Important Intermediate Files

| File | Description |
|------|-------------|
| `merged_R1/R2.fastq.gz` | Merged lanes/runs per sample (Stage -1) |
| `contigs.fasta` | Assembly contigs (Stage 1) |
| `plasmid_results/` | Plasmid detection outputs (Stage 2) |
| `dastool_DASTool_bins/*.fa` | Refined bins from DAS Tool (Stage 4) |
| `reassembled_bins/*.fa` | Reassembled bin versions (Stage 5) |
| `purified_bins/*.fa` | MAGpurify cleaned bins (Stage 6) |
| `selected_bins/*.fa` | Final selected high-quality bins (Stage 7.5) |
| `all_bins/` | Collected bins for GTDB-Tk (Stage 9) |

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

1. **Lane/Run merge**: Always run stage -1 if you have multi-lane or multi-run data
2. **Quality filtering**: Don't skip stages 0 and 0.5 (quality filtering and validation)
3. **Plasmid detection**: Stage 2 provides valuable information about mobile genetic elements
4. **Reassembly versions**: Always run stage 7.5 to select best bin versions
5. **Taxonomic classification**: Stage 9 (GTDB-Tk) is essential for biological interpretation
6. **Final visualization**: Stage 10 generates publication-ready plots with taxonomy labels
7. **Checkpoints**: Don't delete checkpoint files unless reprocessing from scratch
8. **Log files**: Check logs in `${OUTPUT_DIR}/logs/` for troubleshooting
9. **Master script**: Use `run_pipeline.sh` for automated execution of all stages

---

## Citation

If you use this pipeline, please cite the following tools:

- **MetaWRAP**: Uritskiy et al. (2018) Microbiome
- **DAS Tool**: Sieber et al. (2018) Nature Microbiology
- **CheckM2**: Chklovski et al. (2023) Nature Methods
- **MAGpurify**: Nayfach et al. (2019) Bioinformatics
- **GTDB-Tk**: Chaumeil et al. (2022) Bioinformatics
- **CoverM**: https://github.com/wwood/CoverM
- **PlasClass**: Pellow et al. (2020) NAR Genomics and Bioinformatics
- **MOB-suite**: Robertson & Nash (2018) Microbial Genomics
- **MEGAHIT**: Li et al. (2015) Bioinformatics
- **MetaSPAdes**: Nurk et al. (2017) Genome Research
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

**Last Updated**: 2025-11-02
**Pipeline Version**: 2.1

## Recent Updates

### Version 2.1 (2025-11-02)
- Added Stage -1: Lane/Run detection and merging
- Added Stage 2: Plasmid detection with PlasClass and MOB-suite
- Added Stage 9: Bin collection with GTDB-Tk taxonomic classification
- Added Stage 10: Final report generation with taxonomy-labeled plots
- Updated stage numbering to reflect actual pipeline architecture
- Enhanced support for multi-run sample sheets
- Improved workflow documentation and output structure descriptions
