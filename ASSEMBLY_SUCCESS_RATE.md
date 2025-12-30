# Assembly Success Rate Calculation

## Overview

The assembly success rate measures the percentage of input reads that successfully mapped back to the assembled contigs. This metric provides important quality information about your metaSPAdes assemblies.

## How It Works

1. **Reads Mapping**: Uses BBMap to align the input reads back to the assembled contigs
2. **Success Calculation**: Calculates the percentage of reads that successfully mapped
3. **Output**: Saves the rate to `assembly_success_rate.txt` and updates assembly statistics

## Automatic Calculation (New Assemblies)

For new assemblies, the success rate is **automatically calculated** after assembly completes in:
- `02_assembly.sh` - Individual sample assemblies
- `02_coassembly.sh` - Co-assemblies per treatment

The rate is included in the log output and saved to the assembly directory.

## Manual Calculation (Existing Assemblies)

If you've already completed your pipeline and want to calculate success rates retroactively, use the standalone script.

### Required Environment Variables

Before running the script, you need to set two environment variables:

1. **`OUTPUT_DIR`**: Path to your pipeline output directory (where assemblies are located)
2. **`PIPELINE_DIR`**: Path to the pipeline scripts directory (where `00_config_utilities.sh` is located)

### Basic Usage

```bash
# Set required environment variables
export OUTPUT_DIR=/path/to/your/output/directory
export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline

# Calculate for all assemblies (both individual and coassembly)
./calculate_assembly_success_rates.sh

# Calculate only for individual assemblies
./calculate_assembly_success_rates.sh --mode individual

# Calculate only for coassemblies
./calculate_assembly_success_rates.sh --mode coassembly

# Calculate for specific treatment
./calculate_assembly_success_rates.sh --treatment treatment1

# Calculate for specific sample
./calculate_assembly_success_rates.sh --mode individual --sample sample1 --treatment treatment1
```

### Setting Variables Inline

```bash
# Set variables and run in one command
OUTPUT_DIR=/path/to/output PIPELINE_DIR=/path/to/scripts ./calculate_assembly_success_rates.sh
```

### Submit as SLURM Job

The script supports two execution modes when submitted via SLURM:

#### **Array Job Mode (Recommended - Parallel Processing)**

Each treatment runs as a separate job for maximum parallelism:

```bash
# Set environment variables first
export OUTPUT_DIR=/path/to/your/output/directory
export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline

# Submit array job (processes each treatment in parallel)
sbatch calculate_assembly_success_rates.sh

# The script automatically uses array job mode and will:
# - Detect all treatments
# - Launch one job per treatment
# - Process them in parallel (up to 10 simultaneous jobs by default)
```

#### **Sequential Mode (Single Job)**

Process all treatments sequentially in one job:

```bash
# Disable array job by removing array directive
sbatch --array= calculate_assembly_success_rates.sh

# Or run directly (not as batch job)
./calculate_assembly_success_rates.sh
```

#### **For Specific Mode or Treatment**

```bash
# Array job for specific mode
sbatch --export=ALL calculate_assembly_success_rates.sh --mode individual

# Sequential for specific treatment
sbatch --array= --export=ALL calculate_assembly_success_rates.sh --treatment treatment1
```

#### **Setting Variables Inline**

```bash
# Array job mode with inline variables
sbatch --export=OUTPUT_DIR=/path/to/output,PIPELINE_DIR=/path/to/scripts calculate_assembly_success_rates.sh

# Sequential mode with inline variables
sbatch --array= --export=OUTPUT_DIR=/path/to/output,PIPELINE_DIR=/path/to/scripts calculate_assembly_success_rates.sh
```

### Why Use Array Job Mode?

**Benefits:**
- **Much faster**: All treatments processed in parallel instead of sequentially
- **Better resource utilization**: Each treatment gets dedicated CPUs and memory
- **Fault tolerant**: If one treatment fails, others continue
- **Easy monitoring**: Use `squeue` to see progress of each treatment

**When to use:**
- You have multiple treatments (recommended for 2+ treatments)
- You want results quickly
- Your cluster allows array jobs

**When to use sequential mode:**
- You have only one treatment
- Your cluster restricts array jobs
- You prefer simpler job management

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-m, --mode MODE` | Assembly mode: `individual`, `coassembly`, or `both` | `both` |
| `-s, --sample NAME` | Process specific sample (individual mode only) | All samples |
| `-t, --treatment NAME` | Process specific treatment | All treatments |
| `-o, --output FILE` | Output summary report filename | `assembly_success_rates_summary.tsv` |
| `-h, --help` | Show help message | - |

**Note:** In array job mode, the `--treatment` option is ignored since each array task processes one treatment automatically.

## Output Files

### Per-Assembly Output

For each assembly, the script creates/updates:

1. **`assembly_success_rate.txt`**
   - Contains just the percentage value (e.g., `85.23`)
   - Located in the assembly directory

2. **`assembly_statistics.txt`** (updated)
   - Adds an "Assembly Quality" section with the success rate
   - Includes explanation of the metric

3. **Mapping files** (in assembly directory)
   - `assembly_mapping_stats.txt` - BBMap statistics
   - `assembly_mapping_scafstats.txt` - Per-scaffold statistics
   - `assembly_mapping_covstats.txt` - Coverage statistics
   - `assembly_mapping_rpkm.txt` - RPKM values

### Summary Report

The standalone script generates a TSV summary file (default: `assembly_success_rates_summary.tsv`) with columns:

- **Type**: `individual` or `coassembly`
- **Treatment**: Treatment name
- **Sample**: Sample name (or `ALL` for coassembly)
- **Total_Reads**: Number of input reads
- **Assembly_Success_Rate(%)**: Calculated success rate
- **Contigs_File**: Path to contigs file
- **Status**: Calculation status
  - `success` - Newly calculated
  - `already_calculated` - Previously calculated (skipped)
  - `failed` - Calculation failed
  - `no_assembly` - Assembly file not found
  - `no_reads` - Input reads not found

## Example Summary Output

```tsv
Type        Treatment  Sample    Total_Reads  Assembly_Success_Rate(%)  Contigs_File                              Status
individual  control    sample1   5000000      87.45                    /path/to/assembly/control/sample1/...     success
individual  control    sample2   4800000      82.31                    /path/to/assembly/control/sample2/...     success
coassembly  control    ALL       12000000     89.67                    /path/to/coassembly/control/...           success
```

## Interpreting Results

### What is a "good" assembly success rate?

- **>80%**: Excellent assembly, most reads incorporated
- **60-80%**: Good assembly, typical for complex metagenomes
- **40-60%**: Moderate assembly, may indicate high diversity or low coverage
- **<40%**: Poor assembly, check input quality and assembly parameters

### Factors affecting success rate:

1. **Community complexity**: More diverse communities → lower rates
2. **Sequencing depth**: Higher coverage → higher rates
3. **Read quality**: Better quality → higher rates
4. **Assembly parameters**: Different k-mer sizes affect assembly
5. **Strain heterogeneity**: High variation → lower rates

## Troubleshooting

### Script can't find 00_config_utilities.sh

**Error message:**
```
ERROR: Cannot find 00_config_utilities.sh
```

**Solution:**
Set the `PIPELINE_DIR` environment variable to point to the directory containing the pipeline scripts:

```bash
export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline
```

If submitting via SLURM, use:
```bash
sbatch --export=PIPELINE_DIR=/path/to/scripts,OUTPUT_DIR=/path/to/output calculate_assembly_success_rates.sh
```

### Script can't find input reads

The script looks for reads in this order:
1. `validated/` directory (after validation stage)
2. `quality_filtering/` directory (after quality filtering stage)

If neither exists, the calculation will fail.

### Script can't find merged reads (coassembly)

For coassembly, the script requires merged reads saved in:
```
OUTPUT_DIR/coassembly/{treatment}/merged_reads/
  ├── merged_R1.fastq.gz
  └── merged_R2.fastq.gz
```

These are automatically saved by `02_coassembly.sh` (line 263-276).

### Calculation is slow

The script uses BBMap with fast mode enabled. For large assemblies:
- Uses 16 CPUs by default (adjust `#SBATCH --cpus-per-task`)
- Requires ~64GB RAM (adjust `#SBATCH --mem` if needed)
- Time scales with read count and assembly size

### Re-calculating existing rates

The script automatically skips assemblies that already have `assembly_success_rate.txt`. To force recalculation:

```bash
# Remove existing success rate files
find ${OUTPUT_DIR}/assembly -name "assembly_success_rate.txt" -delete
find ${OUTPUT_DIR}/coassembly -name "assembly_success_rate.txt" -delete

# Then run the script again
./calculate_assembly_success_rates.sh
```

## Integration with Pipeline

The assembly success rate is now part of the standard pipeline workflow:

1. **Stage 2 (Assembly)**: Automatically calculated for individual assemblies
2. **Stage 2b (Coassembly)**: Automatically calculated for coassemblies
3. **Standalone Script**: Can calculate for completed assemblies

All three methods use the same underlying function (`calculate_assembly_success_rate()` in `00_config_utilities.sh`).

## Requirements

- **BBMap** conda environment
- Input reads must be available in validated/ or quality_filtering/ directories
- For coassembly: merged reads must exist in coassembly/{treatment}/merged_reads/

## Questions or Issues?

If you encounter problems or have questions about assembly success rates, please open an issue on the GitHub repository.
