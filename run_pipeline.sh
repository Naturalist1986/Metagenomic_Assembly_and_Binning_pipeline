#!/bin/bash

# run_pipeline.sh - Master script to run the complete metagenomic pipeline
# UPDATED: Now supports sample sheets with multiple runs per sample

# Default values
START_STAGE=0  # Start from quality filtering by default (skip optional merge/plasmid stages)
END_STAGE=9    # Updated: removed stage 8 (metawrap quant)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$SCRIPT_DIR"  # Export for SLURM jobs

SPECIFIC_TREATMENTS=()
SPECIFIC_SAMPLES=()
ASSEMBLY_MODE="individual"
TREATMENT_LEVEL_BINNING=false
USE_COMEBIN=false  # Use COMEBin for binning instead of MetaWRAP
SKIP_MERGE_LANES=true  # Skip lane merging by default (optional stage)
SKIP_PLASMID_DETECTION=true  # Skip plasmid detection by default (optional stage)

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Master script to run the complete metagenomic pipeline.
Supports multiple sample sheet formats including multi-run datasets.

OPTIONS:
    -s, --start-stage NUM         Start from stage NUM (0 to 9) [default: 0]
    -e, --end-stage NUM           End at stage NUM (0 to 9) [default: 9]
    -a, --assembly-mode MODE      Assembly mode: 'individual' or 'coassembly' [default: individual]
    -b, --treatment-level-binning Use treatment-level binning instead of sample-level
    --comebin                     Use COMEBin for binning instead of MetaWRAP (requires comebin conda env)
    -t, --treatment NAME          Run only for specific treatment/group (can be used multiple times)
    -m, --sample NAME             Run only for specific sample (can be used multiple times)
    -i, --input-dir PATH          Input directory containing FASTQ files
    -o, --output-dir PATH         Output directory for results
    -S, --sample-sheet PATH       Sample sheet (Excel, CSV, or TSV)
    -A, --account NAME            SLURM account name for job submission
    --merge-lanes                 Include lane/run merging stage (optional, off by default)
    --plasmid-detection           Include plasmid detection stage (optional, off by default)
    --assembly-threads NUM        Number of threads for assembly [default: SLURM_CPUS_PER_TASK or 32]
    --assembly-memory NUM         Memory in GB for individual assembly [default: 250]
    --coassembly-memory NUM       Memory in GB for coassembly [default: 500]
    --bin-reassembly-threads NUM  Threads for bin reassembly [default: 8]
    --bin-reassembly-memory NUM   Memory in GB for bin reassembly [default: 32]
    -c, --create-template         Create sample sheet template and exit
    -d, --dry-run                 Show what would be run without executing
    -h, --help                    Show this help message

SAMPLE SHEET FORMATS:
    The pipeline supports two sample sheet formats:

    FORMAT 1 - Legacy (TSV/Excel):
        Columns: Sample_Name, Treatment, R1_File, R2_File
        Example:
            Sample_Name    Treatment    R1_File              R2_File
            sample1        control      sample1_R1.fq.gz     sample1_R2.fq.gz
            sample2        treated      sample2_R1.fq.gz     sample2_R2.fq.gz

    FORMAT 2 - Multi-Run (CSV/TSV/Excel):
        Columns: sample, run, group, short_reads_1, short_reads_2, [short_reads_platform]
        Example:
            sample,run,group,short_reads_1,short_reads_2,short_reads_platform
            sample1,1,control,s1_run1_R1.fq.gz,s1_run1_R2.fq.gz,ILLUMINA
            sample1,2,control,s1_run2_R1.fq.gz,s1_run2_R2.fq.gz,ILLUMINA
            sample2,1,treated,s2_run1_R1.fq.gz,s2_run1_R2.fq.gz,ILLUMINA

        - Multiple runs of the same sample are automatically detected and merged in stage -1
        - 'group' column is treated as 'treatment' in the pipeline
        - 'run' column is used to group multiple sequencing runs
        - Format is auto-detected based on column headers

STAGES:
    0  - Quality Filtering (Trimmomatic)
    1  - Validation & Repair (validates paired-end sync, repairs if needed)
    2  - Assembly (MetaSPAdes - individual or co-assembly)
    3  - Binning (MetaWRAP or COMEBin - sample/treatment-level based on flags)
    4  - Bin Refinement (DAS Tool - sample-level or treatment-level based on -b flag)
    5  - Bin Reassembly (MetaWRAP - sample-level or treatment-level based on -b flag)
    6  - MAGpurify (contamination removal - sample-level or treatment-level based on -b flag)
    7  - CheckM2 (quality assessment - sample-level or treatment-level based on -b flag)
    8  - Bin Selection (selects best reassembly version: orig/strict/permissive)
    9  - Bin Collection (per treatment: collects selected bins, runs CoverM & GTDB-Tk)

OPTIONAL STAGES (use flags to enable):
    --merge-lanes        Lane/Run Detection & Merge (auto-detects and merges multiple lanes/runs)
    --plasmid-detection  Plasmid Detection (PlasClass and MOB-suite - runs after assembly)

ADDITIONAL TOOLS (run separately after pipeline completes):
    final_report.sh      Generate cross-treatment comparison plots with taxonomy labels
                         (Requires all treatments to complete stage 9 first)

BINNING MODES:
    By default, binning is performed at the sample level using MetaWRAP (MetaBat2, MaxBin2, Concoct).

    Use --treatment-level-binning to bin at the treatment level instead:
    - All samples in a treatment are binned together
    - Useful for co-assembly workflows or when you want bins across replicates
    - Uses 03_binning_treatment_level.sh instead of 03_binning.sh

    Use --comebin to use COMEBin for binning instead of MetaWRAP:
    - Uses contrastive multi-view representation learning
    - Requires 'comebin' conda environment
    - Uses BAM files from MetaWRAP work_files or creates them with bowtie2
    - Uses 03b_comebin.sh instead of standard binning scripts

EXAMPLES:
    # Run complete pipeline (default: stages 0-9, no lane merge or plasmid detection)
    $0 -i /path/to/fastq -o /path/to/output

    # Run with lane merging for multi-run samples
    $0 --merge-lanes -i /path/to/fastq -o /path/to/output

    # Run with plasmid detection
    $0 --plasmid-detection -i /path/to/fastq -o /path/to/output

    # Run with all optional stages
    $0 --merge-lanes --plasmid-detection -i /path/to/fastq -o /path/to/output

    # Using sample sheet with SLURM account
    $0 -S samples.csv -i /path/to/fastq -o /path/to/output -A my_account

    # Co-assembly with treatment-level binning
    $0 --assembly-mode coassembly --treatment-level-binning \\
       -i /path/to/fastq -o /path/to/output

    # Use COMEBin for binning instead of MetaWRAP
    $0 --comebin -i /path/to/fastq -o /path/to/output

    # Run specific treatment only
    $0 --treatment control -i /path/to/fastq -o /path/to/output

    # Run specific samples only
    $0 --sample sample1 --sample sample2 -i /path/to/fastq -o /path/to/output

    # Run only quality filtering and validation (stages 0-1)
    $0 --start-stage 0 --end-stage 1 -i /path/to/fastq -o /path/to/output

    # Dry run to preview execution
    $0 --dry-run -i /path/to/fastq -o /path/to/output

    # Custom assembly resources (16 threads, 128GB memory)
    $0 --assembly-threads 16 --assembly-memory 128 \\
       -i /path/to/fastq -o /path/to/output

    # Custom coassembly with different memory allocation
    $0 --assembly-mode coassembly --coassembly-memory 400 \\
       -i /path/to/fastq -o /path/to/output

    # Create sample sheet template
    $0 --create-template -i /path/to/fastq

NOTES:
    - Optional stages (lane merging, plasmid detection) are OFF by default
    - Use --merge-lanes if you have multiple lanes/runs per sample that need merging
    - Use --plasmid-detection to identify plasmid contigs and mobile elements
    - File paths in sample sheets should be relative to INPUT_DIR
    - The pipeline auto-detects sample sheet format from column headers
    - Set INPUT_DIR and OUTPUT_DIR or they will be prompted

EOF
}

# Parse arguments
DRY_RUN=false
CREATE_TEMPLATE=false
INPUT_DIR_ARG=""
OUTPUT_DIR_ARG=""
SAMPLE_SHEET_ARG=""
SLURM_ACCOUNT_ARG=""
ASSEMBLY_THREADS_ARG=""
ASSEMBLY_MEMORY_ARG=""
COASSEMBLY_MEMORY_ARG=""
BIN_REASSEMBLY_THREADS_ARG=""
BIN_REASSEMBLY_MEMORY_ARG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--start-stage)
            START_STAGE="$2"
            shift 2
            ;;
        -e|--end-stage)
            END_STAGE="$2"
            shift 2
            ;;
        -a|--assembly-mode)
            ASSEMBLY_MODE="$2"
            if [ "$ASSEMBLY_MODE" != "individual" ] && [ "$ASSEMBLY_MODE" != "coassembly" ]; then
                echo "Error: Assembly mode must be 'individual' or 'coassembly'"
                exit 1
            fi
            shift 2
            ;;
        -b|--treatment-level-binning)
            TREATMENT_LEVEL_BINNING=true
            shift
            ;;
        --comebin)
            USE_COMEBIN=true
            shift
            ;;
        -t|--treatment)
            SPECIFIC_TREATMENTS+=("$2")
            shift 2
            ;;
        -m|--sample)
            SPECIFIC_SAMPLES+=("$2")
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR_ARG="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR_ARG="$2"
            shift 2
            ;;
        -S|--sample-sheet)
            SAMPLE_SHEET_ARG="$2"
            shift 2
            ;;
        -A|--account)
            SLURM_ACCOUNT_ARG="$2"
            shift 2
            ;;
        --merge-lanes)
            SKIP_MERGE_LANES=false
            shift
            ;;
        --plasmid-detection)
            SKIP_PLASMID_DETECTION=false
            shift
            ;;
        --assembly-threads)
            ASSEMBLY_THREADS_ARG="$2"
            shift 2
            ;;
        --assembly-memory)
            ASSEMBLY_MEMORY_ARG="$2"
            shift 2
            ;;
        --coassembly-memory)
            COASSEMBLY_MEMORY_ARG="$2"
            shift 2
            ;;
        --bin-reassembly-threads)
            BIN_REASSEMBLY_THREADS_ARG="$2"
            shift 2
            ;;
        --bin-reassembly-memory)
            BIN_REASSEMBLY_MEMORY_ARG="$2"
            shift 2
            ;;
        -c|--create-template)
            CREATE_TEMPLATE=true
            shift
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Set environment variables (remove trailing slashes)
if [ -n "$INPUT_DIR_ARG" ]; then
    export INPUT_DIR="${INPUT_DIR_ARG%/}"
fi

if [ -n "$OUTPUT_DIR_ARG" ]; then
    export OUTPUT_DIR="${OUTPUT_DIR_ARG%/}"
fi

if [ -n "$SAMPLE_SHEET_ARG" ]; then
    export SAMPLE_SHEET="$SAMPLE_SHEET_ARG"
fi

if [ -n "$SLURM_ACCOUNT_ARG" ]; then
    export SLURM_ACCOUNT="$SLURM_ACCOUNT_ARG"
fi

# Export assembly parameters if provided
if [ -n "$ASSEMBLY_THREADS_ARG" ]; then
    export ASSEMBLY_THREADS="$ASSEMBLY_THREADS_ARG"
fi

if [ -n "$ASSEMBLY_MEMORY_ARG" ]; then
    export ASSEMBLY_MEMORY="$ASSEMBLY_MEMORY_ARG"
fi

if [ -n "$COASSEMBLY_MEMORY_ARG" ]; then
    export COASSEMBLY_MEMORY="$COASSEMBLY_MEMORY_ARG"
fi

if [ -n "$BIN_REASSEMBLY_THREADS_ARG" ]; then
    export BIN_REASSEMBLY_THREADS="$BIN_REASSEMBLY_THREADS_ARG"
fi

if [ -n "$BIN_REASSEMBLY_MEMORY_ARG" ]; then
    export BIN_REASSEMBLY_MEMORY="$BIN_REASSEMBLY_MEMORY_ARG"
fi

export ASSEMBLY_MODE
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"

# Source configuration after setting environment variables
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Helper function to check if a value is in an array
is_in_array() {
    local value="$1"
    shift
    local array=("$@")
    
    # If array is empty, return true (no filter applied)
    if [ ${#array[@]} -eq 0 ]; then
        return 0
    fi
    
    # Check if value is in array
    for item in "${array[@]}"; do
        if [ "$item" = "$value" ]; then
            return 0
        fi
    done
    
    return 1
}

# Function to validate stage numbers (0-9)
validate_stage() {
    local stage="$1"
    local stage_name="$2"

    # Check if it's a valid stage number (0-9 only)
    if ! [[ "$stage" =~ ^[0-9]$ ]]; then
        echo "Error: $stage_name must be 0-9"
        return 1
    fi
    return 0
}

# Validate inputs
if ! validate_stage "$START_STAGE" "Start stage"; then
    exit 1
fi

if ! validate_stage "$END_STAGE" "End stage"; then
    exit 1
fi

# Compare stages (simple integer comparison)
compare_stages() {
    local stage1="$1"
    local stage2="$2"
    local op="$3"  # gt, ge, lt, le, eq

    case "$op" in
        gt) [ "$stage1" -gt "$stage2" ] ;;
        ge) [ "$stage1" -ge "$stage2" ] ;;
        lt) [ "$stage1" -lt "$stage2" ] ;;
        le) [ "$stage1" -le "$stage2" ] ;;
        eq) [ "$stage1" -eq "$stage2" ] ;;
    esac
}

if compare_stages "$START_STAGE" "$END_STAGE" "gt"; then
    echo "Error: Start stage cannot be greater than end stage"
    exit 1
fi

# Check required directories
if [ -z "$INPUT_DIR" ]; then
    echo "Error: Input directory must be specified with -i or in config"
    exit 1
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Output directory must be specified with -o or in config"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create sample sheet template if requested
if [ "$CREATE_TEMPLATE" = true ]; then
    echo "Creating sample sheet templates..."
    echo ""
    echo "=== Legacy Format Template ==="
    create_sample_sheet_template "${INPUT_DIR}/sample_sheet_legacy_template.tsv"
    echo ""
    echo "=== Multi-Run Format Template ==="
    cat > "${INPUT_DIR}/sample_sheet_multirun_template.csv" << 'EOF'
sample,run,group,short_reads_1,short_reads_2,short_reads_platform
sample1,1,control,sample1_run1_R1.fq.gz,sample1_run1_R2.fq.gz,ILLUMINA
sample1,2,control,sample1_run2_R1.fq.gz,sample1_run2_R2.fq.gz,ILLUMINA
sample2,1,treated,sample2_run1_R1.fq.gz,sample2_run1_R2.fq.gz,ILLUMINA
sample2,2,treated,sample2_run2_R1.fq.gz,sample2_run2_R2.fq.gz,ILLUMINA
EOF
    echo "Multi-run template created: ${INPUT_DIR}/sample_sheet_multirun_template.csv"
    echo ""
    echo "Edit the appropriate template and run the pipeline with -S option."
    echo "The pipeline will auto-detect the format based on column headers."
    exit 0
fi

# Initialize sample information
echo "Initializing sample information..."
init_sample_info

if [ $? -ne 0 ]; then
    echo "Error: Failed to initialize sample information"
    exit 1
fi

# Set file paths for export to SLURM jobs
export TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
export SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"

# Define stage names and scripts (renumbered, removed stage 8 metawrap quant)
declare -A STAGE_NAMES=(
    [0]="Quality Filtering"
    [1]="Validation & Repair"
    [2]="Assembly"
    [3]="Binning"
    [4]="Bin Refinement"
    [5]="Bin Reassembly"
    [6]="MAGpurify"
    [7]="CheckM2"
    [8]="Bin Selection"
    [9]="Bin Collection (per treatment: collect bins, CoverM abundance, GTDB-Tk)"
)

declare -A STAGE_SCRIPTS=(
    [0]="00_quality_filtering.sh"
    [1]="01_validate_repair.sh"
    [2]="02_assembly.sh"  # Will be updated based on assembly mode
    [3]="03_binning.sh"  # Will be updated based on binning mode
    [4]="04_bin_refinement.sh"
    [5]="05_bin_reassembly.sh"
    [6]="06_magpurify.sh"
    [7]="07_checkm2.sh"
    [8]="08_bin_selection.sh"
    [9]="09_bin_collection.sh"
)

# Optional stage scripts (only run if flags are set)
declare -A OPTIONAL_STAGE_NAMES=(
    [merge_lanes]="Lane/Run Detection & Merge"
    [plasmid_detection]="Plasmid Detection"
)

declare -A OPTIONAL_STAGE_SCRIPTS=(
    [merge_lanes]="-01_merge_lanes.sh"
    [plasmid_detection]="plasmid_detection.sh"
)

# Function to get the correct assembly script based on mode
get_assembly_script() {
    if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        echo "02_coassembly.sh"
    else
        echo "02_assembly.sh"
    fi
}

# Function to get the correct binning script based on mode
get_binning_script() {
    if [ "$USE_COMEBIN" = true ]; then
        echo "03b_comebin.sh"
    elif [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        echo "03_binning_treatment_level.sh"
    else
        echo "03_binning.sh"
    fi
}

# Update assembly script based on mode
STAGE_SCRIPTS[2]=$(get_assembly_script)

# Update binning script based on mode
STAGE_SCRIPTS[3]=$(get_binning_script)

# Update stage names to reflect modes
if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
    STAGE_NAMES[2]="Co-Assembly (per treatment)"
else
    STAGE_NAMES[2]="Assembly (per sample)"
fi

if [ "$USE_COMEBIN" = true ]; then
    STAGE_NAMES[3]="Binning (COMEBin)"
    STAGE_NAMES[4]="Bin Refinement (sample-level)"
    STAGE_NAMES[5]="Bin Reassembly (sample-level)"
    STAGE_NAMES[6]="MAGpurify (sample-level)"
    STAGE_NAMES[7]="CheckM2 (sample-level)"
    STAGE_NAMES[8]="Bin Selection (sample-level)"
elif [ "$TREATMENT_LEVEL_BINNING" = true ]; then
    STAGE_NAMES[3]="Binning (treatment-level)"
    STAGE_NAMES[4]="Bin Refinement (treatment-level)"
    STAGE_NAMES[5]="Bin Reassembly (treatment-level)"
    STAGE_NAMES[6]="MAGpurify (treatment-level)"
    STAGE_NAMES[7]="CheckM2 (treatment-level)"
    STAGE_NAMES[8]="Bin Selection (treatment-level)"
else
    STAGE_NAMES[3]="Binning (sample-level)"
    STAGE_NAMES[4]="Bin Refinement (sample-level)"
    STAGE_NAMES[5]="Bin Reassembly (sample-level)"
    STAGE_NAMES[6]="MAGpurify (sample-level)"
    STAGE_NAMES[7]="CheckM2 (sample-level)"
    STAGE_NAMES[8]="Bin Selection (sample-level)"
fi

# Function to check if filters match any samples
validate_filters() {
    local total_samples=$(get_total_samples)
    local matching_samples=()
    local all_treatments=()
    local all_samples=()
    
    # Collect all available samples and treatments
    for i in $(seq 0 $((total_samples - 1))); do
        local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
        if [ -n "$sample_info" ]; then
            IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
            all_treatments+=("$treatment")
            all_samples+=("$sample_name")
            
            # Check if this sample matches filters
            local matches=true
            if ! is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                matches=false
            fi
            if ! is_in_array "$sample_name" "${SPECIFIC_SAMPLES[@]}"; then
                matches=false
            fi
            
            if [ "$matches" = true ]; then
                matching_samples+=("$sample_name ($treatment)")
            fi
        fi
    done
    
    # Remove duplicates from treatments and samples arrays
    local unique_treatments=($(printf '%s\n' "${all_treatments[@]}" | sort -u))
    local unique_samples=($(printf '%s\n' "${all_samples[@]}" | sort -u))
    
    if [ ${#matching_samples[@]} -eq 0 ]; then
        echo ""
        echo "❌ ERROR: No samples match your filters!"
        echo ""
        
        if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
            echo "🏷️  Treatment filter(s): ${SPECIFIC_TREATMENTS[*]}"
            echo "📋 Available treatments:"
            printf '   • %s\n' "${unique_treatments[@]}"
            echo ""
            
            # Look for case-insensitive matches
            local suggestions=()
            for filter in "${SPECIFIC_TREATMENTS[@]}"; do
                for treatment in "${unique_treatments[@]}"; do
                    if [[ "${treatment,,}" == "${filter,,}" ]]; then
                        suggestions+=("$treatment")
                    fi
                done
            done
            
            if [ ${#suggestions[@]} -gt 0 ]; then
                echo "💡 Did you mean one of these? (case-sensitive)"
                printf '   • %s\n' "${suggestions[@]}"
                echo ""
            fi
        fi
        
        if [ ${#SPECIFIC_SAMPLES[@]} -gt 0 ]; then
            echo "🏷️  Sample filter(s): ${SPECIFIC_SAMPLES[*]}"
            echo "📋 Available samples:"
            printf '   • %s\n' "${unique_samples[@]}"
            echo ""
            
            # Look for case-insensitive matches
            local suggestions=()
            for filter in "${SPECIFIC_SAMPLES[@]}"; do
                for sample in "${unique_samples[@]}"; do
                    if [[ "${sample,,}" == "${filter,,}" ]]; then
                        suggestions+=("$sample")
                    fi
                done
            done
            
            if [ ${#suggestions[@]} -gt 0 ]; then
                echo "💡 Did you mean one of these? (case-sensitive)"
                printf '   • %s\n' "${suggestions[@]}"
                echo ""
            fi
        fi
        
        echo "💡 Suggestions:"
        echo "   • Remove filters to process all samples"
        echo "   • Check spelling and case sensitivity"
        echo "   • Use --dry-run to preview what would be processed"
        echo ""
        return 1
    fi
    
    return 0
}

# Function to calculate array size
calculate_array_size() {
    local stage=$1
    
    local total_samples=$(get_total_samples)
    
    # Special handling for stage 2 (Assembly) in coassembly mode
    if [ "$stage" = "2" ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        # Co-assembly: one job per treatment
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
                # Count only filtered treatments
                local count=0
                for treatment in $treatments_list; do
                    if is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        ((count++))
                    fi
                done
                echo "$count"
            else
                echo "$treatments_list" | wc -w
            fi
        else
            echo "0"
        fi
        return
    fi
    
    # Special handling for stage 3 in treatment-level binning mode
    if [ "$stage" = "3" ] && [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        # Treatment-level binning: one job per treatment
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
                # Count only filtered treatments
                local count=0
                for treatment in $treatments_list; do
                    if is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        ((count++))
                    fi
                done
                echo "$count"
            else
                echo "$treatments_list" | wc -w
            fi
        else
            echo "0"
        fi
        return
    fi
    
    # Special handling for stage 4 in treatment-level binning mode
    if [ "$stage" = "4" ] && [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        # Treatment-level bin refinement follows treatment-level binning
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
                # Count only filtered treatments
                local count=0
                for treatment in $treatments_list; do
                    if is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        ((count++))
                    fi
                done
                echo "$count"
            else
                echo "$treatments_list" | wc -w
            fi
        else
            echo "0"
        fi
        return
    fi

    # Special handling for stages 5, 6, 7, 8 in treatment-level binning mode
    if [ "$TREATMENT_LEVEL_BINNING" = true ] && \
       ( [ "$stage" = "5" ] || [ "$stage" = "6" ] || [ "$stage" = "7" ] || [ "$stage" = "8" ] ); then
        # Treatment-level: run once per treatment
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
                # Count only filtered treatments
                local count=0
                for treatment in $treatments_list; do
                    if is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        ((count++))
                    fi
                done
                echo "$count"
            else
                echo "$treatments_list" | wc -w
            fi
        else
            echo "0"
        fi
        return
    fi

    # Sample-level stages (0-8 in individual mode)
    if compare_stages "$stage" "8" "le"; then
        # Filter samples if specific treatment or sample requested
        if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ] || [ ${#SPECIFIC_SAMPLES[@]} -gt 0 ]; then
            local count=0
            for i in $(seq 0 $((total_samples - 1))); do
                local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
                if [ -n "$sample_info" ]; then
                    IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                    
                    # Check filters
                    if ! is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        continue
                    fi
                    if ! is_in_array "$sample_name" "${SPECIFIC_SAMPLES[@]}"; then
                        continue
                    fi
                    
                    ((count++))
                fi
            done
            echo "$count"
        else
            echo "$total_samples"
        fi
    elif [ "$stage" = "9" ]; then
        # Stage 9: Final Report - treatment-level (one job per treatment)
        # This includes bin collection, CoverM abundance, GTDB-Tk, and taxonomy plots
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
                # Count only filtered treatments
                local count=0
                for treatment in $treatments_list; do
                    if is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                        ((count++))
                    fi
                done
                echo "$count"
            else
                echo "$treatments_list" | wc -w
            fi
        else
            echo "0"
        fi
    else
        # Unknown stage
        echo "0"
    fi
}

# Function to get array indices for filtered samples
get_filtered_indices() {
    local stage=$1
    total_samples=$(get_total_samples)
    indices=()
    
    # Special handling for stage 2 (Assembly) in coassembly mode
    if [ "$stage" = "2" ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        # Co-assembly: return treatment indices
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if ! is_in_array "${treatments[$i]}" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                indices+=($i)
            done
        fi
        echo "${indices[@]}"
        return
    fi
    
    # Special handling for stage 3 in treatment-level binning mode
    if [ "$stage" = "3" ] && [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        # Treatment-level binning: return treatment indices
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if ! is_in_array "${treatments[$i]}" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                indices+=($i)
            done
        fi
        echo "${indices[@]}"
        return
    fi
    
    # Special handling for stage 4 in treatment-level binning mode
    if [ "$stage" = "4" ] && [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        # Treatment-level bin refinement: return treatment indices
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if ! is_in_array "${treatments[$i]}" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                indices+=($i)
            done
        fi
        echo "${indices[@]}"
        return
    fi

    # Special handling for stages 5, 6, 7, 8 in treatment-level binning mode
    if [ "$TREATMENT_LEVEL_BINNING" = true ] && \
       ( [ "$stage" = "5" ] || [ "$stage" = "6" ] || [ "$stage" = "7" ] || [ "$stage" = "8" ] ); then
        # Treatment-level: return treatment indices
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if ! is_in_array "${treatments[$i]}" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                indices+=($i)
            done
        fi
        echo "${indices[@]}"
        return
    fi

    # Sample-level stages (0-8)
    if compare_stages "$stage" "8" "le"; then
        for i in $(seq 0 $((total_samples - 1))); do
            sample_info=$(get_sample_info_by_index $i 2>/dev/null)
            if [ -n "$sample_info" ]; then
                IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"

                # Check filters
                if ! is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                if ! is_in_array "$sample_name" "${SPECIFIC_SAMPLES[@]}"; then
                    continue
                fi

                indices+=($i)
            fi
        done
    elif [ "$stage" = "9" ]; then
        # Stage 9: Final Report - treatment-level (includes bin collection, abundance, GTDB-Tk, plots)
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if ! is_in_array "${treatments[$i]}" "${SPECIFIC_TREATMENTS[@]}"; then
                    continue
                fi
                indices+=($i)
            done
        fi
    fi

    echo "${indices[@]}"
}

# Function to submit job
submit_job() {
    local stage=$1
    local script="${STAGE_SCRIPTS[$stage]}"
    local stage_name="${STAGE_NAMES[$stage]}"
    local dependency=$2
    
    if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
        echo "❌ Error: Script not found: ${SCRIPT_DIR}/${script}"
        return 1
    fi
    
    # Create SLURM log directory
    local slurm_log_dir="${OUTPUT_DIR}/logs/slurm"
    mkdir -p "$slurm_log_dir"
    
    # Calculate array size
    array_size=$(calculate_array_size $stage)
    
    if [ $array_size -eq 0 ]; then
        echo "⚠️  No samples/treatments to process for stage $stage"
        echo "   This usually means your filters don't match any samples"
        return 0
    fi
    
    # Build sbatch command
    local cmd="sbatch"
    
    # CRITICAL: Export ALL environment variables to SLURM job
    cmd+=" --export=ALL"
    cmd+=",OUTPUT_DIR=${OUTPUT_DIR}"
    cmd+=",INPUT_DIR=${INPUT_DIR}"
    cmd+=",WORK_DIR=${WORK_DIR}"
    cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR}"
    cmd+=",ASSEMBLY_MODE=${ASSEMBLY_MODE}"
    cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING}"
    cmd+=",USE_COMEBIN=${USE_COMEBIN}"
    cmd+=",TREATMENTS_FILE=${TREATMENTS_FILE}"
    cmd+=",SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE}"
    cmd+=",SLURM_ACCOUNT=${SLURM_ACCOUNT}"

    # Export assembly parameters if set
    if [ -n "$ASSEMBLY_THREADS" ]; then
        cmd+=",ASSEMBLY_THREADS=${ASSEMBLY_THREADS}"
    fi
    if [ -n "$ASSEMBLY_MEMORY" ]; then
        cmd+=",ASSEMBLY_MEMORY=${ASSEMBLY_MEMORY}"
    fi
    if [ -n "$COASSEMBLY_MEMORY" ]; then
        cmd+=",COASSEMBLY_MEMORY=${COASSEMBLY_MEMORY}"
    fi
    if [ -n "$BIN_REASSEMBLY_THREADS" ]; then
        cmd+=",BIN_REASSEMBLY_THREADS=${BIN_REASSEMBLY_THREADS}"
    fi
    if [ -n "$BIN_REASSEMBLY_MEMORY" ]; then
        cmd+=",BIN_REASSEMBLY_MEMORY=${BIN_REASSEMBLY_MEMORY}"
    fi
    
    # Add log file specifications
    local script_basename=$(basename -- "$script" .sh)
    cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
    cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"
    
    # Add dependency if provided
    if [ -n "$dependency" ]; then
        cmd+=" --dependency=afterok:${dependency}"
    fi

    # Add account if provided
    if [ -n "$SLURM_ACCOUNT" ]; then
        cmd+=" --account=${SLURM_ACCOUNT}"
    fi

    # Add array specification
    if [ $array_size -gt 1 ]; then
        indices=($(get_filtered_indices $stage))
        if [ ${#indices[@]} -gt 0 ]; then
            # Create array string
            local array_str=""
            for i in "${indices[@]}"; do
                if [ -z "$array_str" ]; then
                    array_str="$i"
                else
                    array_str="$array_str,$i"
                fi
            done
            cmd+=" --array=${array_str}"
        fi
    elif [ $array_size -eq 1 ]; then
        indices=($(get_filtered_indices $stage))
        if [ ${#indices[@]} -gt 0 ]; then
            cmd+=" --array=${indices[0]}"
        fi
    fi
    
    cmd+=" ${SCRIPT_DIR}/${script}"
    
    if [ "$DRY_RUN" = true ]; then
        echo "🔍 [DRY RUN] Would execute: $cmd"
        echo "   Array size: $array_size"
        echo "   Indices: $(get_filtered_indices $stage)"
        echo "999999"  # Fake job ID for dry run
    else
        echo "🚀 Submitting stage $stage: $stage_name"
        echo "   Array size: $array_size samples/treatments"
        echo "   Command: $cmd"
        
        # Submit job and capture job ID
        result=$($cmd 2>&1)
        if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
            job_id="${BASH_REMATCH[1]}"
            echo "✅ Submitted job $job_id for stage $stage ($stage_name)"
            echo "$job_id"
        else
            echo "❌ Error submitting job: $result"
            echo ""
        fi
    fi
}

# Main execution
echo "========================================"
echo "🧬 Metagenomic Pipeline Execution"
echo "========================================"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Sample sheet: ${SAMPLE_SHEET:-Auto-discovery}"
echo "Assembly mode: $ASSEMBLY_MODE"
if [ "$USE_COMEBIN" = true ]; then
    echo "Binning mode: COMEBin"
else
    echo "Binning mode: $([ "$TREATMENT_LEVEL_BINNING" = true ] && echo "treatment-level" || echo "sample-level")"
fi
echo "SLURM account: ${SLURM_ACCOUNT:-default}"
echo "Start stage: $START_STAGE (${STAGE_NAMES[$START_STAGE]})"
echo "End stage: $END_STAGE (${STAGE_NAMES[$END_STAGE]})"
echo "Total samples: $(get_total_samples)"

if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ]; then
    echo "Specific treatment(s): ${SPECIFIC_TREATMENTS[*]}"
fi

if [ ${#SPECIFIC_SAMPLES[@]} -gt 0 ]; then
    echo "Specific sample(s): ${SPECIFIC_SAMPLES[*]}"
fi

echo "Dry run: $DRY_RUN"
echo "========================================"
echo

# Validate filters before proceeding
if ! validate_filters; then
    exit 1
fi

# Display sample information
echo "📋 Samples to Process:"
echo "======================"
total_samples=$(get_total_samples)
processed_count=0
for i in $(seq 0 $((total_samples - 1))); do
    sample_info=$(get_sample_info_by_index $i 2>/dev/null)
    if [ -n "$sample_info" ]; then
        IFS='|' read -r sample_name treatment r1_path r2_path <<< "$sample_info"
        
        # Check filters
        if ! is_in_array "$treatment" "${SPECIFIC_TREATMENTS[@]}"; then
            continue
        fi
        if ! is_in_array "$sample_name" "${SPECIFIC_SAMPLES[@]}"; then
            continue
        fi
        
        # Count runs/lanes (removed 'local' keyword - not in a function)
        num_runs=$(echo "$r1_path" | tr ',' '\n' | wc -l)
        if [ $num_runs -gt 1 ]; then
            echo "✅ [$i] $sample_name ($treatment) - $num_runs runs/lanes"
        else
            echo "✅ [$i] $sample_name ($treatment)"
        fi
        ((processed_count++))
    fi
done

if [ $processed_count -eq 0 ]; then
    echo "❌ No samples match your filters!"
    exit 1
fi

echo ""
echo "📊 Will process $processed_count sample(s)"
echo ""

# Submit jobs with dependencies
previous_job_id=""
jobs_submitted=0

# Handle optional stages first (if enabled)
OPTIONAL_STAGES_TO_RUN=()

if [ "$SKIP_MERGE_LANES" = false ]; then
    OPTIONAL_STAGES_TO_RUN+=("merge_lanes")
fi

# Submit optional stages
for opt_stage in "${OPTIONAL_STAGES_TO_RUN[@]}"; do
    script="${OPTIONAL_STAGE_SCRIPTS[$opt_stage]}"
    stage_name="${OPTIONAL_STAGE_NAMES[$opt_stage]}"

    echo "📄 Processing optional stage: $stage_name"
    echo "🚀 Submitting optional stage: $stage_name"

    # Build and submit sbatch command similar to submit_job function
    # (simplified version, assumes sample-level processing)
    slurm_log_dir="${OUTPUT_DIR}/logs/slurm"
    mkdir -p "$slurm_log_dir"

    total_samples=$(get_total_samples)
    if [ $total_samples -eq 0 ]; then
        echo "⚠️  No samples to process for optional stage $opt_stage"
        continue
    fi

    cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
    cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
    cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING},USE_COMEBIN=${USE_COMEBIN},TREATMENTS_FILE=${TREATMENTS_FILE}"
    cmd+=",SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE},SLURM_ACCOUNT=${SLURM_ACCOUNT}"

    # Export assembly parameters if set
    if [ -n "$ASSEMBLY_THREADS" ]; then
        cmd+=",ASSEMBLY_THREADS=${ASSEMBLY_THREADS}"
    fi
    if [ -n "$ASSEMBLY_MEMORY" ]; then
        cmd+=",ASSEMBLY_MEMORY=${ASSEMBLY_MEMORY}"
    fi
    if [ -n "$COASSEMBLY_MEMORY" ]; then
        cmd+=",COASSEMBLY_MEMORY=${COASSEMBLY_MEMORY}"
    fi
    if [ -n "$BIN_REASSEMBLY_THREADS" ]; then
        cmd+=",BIN_REASSEMBLY_THREADS=${BIN_REASSEMBLY_THREADS}"
    fi
    if [ -n "$BIN_REASSEMBLY_MEMORY" ]; then
        cmd+=",BIN_REASSEMBLY_MEMORY=${BIN_REASSEMBLY_MEMORY}"
    fi

    if [ -n "$previous_job_id" ]; then
        cmd+=" --dependency=afterok:${previous_job_id}"
    fi

    if [ -n "$SLURM_ACCOUNT" ]; then
        cmd+=" --account=${SLURM_ACCOUNT}"
    fi

    script_basename=$(basename -- "$script" .sh)
    cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
    cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"

    if [ $total_samples -gt 1 ]; then
        cmd+=" --array=0-$((total_samples - 1))"
    elif [ $total_samples -eq 1 ]; then
        cmd+=" --array=0"
    fi

    cmd+=" ${SCRIPT_DIR}/${script}"

    if [ "$DRY_RUN" = true ]; then
        echo "🔍 [DRY RUN] Would execute: $cmd"
        previous_job_id="999999"
    else
        echo "   Command: $cmd"
        result=$($cmd 2>&1)
        if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
            previous_job_id="${BASH_REMATCH[1]}"
            echo "✅ Submitted job $previous_job_id for optional stage ($stage_name)"
            ((jobs_submitted++))
        else
            echo "❌ Error submitting optional stage job: $result"
        fi
    fi
    echo
done

# Now handle main pipeline stages (0-9)
STAGES_TO_RUN=()
for stage in $(seq $START_STAGE $END_STAGE); do
    # Skip plasmid detection stage if it's between other stages
    # Actually, plasmid detection is now optional, so we handle it separately below
    STAGES_TO_RUN+=("$stage")
done

# Submit main pipeline stages
for stage in "${STAGES_TO_RUN[@]}"; do
    echo "📄 Processing stage $stage: ${STAGE_NAMES[$stage]}"

    job_id=$(submit_job "$stage" "$previous_job_id")

    if [ -z "$job_id" ] || [ "$job_id" = "0" ]; then
        echo "❌ Error: Failed to submit job for stage $stage"
        exit 1
    fi

    if [ "$job_id" != "999999" ]; then  # Not a dry run fake ID
        ((jobs_submitted++))
    fi

    previous_job_id="${job_id}"
    echo

    # Insert plasmid detection after assembly if enabled
    if [ "$stage" = "2" ] && [ "$SKIP_PLASMID_DETECTION" = false ]; then
        echo "📄 Processing optional stage: Plasmid Detection"

        script="${OPTIONAL_STAGE_SCRIPTS[plasmid_detection]}"
        stage_name="${OPTIONAL_STAGE_NAMES[plasmid_detection]}"
        slurm_log_dir="${OUTPUT_DIR}/logs/slurm"

        total_samples=$(get_total_samples)
        if [ $total_samples -gt 0 ]; then
            cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
            cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
            cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING},USE_COMEBIN=${USE_COMEBIN},TREATMENTS_FILE=${TREATMENTS_FILE}"
            cmd+=",SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE},SLURM_ACCOUNT=${SLURM_ACCOUNT}"

            # Export assembly parameters if set
            if [ -n "$ASSEMBLY_THREADS" ]; then
                cmd+=",ASSEMBLY_THREADS=${ASSEMBLY_THREADS}"
            fi
            if [ -n "$ASSEMBLY_MEMORY" ]; then
                cmd+=",ASSEMBLY_MEMORY=${ASSEMBLY_MEMORY}"
            fi
            if [ -n "$COASSEMBLY_MEMORY" ]; then
                cmd+=",COASSEMBLY_MEMORY=${COASSEMBLY_MEMORY}"
            fi
            if [ -n "$BIN_REASSEMBLY_THREADS" ]; then
                cmd+=",BIN_REASSEMBLY_THREADS=${BIN_REASSEMBLY_THREADS}"
            fi
            if [ -n "$BIN_REASSEMBLY_MEMORY" ]; then
                cmd+=",BIN_REASSEMBLY_MEMORY=${BIN_REASSEMBLY_MEMORY}"
            fi

            if [ -n "$previous_job_id" ]; then
                cmd+=" --dependency=afterok:${previous_job_id}"
            fi

            if [ -n "$SLURM_ACCOUNT" ]; then
                cmd+=" --account=${SLURM_ACCOUNT}"
            fi

            script_basename=$(basename -- "$script" .sh)
            cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
            cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"

            if [ $total_samples -gt 1 ]; then
                cmd+=" --array=0-$((total_samples - 1))"
            else
                cmd+=" --array=0"
            fi

            cmd+=" ${SCRIPT_DIR}/${script}"

            if [ "$DRY_RUN" = true ]; then
                echo "🔍 [DRY RUN] Would execute: $cmd"
                previous_job_id="999999"
            else
                echo "🚀 Submitting optional stage: $stage_name"
                echo "   Command: $cmd"
                result=$($cmd 2>&1)
                if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
                    previous_job_id="${BASH_REMATCH[1]}"
                    echo "✅ Submitted job $previous_job_id for optional stage ($stage_name)"
                    ((jobs_submitted++))
                else
                    echo "❌ Error submitting optional stage job: $result"
                fi
            fi
            echo
        fi
    fi
done

echo "========================================"
echo "🎉 Pipeline submission complete!"
echo "========================================"

if [ "$DRY_RUN" = false ]; then
    echo "📈 Jobs submitted: $jobs_submitted"
    echo ""
    echo "👀 Monitor job progress with:"
    echo "   squeue -u \$USER"
    echo ""
    echo "📄 View job details with:"
    echo "   sacct -j <job_id>"
    echo ""
    echo "📂 Check logs in:"
    echo "   ${OUTPUT_DIR}/logs/"
    echo ""
    echo "✅ Check progress with:"
    echo "   ls ${OUTPUT_DIR}/checkpoints/*/"
else
    echo "👀 This was a dry run - no jobs were actually submitted"
    echo "   Remove --dry-run to submit jobs for real"
fi