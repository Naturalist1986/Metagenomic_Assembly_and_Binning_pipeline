#!/bin/bash

# run_pipeline.sh - Master script to run the complete metagenomic pipeline
# UPDATED: Now supports sample sheets with multiple runs per sample

# Default values
START_STAGE=0  # Start from quality filtering by default (skip optional merge/plasmid stages)
END_STAGE=6    # Updated: removed reassembly, magpurify, and checkm2 stages (Binette runs checkm2)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$SCRIPT_DIR"  # Export for SLURM jobs

SPECIFIC_TREATMENTS=()
SPECIFIC_SAMPLES=()
ASSEMBLY_MODE="individual"
TREATMENT_LEVEL_BINNING=false
# Note: All 5 binners (MetaBAT2, MaxBin2, CONCOCT, COMEBin, SemiBin) now run automatically in stage 3
PER_SAMPLE_CROSS_MAPPING=false  # Map all samples in treatment to each sample's assembly
USE_BINSPREADER=false  # Use BinSPreader for graph-aware bin refinement
USE_BINETTE=false  # Use Binette for consensus binning (replaces DAS Tool)
BINETTE_USE_REFINED=false  # Make Binette use BinSPreader-refined bins instead of original bins
USE_GUNC=false  # Run GUNC for chimerism detection after CheckM2
SKIP_MERGE_LANES=true  # Skip lane merging by default (optional stage)
SKIP_PLASMID_DETECTION=true  # Skip plasmid detection by default (optional stage)
FORCE_RUN=false  # Force re-run stages even if checkpoints exist

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Master script to run the complete metagenomic pipeline.
Supports multiple sample sheet formats including multi-run datasets.

OPTIONS:
    -s, --start-stage NUM         Start from stage NUM (0 to 6) [default: 0]
    -e, --end-stage NUM           End at stage NUM (0 to 6) [default: 6]
    -a, --assembly-mode MODE      Assembly mode: 'individual' or 'coassembly' [default: individual]
    -b, --treatment-level-binning Use treatment-level binning instead of sample-level
    --per-sample-cross-mapping    Map all treatment samples to each sample's assembly (individual mode only)
    --binspreader                 Use BinSPreader for graph-aware bin refinement (requires assembly graph)
    --binette                     Use Binette for consensus binning instead of DAS Tool (requires binette env)
    --binette-use-refined         Make Binette use BinSPreader-refined bins instead of original bins
    --gunc                        Run GUNC chimerism detection on final bins (requires gunc env)
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
    -f, --force                   Force re-run stages by deleting checkpoints (ignores completion status)
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
    3  - Unified Binning (ALL 5 binners: MetaBAT2, MaxBin2, CONCOCT, COMEBin, SemiBin)
    4  - Consensus Binning (Binette - combines all 5 bin sets, runs CheckM2)
    5  - Bin Selection (selects high-quality bins based on Binette's CheckM2 results)
    6  - Bin Collection (per treatment: collects selected bins, runs CoverM & GTDB-Tk)

OPTIONAL STAGES (use flags to enable):
    --merge-lanes        Lane/Run Detection & Merge (auto-detects and merges multiple lanes/runs)
    --plasmid-detection  Plasmid Detection (PlasClass and MOB-suite - runs after assembly)

ADDITIONAL TOOLS (run separately after pipeline completes):
    final_report.sh      Generate cross-treatment comparison plots with taxonomy labels
                         (Requires all treatments to complete stage 6 first)

BINNING MODES:
    By default, binning is performed at the sample level using MetaWRAP (MetaBat2, MaxBin2, Concoct).

    Use --treatment-level-binning to bin at the treatment level instead:
    - All samples in a treatment are binned together
    - Useful for co-assembly workflows or when you want bins across replicates
    - Uses 03_binning_treatment_level.sh instead of 03_binning.sh

    Use --per-sample-cross-mapping with -a individual:
    - Performs individual assembly per sample
    - Maps ALL samples in treatment to EACH sample's assembly
    - Creates bins using multi-sample coverage information
    - Useful for leveraging cross-sample coverage when you want per-sample assemblies

    Use --binspreader for graph-aware bin refinement:
    - Refines bins using assembly graph structure from SPAdes
    - Uses multiple assignment mode (-m) for shared genomic regions
    - Particularly useful for closely related strains sharing core genes
    - Runs after initial binning, before consensus/refinement

    Use --binette for consensus binning:
    - Combines multiple binning results (initial binners + BinSPreader if used)
    - Replaces DAS Tool with more advanced consensus approach
    - Uses CheckM2 for quality-guided bin selection
    - Produces final high-quality bin set

    Use --gunc for chimerism detection:
    - Detects chimeric bins from contamination or strain mixing
    - Runs after CheckM2 quality assessment
    - Flags bins with clade_separation_score > 0.45
    - Essential QC for closely related community members

EXAMPLES:
    # Run complete pipeline (default: stages 0-6, no lane merge or plasmid detection)
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

    # Per-sample assembly with cross-sample mapping (all samples map to each assembly)
    $0 -a individual --per-sample-cross-mapping \\
       -i /path/to/fastq -o /path/to/output

    # Advanced binning with graph-aware refinement
    $0 --binspreader --binette --gunc \\
       -i /path/to/fastq -o /path/to/output

    # Nodule metagenomics: cross-sample mapping + full QC pipeline
    $0 -a individual --per-sample-cross-mapping \\
       --binspreader --binette --gunc \\
       -i /path/to/fastq -o /path/to/output

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

    # Force re-run stages 4-5 (ignoring checkpoints)
    $0 -s 4 -e 5 --force -i /path/to/fastq -o /path/to/output

NOTES:
    - Optional stages (lane merging, plasmid detection) are OFF by default
    - Use --merge-lanes if you have multiple lanes/runs per sample that need merging
    - Use --plasmid-detection to identify plasmid contigs and mobile elements
    - Use --force to re-run stages even if they've already completed (deletes checkpoints)
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
        --semibin)
            USE_SEMIBIN=true
            shift
            ;;
        --per-sample-cross-mapping)
            PER_SAMPLE_CROSS_MAPPING=true
            shift
            ;;
        --binspreader)
            USE_BINSPREADER=true
            shift
            ;;
        --binette)
            USE_BINETTE=true
            shift
            ;;
        --binette-use-refined)
            BINETTE_USE_REFINED=true
            shift
            ;;
        --gunc)
            USE_GUNC=true
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
        -f|--force)
            FORCE_RUN=true
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
export PER_SAMPLE_CROSS_MAPPING
export USE_SEMIBIN
export USE_BINSPREADER
export USE_BINETTE
export BINETTE_USE_REFINED
export USE_GUNC
export FORCE_RUN
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

# Function to validate stage numbers (0-6)
validate_stage() {
    local stage="$1"
    local stage_name="$2"

    # Check if it's a valid stage number (0-6 only)
    if ! [[ "$stage" =~ ^[0-6]$ ]]; then
        echo "Error: $stage_name must be 0-6"
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

# Define stage names and scripts (streamlined: removed reassembly, magpurify, checkm2)
declare -A STAGE_NAMES=(
    [0]="Quality Filtering"
    [1]="Validation & Repair"
    [2]="Assembly"
    [3]="Binning"
    [4]="Bin Refinement"
    [5]="Bin Selection"
    [6]="Bin Collection (per treatment: collect bins, CoverM abundance, GTDB-Tk)"
)

declare -A STAGE_SCRIPTS=(
    [0]="00_quality_filtering.sh"
    [1]="01_validate_repair.sh"
    [2]="02_assembly.sh"  # Will be updated based on assembly mode
    [3]="03_binning.sh"  # Will be updated based on binning mode
    [4]="04_bin_refinement.sh"
    [5]="08_bin_selection.sh"
    [6]="09_bin_collection.sh"
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
    elif [ "$USE_SEMIBIN" = true ]; then
        echo "03c_semibin.sh"
    elif [ "$TREATMENT_LEVEL_BINNING" = true ]; then
        echo "03_binning_treatment_level.sh"
    else
        echo "03_binning.sh"
    fi
}

# Function to get the correct refinement script based on flags
get_refinement_script() {
    # Binette is now the default for stage 4 (uses all 5 bin sets)
    # BinSPreader is optional and handled separately if enabled
    if [ "$USE_BINSPREADER" = true ]; then
        echo "04c_binspreader.sh"
    else
        # Default to Binette for consensus binning of all 5 binners
        echo "04b_binette.sh"
    fi
}

# Update assembly script based on mode
STAGE_SCRIPTS[2]=$(get_assembly_script)

# Update binning script based on mode
STAGE_SCRIPTS[3]=$(get_binning_script)

# Update refinement script based on flags
STAGE_SCRIPTS[4]=$(get_refinement_script)

# Update stage names to reflect modes
if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
    STAGE_NAMES[2]="Co-Assembly (per treatment)"
    STAGE_NAMES[3]="Unified Binning (5 binners - treatment-level)"
    STAGE_NAMES[4]="Consensus Binning (Binette - treatment-level)"
    STAGE_NAMES[5]="Bin Reassembly (treatment-level)"
    STAGE_NAMES[6]="MAGpurify (treatment-level)"
    STAGE_NAMES[7]="CheckM2 (treatment-level)"
    STAGE_NAMES[8]="Bin Selection (treatment-level)"
else
    STAGE_NAMES[2]="Assembly (per sample)"
    STAGE_NAMES[3]="Unified Binning (5 binners - sample-level)"
    STAGE_NAMES[4]="Consensus Binning (Binette - sample-level)"
    STAGE_NAMES[5]="Bin Reassembly (sample-level)"
    STAGE_NAMES[6]="MAGpurify (sample-level)"
    STAGE_NAMES[7]="CheckM2 (sample-level)"
    STAGE_NAMES[8]="Bin Selection (sample-level)"
fi

# Update stage 4 name if BinSPreader is used (optional refinement step)
if [ "$USE_BINSPREADER" = true ]; then
    if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        STAGE_NAMES[4]="Bin Refinement (BinSPreader + Binette - treatment-level)"
    else
        STAGE_NAMES[4]="Bin Refinement (BinSPreader + Binette - sample-level)"
    fi
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
    
    # Special handling for stage 3 in treatment-level binning mode OR COMEBin with coassembly
    if [ "$stage" = "3" ] && { [ "$TREATMENT_LEVEL_BINNING" = true ] || { [ "$USE_COMEBIN" = true ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; }; }; then
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
    
    # Special handling for stage 4 in coassembly or treatment-level binning mode
    if [ "$stage" = "4" ] && ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ); then
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

    # Special handling for stage 5 (bin selection) in coassembly or treatment-level binning mode
    if ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ) && [ "$stage" = "5" ]; then
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

    # Sample-level stages (0-5 in individual mode)
    if compare_stages "$stage" "5" "le"; then
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
    elif [ "$stage" = "6" ]; then
        # Stage 6: Bin Collection - treatment-level (one job per treatment)
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
    
    # Special handling for stage 3 in treatment-level binning mode OR COMEBin with coassembly
    if [ "$stage" = "3" ] && { [ "$TREATMENT_LEVEL_BINNING" = true ] || { [ "$USE_COMEBIN" = true ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; }; }; then
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
    
    # Special handling for stage 4 in coassembly or treatment-level binning mode
    if [ "$stage" = "4" ] && ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ); then
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

    # Special handling for stage 5 (bin selection) in coassembly or treatment-level binning mode
    if ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ) && [ "$stage" = "5" ]; then
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

    # Sample-level stages (0-5)
    if compare_stages "$stage" "5" "le"; then
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
    elif [ "$stage" = "6" ]; then
        # Stage 6: Bin Collection - treatment-level (includes bin collection, abundance, GTDB-Tk)
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
    cmd+=",USE_SEMIBIN=${USE_SEMIBIN}"
    cmd+=",USE_BINSPREADER=${USE_BINSPREADER}"
    cmd+=",USE_BINETTE=${USE_BINETTE}"
    cmd+=",BINETTE_USE_REFINED=${BINETTE_USE_REFINED}"
    cmd+=",USE_GUNC=${USE_GUNC}"
    cmd+=",PER_SAMPLE_CROSS_MAPPING=${PER_SAMPLE_CROSS_MAPPING}"
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
        echo "🔍 [DRY RUN] Would execute: $cmd" >&2
        echo "   Array size: $array_size" >&2
        echo "   Indices: $(get_filtered_indices $stage)" >&2
        echo "999999"  # Fake job ID for dry run
    else
        echo "🚀 Submitting stage $stage: $stage_name" >&2
        echo "   Array size: $array_size samples/treatments" >&2
        echo "   Command: $cmd" >&2

        # Submit job and capture job ID
        result=$($cmd 2>&1)
        if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
            job_id="${BASH_REMATCH[1]}"
            echo "✅ Submitted job $job_id for stage $stage ($stage_name)" >&2
            echo "$job_id"
        else
            echo "❌ Error submitting job: $result" >&2
            echo "" >&2
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
elif [ "$USE_SEMIBIN" = true ]; then
    echo "Binning mode: SemiBin2"
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
    cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING},USE_COMEBIN=${USE_COMEBIN},USE_SEMIBIN=${USE_SEMIBIN}"
    cmd+=",USE_BINSPREADER=${USE_BINSPREADER},USE_BINETTE=${USE_BINETTE},USE_GUNC=${USE_GUNC}"
    cmd+=",PER_SAMPLE_CROSS_MAPPING=${PER_SAMPLE_CROSS_MAPPING},TREATMENTS_FILE=${TREATMENTS_FILE}"
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

# Now handle main pipeline stages (0-6)
STAGES_TO_RUN=()
for stage in $(seq $START_STAGE $END_STAGE); do
    # Skip plasmid detection stage if it's between other stages
    # Actually, plasmid detection is now optional, so we handle it separately below
    STAGES_TO_RUN+=("$stage")
done

# Submit main pipeline stages
for stage in "${STAGES_TO_RUN[@]}"; do
    echo "📄 Processing stage $stage: ${STAGE_NAMES[$stage]}"

    # Delete checkpoints if --force flag is used
    if [ "$FORCE_RUN" = true ]; then
        echo "  ⚠️  Force mode enabled: Deleting existing checkpoints for stage $stage"

        # Determine if this stage runs at treatment-level or sample-level
        if ( [ "$stage" = "2" ] && [ "$ASSEMBLY_MODE" = "coassembly" ] ) || \
           ( [ "$stage" = "3" ] && [ "$TREATMENT_LEVEL_BINNING" = true ] ) || \
           ( [ "$stage" = "4" ] && ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ) ) || \
           ( [ "$stage" = "5" ] && ( [ "$TREATMENT_LEVEL_BINNING" = true ] || [ "$ASSEMBLY_MODE" = "coassembly" ] ) ) || \
           [ "$stage" = "6" ]; then
            # Treatment-level: delete treatment checkpoints
            local treatments_list=$(get_treatments)
            for treatment in $treatments_list; do
                local checkpoint_pattern="${OUTPUT_DIR}/checkpoints/${treatment}/"
                local stage_name="${STAGE_SCRIPTS[$stage]%.sh}"
                stage_name="${stage_name#??_}"  # Remove number prefix (e.g., "03_")

                # Delete various checkpoint patterns for this stage
                rm -f "${checkpoint_pattern}${stage_name}_complete" 2>/dev/null
                rm -f "${checkpoint_pattern}${treatment}_${stage_name}_complete" 2>/dev/null

                echo "    Deleted checkpoints: ${checkpoint_pattern}*${stage_name}*"
            done
        else
            # Sample-level: delete sample checkpoints
            for i in $(seq 0 $(($(get_total_samples) - 1))); do
                local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
                if [ -n "$sample_info" ]; then
                    IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                    local checkpoint_pattern="${OUTPUT_DIR}/checkpoints/${treatment}/"
                    local stage_name="${STAGE_SCRIPTS[$stage]%.sh}"
                    stage_name="${stage_name#??_}"  # Remove number prefix

                    # Delete various checkpoint patterns for this sample
                    rm -f "${checkpoint_pattern}${sample_name}_${stage_name}_complete" 2>/dev/null

                    echo "    Deleted checkpoint: ${checkpoint_pattern}${sample_name}_${stage_name}_complete"
                fi
            done
        fi
    fi

    # Special handling for stage 3: Run unified binning workflow
    if [ "$stage" = "3" ]; then
        # Stage 3a: Create shared BAM files
        echo "  Step 1: Creating shared BAM files (03a_create_shared_bams.sh)"

        bam_script="03a_create_shared_bams.sh"
        slurm_log_dir="${OUTPUT_DIR}/logs/slurm"

        array_size=$(calculate_array_size 3)

        if [ $array_size -gt 0 ]; then
            cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
            cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
            cmd+=",TREATMENTS_FILE=${TREATMENTS_FILE},SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE}"
            cmd+=",SLURM_ACCOUNT=${SLURM_ACCOUNT}"

            if [ -n "$previous_job_id" ]; then
                cmd+=" --dependency=afterok:${previous_job_id}"
            fi

            if [ -n "$SLURM_ACCOUNT" ]; then
                cmd+=" --account=${SLURM_ACCOUNT}"
            fi

            script_basename=$(basename -- "$bam_script" .sh)
            cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
            cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"

            if [ $array_size -gt 1 ]; then
                indices=($(get_filtered_indices 3))
                if [ ${#indices[@]} -gt 0 ]; then
                    array_str=""
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
                indices=($(get_filtered_indices 3))
                if [ ${#indices[@]} -gt 0 ]; then
                    cmd+=" --array=${indices[0]}"
                fi
            fi

            cmd+=" ${SCRIPT_DIR}/${bam_script}"

            if [ "$DRY_RUN" = true ]; then
                echo "  🔍 [DRY RUN] Would execute: $cmd"
                bam_job_id="999999"
            else
                echo "  🚀 Submitting BAM creation job"
                echo "     Command: $cmd"
                result=$($cmd 2>&1)
                if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
                    bam_job_id="${BASH_REMATCH[1]}"
                    echo "  ✅ Submitted BAM creation job $bam_job_id"
                    ((jobs_submitted++))
                else
                    echo "  ❌ Error submitting BAM creation job: $result"
                    exit 1
                fi
            fi
        else
            echo "  ⚠️  No samples/treatments to process for BAM creation"
            bam_job_id=""
        fi

        # Stage 3b: Run all 3 binners in parallel (depend on BAM creation)
        echo "  Step 2: Running all 5 binners (MetaBAT2, MaxBin2, CONCOCT, COMEBin, SemiBin)"

        binner_scripts=("03_binning.sh" "03b_comebin.sh" "03c_semibin.sh")
        binner_names=("MetaBAT2/MaxBin2/CONCOCT" "COMEBin" "SemiBin")
        binner_job_ids=()

        for i in "${!binner_scripts[@]}"; do
            binner_script="${binner_scripts[$i]}"
            binner_name="${binner_names[$i]}"

            echo "    Running $binner_name ($binner_script)..."

            if [ $array_size -gt 0 ]; then
                cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
                cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
                cmd+=",TREATMENTS_FILE=${TREATMENTS_FILE},SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE}"
                cmd+=",SLURM_ACCOUNT=${SLURM_ACCOUNT}"

                # Depend on BAM creation job
                if [ -n "$bam_job_id" ]; then
                    cmd+=" --dependency=afterok:${bam_job_id}"
                fi

                if [ -n "$SLURM_ACCOUNT" ]; then
                    cmd+=" --account=${SLURM_ACCOUNT}"
                fi

                script_basename=$(basename -- "$binner_script" .sh)
                cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
                cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"

                if [ $array_size -gt 1 ]; then
                    indices=($(get_filtered_indices 3))
                    if [ ${#indices[@]} -gt 0 ]; then
                        array_str=""
                        for idx in "${indices[@]}"; do
                            if [ -z "$array_str" ]; then
                                array_str="$idx"
                            else
                                array_str="$array_str,$idx"
                            fi
                        done
                        cmd+=" --array=${array_str}"
                    fi
                elif [ $array_size -eq 1 ]; then
                    indices=($(get_filtered_indices 3))
                    if [ ${#indices[@]} -gt 0 ]; then
                        cmd+=" --array=${indices[0]}"
                    fi
                fi

                cmd+=" ${SCRIPT_DIR}/${binner_script}"

                if [ "$DRY_RUN" = true ]; then
                    echo "    🔍 [DRY RUN] Would execute: $cmd"
                    binner_job_ids+=("999999")
                else
                    echo "    🚀 Submitting $binner_name"
                    result=$($cmd 2>&1)
                    if [[ $result =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
                        job_id="${BASH_REMATCH[1]}"
                        echo "    ✅ Submitted $binner_name job $job_id"
                        binner_job_ids+=("$job_id")
                        ((jobs_submitted++))
                    else
                        echo "    ❌ Error submitting $binner_name job: $result"
                        exit 1
                    fi
                fi
            fi
        done

        # Set previous_job_id to all binner jobs (for stage 4 dependency)
        if [ ${#binner_job_ids[@]} -gt 0 ]; then
            # Join all job IDs with ":"
            previous_job_id="${binner_job_ids[0]}"
            for ((i=1; i<${#binner_job_ids[@]}; i++)); do
                previous_job_id="${previous_job_id}:${binner_job_ids[$i]}"
            done
        fi

        echo ""
        continue
    fi

    # Normal stage submission for non-stage-3 stages
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

    # Insert Binette after BinSPreader if both are enabled
    if [ "$stage" = "4" ] && [ "$USE_BINSPREADER" = true ] && [ "$USE_BINETTE" = true ]; then
        echo "📄 Processing stage 4 (part 2): Binette Consensus Binning"

        binette_script="04b_binette.sh"
        slurm_log_dir="${OUTPUT_DIR}/logs/slurm"

        # Calculate array size for stage 4
        array_size=$(calculate_array_size 4)

        if [ $array_size -gt 0 ]; then
            cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
            cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
            cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING},USE_COMEBIN=${USE_COMEBIN},USE_SEMIBIN=${USE_SEMIBIN}"
            cmd+=",USE_BINSPREADER=${USE_BINSPREADER},USE_BINETTE=${USE_BINETTE},BINETTE_USE_REFINED=${BINETTE_USE_REFINED},USE_GUNC=${USE_GUNC}"
            cmd+=",PER_SAMPLE_CROSS_MAPPING=${PER_SAMPLE_CROSS_MAPPING},TREATMENTS_FILE=${TREATMENTS_FILE}"
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

            # Add dependency on previous BinSPreader job
            if [ -n "$previous_job_id" ]; then
                cmd+=" --dependency=afterok:${previous_job_id}"
            fi

            if [ -n "$SLURM_ACCOUNT" ]; then
                cmd+=" --account=${SLURM_ACCOUNT}"
            fi

            script_basename=$(basename -- "$binette_script" .sh)
            cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
            cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"

            # Build array indices
            if [ $array_size -gt 1 ]; then
                indices=($(get_filtered_indices 4))
                if [ ${#indices[@]} -gt 0 ]; then
                    array_str=""
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
                indices=($(get_filtered_indices 4))
                if [ ${#indices[@]} -gt 0 ]; then
                    cmd+=" --array=${indices[0]}"
                fi
            fi

            cmd+=" ${SCRIPT_DIR}/${binette_script}"

            if [ "$DRY_RUN" = true ]; then
                echo "   [DRY RUN] Would execute: $cmd"
                binette_job_id="999999"
            else
                binette_job_id=$(eval "$cmd" | awk '{print $NF}')
                if [ -n "$binette_job_id" ]; then
                    echo "   ✅ Submitted Binette job: $binette_job_id"
                    ((jobs_submitted++))
                    previous_job_id="$binette_job_id"
                else
                    echo "   ❌ Failed to submit Binette job"
                fi
            fi
        fi
        echo
    fi

    # Insert plasmid detection after assembly if enabled
    if [ "$stage" = "2" ] && [ "$SKIP_PLASMID_DETECTION" = false ]; then
        echo "📄 Processing optional stage: Plasmid Detection"

        script="${OPTIONAL_STAGE_SCRIPTS[plasmid_detection]}"
        stage_name="${OPTIONAL_STAGE_NAMES[plasmid_detection]}"
        slurm_log_dir="${OUTPUT_DIR}/logs/slurm"

        # Calculate filtered sample count (respecting treatment/sample filters)
        total_samples=$(get_total_samples)
        filtered_count=0

        if [ ${#SPECIFIC_TREATMENTS[@]} -gt 0 ] || [ ${#SPECIFIC_SAMPLES[@]} -gt 0 ]; then
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

                    ((filtered_count++))
                fi
            done
        else
            filtered_count=$total_samples
        fi

        if [ $filtered_count -gt 0 ]; then
            cmd="sbatch --export=ALL,OUTPUT_DIR=${OUTPUT_DIR},INPUT_DIR=${INPUT_DIR},WORK_DIR=${WORK_DIR}"
            cmd+=",PIPELINE_SCRIPT_DIR=${PIPELINE_SCRIPT_DIR},ASSEMBLY_MODE=${ASSEMBLY_MODE}"
            cmd+=",TREATMENT_LEVEL_BINNING=${TREATMENT_LEVEL_BINNING},USE_COMEBIN=${USE_COMEBIN},USE_SEMIBIN=${USE_SEMIBIN}"
            cmd+=",USE_BINSPREADER=${USE_BINSPREADER},USE_BINETTE=${USE_BINETTE},USE_GUNC=${USE_GUNC}"
            cmd+=",PER_SAMPLE_CROSS_MAPPING=${PER_SAMPLE_CROSS_MAPPING},TREATMENTS_FILE=${TREATMENTS_FILE}"
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

            # Build filtered array indices
            if [ $filtered_count -gt 1 ]; then
                # Get filtered indices
                filtered_indices=()
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

                        filtered_indices+=($i)
                    fi
                done

                # Create array string
                array_str=""
                for idx in "${filtered_indices[@]}"; do
                    if [ -z "$array_str" ]; then
                        array_str="$idx"
                    else
                        array_str="$array_str,$idx"
                    fi
                done
                cmd+=" --array=${array_str}"
            elif [ $filtered_count -eq 1 ]; then
                # Find the single filtered index
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

                        cmd+=" --array=$i"
                        break
                    fi
                done
            fi

            cmd+=" ${SCRIPT_DIR}/${script}"

            if [ "$DRY_RUN" = true ]; then
                echo "🔍 [DRY RUN] Would execute: $cmd"
                previous_job_id="999999"
            else
                echo "🚀 Submitting optional stage: $stage_name"
                echo "   Array size: $filtered_count samples"
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