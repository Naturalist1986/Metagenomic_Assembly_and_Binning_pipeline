#!/bin/bash

# run_pipeline.sh - Master script to run the complete metagenomic pipeline

# Default values
START_STAGE=-1  # Start from lane merging by default
END_STAGE=11
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$SCRIPT_DIR"  # Export for SLURM jobs

SPECIFIC_TREATMENT=""
SPECIFIC_SAMPLE=""
ASSEMBLY_MODE="individual"

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Master script to run the complete metagenomic pipeline.

OPTIONS:
    -s, --start-stage NUM    Start from stage NUM (-1 to 11) [default: -1]
    -e, --end-stage NUM      End at stage NUM (-1 to 11) [default: 11]
    -a, --assembly-mode MODE Assembly mode: 'individual' or 'coassembly' [default: individual]
    -t, --treatment NAME     Run only for specific treatment
    -m, --sample NAME        Run only for specific sample
    -i, --input-dir PATH     Input directory containing FASTQ files
    -o, --output-dir PATH    Output directory for results
    -S, --sample-sheet PATH  Sample sheet (Excel or TSV)
    -c, --create-template    Create sample sheet template and exit
    -d, --dry-run           Show what would be run without executing
    -h, --help              Show this help message

STAGES:
   -1  - Lane Detection & Merge (auto-detects and merges multiple sequencing lanes)
    0  - Quality Filtering (Trimmomatic)
  0.5  - Validation & Repair (validates paired-end sync, repairs if needed)
    1  - Assembly (MetaSPAdes - individual or co-assembly)
    2  - Plasmid Detection
    3  - Binning
    4  - Bin Refinement
    5  - Bin Reassembly
    6  - MAGpurify
    7  - CheckM2
    8  - CoverM
    9  - Bin Collection
   10  - Combined Report
   11  - Final Report

EXAMPLES:
    # Run complete pipeline (including lane merge and validation)
    $0 -i /path/to/fastq -o /path/to/output

    # Skip lane merging (if already done or single-lane samples)
    $0 --start-stage 0 -i /path/to/fastq -o /path/to/output

    # Run only quality filtering and validation
    $0 --start-stage 0 --end-stage 0.5 -i /path/to/fastq -o /path/to/output

    # Co-assembly mode
    $0 --assembly-mode coassembly -i /path/to/fastq -o /path/to/output

    # Run for specific treatment with co-assembly
    $0 --treatment carR --assembly-mode coassembly -i /path/to/fastq -o /path/to/output

    # Dry run to preview
    $0 --dry-run -i /path/to/fastq -o /path/to/output

EOF
}

# Parse arguments
DRY_RUN=false
CREATE_TEMPLATE=false
INPUT_DIR_ARG=""
OUTPUT_DIR_ARG=""
SAMPLE_SHEET_ARG=""

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
        -t|--treatment)
            SPECIFIC_TREATMENT="$2"
            shift 2
            ;;
        -m|--sample)
            SPECIFIC_SAMPLE="$2"
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

export ASSEMBLY_MODE
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"

# Source configuration after setting environment variables
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Function to validate stage numbers (supports -1 and 0.5)
validate_stage() {
    local stage="$1"
    local stage_name="$2"
    
    # Check if it's a valid stage number
    if [[ "$stage" != "-1" && "$stage" != "0" && "$stage" != "0.5" ]] && \
       ! [[ "$stage" =~ ^[1-9]$|^1[0-1]$ ]]; then
        echo "Error: $stage_name must be -1, 0, 0.5, or 1-11"
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

# Compare stages (handle fractional and negative)
compare_stages() {
    local stage1="$1"
    local stage2="$2"
    local op="$3"  # gt, ge, lt, le, eq
    
    # Convert to comparable format (multiply by 10 and add 100 to handle negatives)
    # -1 becomes 90, 0 becomes 100, 0.5 becomes 105, 1 becomes 110, etc.
    local s1=$(echo "$stage1" | awk '{print ($1 + 10) * 10}')
    local s2=$(echo "$stage2" | awk '{print ($1 + 10) * 10}')
    
    case "$op" in
        gt) [ $(echo "$s1 > $s2" | bc) -eq 1 ] ;;
        ge) [ $(echo "$s1 >= $s2" | bc) -eq 1 ] ;;
        lt) [ $(echo "$s1 < $s2" | bc) -eq 1 ] ;;
        le) [ $(echo "$s1 <= $s2" | bc) -eq 1 ] ;;
        eq) [ $(echo "$s1 == $s2" | bc) -eq 1 ] ;;
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
    echo "Creating sample sheet template..."
    create_sample_sheet_template "${INPUT_DIR}/sample_sheet_template.tsv"
    echo "Template created. Edit the file and run the pipeline with -S option."
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

# Define stage names and scripts
declare -A STAGE_NAMES=(
    [-1]="Lane Detection & Merge"
    [0]="Quality Filtering"
    [0.5]="Validation & Repair"
    [1]="Assembly"
    [2]="Plasmid Detection"
    [3]="Binning"
    [4]="Bin Refinement"
    [5]="Bin Reassembly"
    [6]="MAGpurify"
    [7]="CheckM2"
    [8]="CoverM"
    [9]="Bin Collection"
    [10]="Combined Report"
    [11]="Final Report"
)

declare -A STAGE_SCRIPTS=(
    [-1]="-01_merge_lanes.sh"
    [0]="00_quality_filtering.sh"
    [0.5]="00b_validate_repair.sh"
    [1]="01_assembly.sh"  # Will be updated based on mode
    [2]="02_plasmid_detection.sh"
    [3]="03_binning.sh"
    [4]="04_bin_refinement.sh"
    [5]="05_bin_reassembly.sh"
    [6]="06_magpurify.sh"
    [7]="07_checkm2.sh"
    [8]="08_coverm.sh"
    [9]="09_bin_collection.sh"
    [10]="10_combined_report.sh"
    [11]="11_final_report.sh"
)

# Function to get the correct assembly script based on mode
get_assembly_script() {
    if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        echo "01b_coassembly.sh"
    else
        echo "01_assembly.sh"
    fi
}

# Update assembly script based on mode
STAGE_SCRIPTS[1]=$(get_assembly_script)

# Update stage name to reflect mode
if [ "$ASSEMBLY_MODE" = "coassembly" ]; then
    STAGE_NAMES[1]="Co-Assembly (per treatment)"
else
    STAGE_NAMES[1]="Assembly (per sample)"
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
            if [ -n "$SPECIFIC_TREATMENT" ] && [ "$treatment" != "$SPECIFIC_TREATMENT" ]; then
                matches=false
            fi
            if [ -n "$SPECIFIC_SAMPLE" ] && [ "$sample_name" != "$SPECIFIC_SAMPLE" ]; then
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
        
        if [ -n "$SPECIFIC_TREATMENT" ]; then
            echo "🏷 Treatment filter: '$SPECIFIC_TREATMENT'"
            echo "📋 Available treatments:"
            printf '   • %s\n' "${unique_treatments[@]}"
            echo ""
            
            # Look for case-insensitive matches
            local suggestions=()
            for treatment in "${unique_treatments[@]}"; do
                if [[ "${treatment,,}" == "${SPECIFIC_TREATMENT,,}" ]]; then
                    suggestions+=("$treatment")
                fi
            done
            
            if [ ${#suggestions[@]} -gt 0 ]; then
                echo "💡 Did you mean one of these? (case-sensitive)"
                printf '   • %s\n' "${suggestions[@]}"
                echo ""
            fi
        fi
        
        if [ -n "$SPECIFIC_SAMPLE" ]; then
            echo "🏷 Sample filter: '$SPECIFIC_SAMPLE'"
            echo "📋 Available samples:"
            printf '   • %s\n' "${unique_samples[@]}"
            echo ""
            
            # Look for case-insensitive matches
            local suggestions=()
            for sample in "${unique_samples[@]}"; do
                if [[ "${sample,,}" == "${SPECIFIC_SAMPLE,,}" ]]; then
                    suggestions+=("$sample")
                fi
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
    
    # Special handling for stage 1 in coassembly mode
    if [ "$stage" = "1" ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        # Co-assembly: one job per treatment
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            echo "$treatments_list" | wc -w
        else
            echo "0"
        fi
        return
    fi
    
    # Sample-level stages (-1, 0, 0.5, 1-8 in individual mode)
    if compare_stages "$stage" "8" "le"; then
        # Filter samples if specific treatment or sample requested
        if [ -n "$SPECIFIC_TREATMENT" ] || [ -n "$SPECIFIC_SAMPLE" ]; then
            local count=0
            for i in $(seq 0 $((total_samples - 1))); do
                local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
                if [ -n "$sample_info" ]; then
                    IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                    
                    # Check filters
                    if [ -n "$SPECIFIC_TREATMENT" ] && [ "$treatment" != "$SPECIFIC_TREATMENT" ]; then
                        continue
                    fi
                    if [ -n "$SPECIFIC_SAMPLE" ] && [ "$sample_name" != "$SPECIFIC_SAMPLE" ]; then
                        continue
                    fi
                    
                    ((count++))
                fi
            done
            echo "$count"
        else
            echo "$total_samples"
        fi
    else
        # Treatment-level stages (9-11) - one job per treatment
        local treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            echo "$treatments_list" | wc -w
        else
            echo "0"
        fi
    fi
}

# Function to get array indices for filtered samples
get_filtered_indices() {
    local stage=$1
    total_samples=$(get_total_samples)
    indices=()
    
    # Special handling for stage 1 in coassembly mode
    if [ "$stage" = "1" ] && [ "$ASSEMBLY_MODE" = "coassembly" ]; then
        # Co-assembly: return treatment indices
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if [ -n "$SPECIFIC_TREATMENT" ] && [ "${treatments[$i]}" != "$SPECIFIC_TREATMENT" ]; then
                    continue
                fi
                indices+=($i)
            done
        fi
        echo "${indices[@]}"
        return
    fi
    
    # Sample-level stages
    if compare_stages "$stage" "8" "le"; then
        for i in $(seq 0 $((total_samples - 1))); do
            sample_info=$(get_sample_info_by_index $i 2>/dev/null)
            if [ -n "$sample_info" ]; then
                IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                
                # Check filters
                if [ -n "$SPECIFIC_TREATMENT" ] && [ "$treatment" != "$SPECIFIC_TREATMENT" ]; then
                    continue
                fi
                if [ -n "$SPECIFIC_SAMPLE" ] && [ "$sample_name" != "$SPECIFIC_SAMPLE" ]; then
                    continue
                fi
                
                indices+=($i)
            fi
        done
    else
        # Treatment-level stages
        treatments_list=$(get_treatments)
        if [ -n "$treatments_list" ]; then
            local treatments=($treatments_list)
            for i in "${!treatments[@]}"; do
                if [ -n "$SPECIFIC_TREATMENT" ] && [ "${treatments[$i]}" != "$SPECIFIC_TREATMENT" ]; then
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
    cmd+=",TREATMENTS_FILE=${TREATMENTS_FILE}"
    cmd+=",SAMPLE_INFO_FILE=${SAMPLE_INFO_FILE}"
    
    # Add log file specifications
    local script_basename=$(basename -- "$script" .sh)
    cmd+=" --output=${slurm_log_dir}/${script_basename}_%A_%a.log"
    cmd+=" --error=${slurm_log_dir}/${script_basename}_%A_%a.err"
    
    # Add dependency if provided
    if [ -n "$dependency" ]; then
        cmd+=" --dependency=afterok:${dependency}"
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
echo "Start stage: $START_STAGE (${STAGE_NAMES[$START_STAGE]})"
echo "End stage: $END_STAGE (${STAGE_NAMES[$END_STAGE]})"
echo "Total samples: $(get_total_samples)"

if [ -n "$SPECIFIC_TREATMENT" ]; then
    echo "Specific treatment: $SPECIFIC_TREATMENT"
fi

if [ -n "$SPECIFIC_SAMPLE" ]; then
    echo "Specific sample: $SPECIFIC_SAMPLE"
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
        if [ -n "$SPECIFIC_TREATMENT" ] && [ "$treatment" != "$SPECIFIC_TREATMENT" ]; then
            continue
        fi
        if [ -n "$SPECIFIC_SAMPLE" ] && [ "$sample_name" != "$SPECIFIC_SAMPLE" ]; then
            continue
        fi
        
        echo "✅ [$i] $sample_name ($treatment)"
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

# Create list of stages to run (handle fractional stages)
STAGES_TO_RUN=()
current_stage=$START_STAGE

while compare_stages "$current_stage" "$END_STAGE" "le"; do
    STAGES_TO_RUN+=("$current_stage")
    
    # Increment stage (handle 0.5 case)
    if [ "$current_stage" = "0" ] && compare_stages "0.5" "$END_STAGE" "le"; then
        current_stage="0.5"
    elif [ "$current_stage" = "0.5" ]; then
        current_stage="1"
    elif [ "$current_stage" = "-1" ]; then
        current_stage="0"
    else
        current_stage=$((current_stage + 1))
    fi
done

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