#!/bin/bash

# clear_checkpoints.sh - Delete checkpoints to allow re-running pipeline stages
#
# Usage:
#   ./clear_checkpoints.sh -s <stage> [-t <treatment>] [-n <sample_name>] [-d] [-a]
#
# Options:
#   -s <stage>       Stage number to clear (required). Use stage numbers like: -1, 0, 0.5, 1, 2, 7.5, 8, 9, 10
#   -t <treatment>   Treatment name (required for treatment-level stages)
#   -n <sample>      Sample name (required for sample-level stages)
#   -d               Also delete output directories for this stage (default: only delete checkpoints)
#   -a               Clear ALL treatments/samples for this stage (use with caution!)
#   -h               Show this help message
#
# Examples:
#   # Clear bin selection checkpoint for treatment 'hok'
#   ./clear_checkpoints.sh -s 7.5 -t hok
#
#   # Clear bin selection and delete outputs for treatment 'hok'
#   ./clear_checkpoints.sh -s 7.5 -t hok -d
#
#   # Clear quality control checkpoint for sample 'sample1'
#   ./clear_checkpoints.sh -s 1 -n sample1
#
#   # Clear all checkpoints for stage 7.5 (all treatments)
#   ./clear_checkpoints.sh -s 7.5 -a
#
#   # Clear final report checkpoint (single execution stage)
#   ./clear_checkpoints.sh -s 10

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Initialize checkpoint directory if not set
if [ -z "$CHECKPOINT_DIR" ]; then
    CHECKPOINT_DIR="${OUTPUT_DIR}/checkpoints"
fi

# Ensure TREATMENTS_FILE and SAMPLE_INFO_FILE are set
if [ -z "$TREATMENTS_FILE" ]; then
    TREATMENTS_FILE="${OUTPUT_DIR}/treatments.txt"
fi

if [ -z "$SAMPLE_INFO_FILE" ]; then
    SAMPLE_INFO_FILE="${OUTPUT_DIR}/sample_info.csv"
fi

# Initialize variables
STAGE=""
TREATMENT=""
SAMPLE_NAME=""
DELETE_OUTPUT=false
CLEAR_ALL=false

# Display help
show_help() {
    grep '^#' "$0" | grep -v '#!/bin/bash' | sed 's/^# \?//'
    exit 0
}

# Parse command line arguments
while getopts "s:t:n:dah" opt; do
    case $opt in
        s) STAGE="$OPTARG" ;;
        t) TREATMENT="$OPTARG" ;;
        n) SAMPLE_NAME="$OPTARG" ;;
        d) DELETE_OUTPUT=true ;;
        a) CLEAR_ALL=true ;;
        h) show_help ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Validate required arguments
if [ -z "$STAGE" ]; then
    echo "ERROR: Stage number is required (-s)"
    echo "Run with -h for help"
    exit 1
fi

# Map stage numbers to checkpoint names and directories
declare -A STAGE_CHECKPOINT_NAMES=(
    [-1]="preqc"
    [0]="qc"
    [0.5]="validate_repair"
    [1]="qc"
    [2]="assembly"
    [3]="assembly_qc"
    [4]="initial_binning"
    [5]="metawrap_refine"
    [6]="magpurify"
    [7]="checkm2"
    [7.5]="bin_selection"
    [8]="metawrap_quant"
    [9]="bin_collection"
    [10]="final_report"
)

declare -A STAGE_OUTPUT_DIRS=(
    [-1]="preqc"
    [0]="qc"
    [0.5]="validate_repair"
    [1]="qc"
    [2]="assembly"
    [3]="quast"
    [4]="initial_binning"
    [5]="metawrap_refine"
    [6]="magpurify"
    [7]="checkm2"
    [7.5]="selected_bins"
    [8]="metawrap_quant"
    [9]="bin_collection"
    [10]="final_report"
)

# Check if stage is valid
if [ -z "${STAGE_CHECKPOINT_NAMES[$STAGE]}" ]; then
    echo "ERROR: Invalid stage number: $STAGE"
    echo "Valid stages: ${!STAGE_CHECKPOINT_NAMES[@]}"
    exit 1
fi

CHECKPOINT_NAME="${STAGE_CHECKPOINT_NAMES[$STAGE]}"
OUTPUT_DIR_NAME="${STAGE_OUTPUT_DIRS[$STAGE]}"

echo "========================================="
echo "Clear Checkpoints Utility"
echo "========================================="
echo "Stage: $STAGE ($CHECKPOINT_NAME)"
echo "Delete outputs: $DELETE_OUTPUT"
echo "Clear all: $CLEAR_ALL"
echo ""

# Determine processing mode based on stage and user flags
# - If user provides -t flag, assume treatment-level
# - If user provides -n flag, assume sample-level
# - Otherwise use stage-based detection
is_treatment_level_stage() {
    local stage="$1"

    # If user explicitly provided treatment flag, it's treatment-level
    if [ -n "$TREATMENT" ]; then
        return 0
    fi

    # If user explicitly provided sample flag, it's sample-level
    if [ -n "$SAMPLE_NAME" ]; then
        return 1
    fi

    # Stage 10 is single execution (neither treatment nor sample level)
    if [ "$stage" = "10" ]; then
        return 0
    fi

    # Stage 9 is always treatment-level
    if [ "$stage" = "9" ]; then
        return 0
    fi

    # Stages 5-8 (including 7.5) are typically treatment-level in modern pipeline
    # Default to treatment-level for these stages when using -a flag
    case "$stage" in
        5|6|7|7.5|8)
            return 0
            ;;
    esac

    # All other stages are sample-level
    return 1
}

# Function to clear checkpoint for a treatment
clear_treatment_checkpoint() {
    local treatment="$1"
    local checkpoint_file="${CHECKPOINT_DIR}/${treatment}/${CHECKPOINT_NAME}_complete"

    echo "Clearing checkpoint for treatment: $treatment"

    if [ -f "$checkpoint_file" ]; then
        rm -v "$checkpoint_file"
        echo "  ✓ Deleted checkpoint: $checkpoint_file"
    else
        echo "  ℹ Checkpoint not found: $checkpoint_file"
    fi

    # Delete output directory if requested
    if [ "$DELETE_OUTPUT" = true ]; then
        local output_dir="${OUTPUT_DIR}/${OUTPUT_DIR_NAME}/${treatment}"
        if [ -d "$output_dir" ]; then
            echo "  Deleting output directory: $output_dir"
            rm -rfv "$output_dir"
            echo "  ✓ Deleted output directory"
        else
            echo "  ℹ Output directory not found: $output_dir"
        fi
    fi

    echo ""
}

# Function to clear checkpoint for a sample
clear_sample_checkpoint() {
    local sample_name="$1"

    # Find treatment for this sample
    local sample_info=$(grep "^${sample_name}," "$SAMPLE_INFO_FILE" 2>/dev/null | head -1)
    if [ -n "$sample_info" ]; then
        local treatment=$(echo "$sample_info" | cut -d',' -f2)
        local checkpoint_file="${CHECKPOINT_DIR}/${treatment}/${sample_name}_${CHECKPOINT_NAME}_complete"

        echo "Clearing checkpoint for sample: $sample_name (treatment: $treatment)"

        if [ -f "$checkpoint_file" ]; then
            rm -v "$checkpoint_file"
            echo "  ✓ Deleted checkpoint: $checkpoint_file"
        else
            echo "  ℹ Checkpoint not found: $checkpoint_file"
        fi

        # Delete output directory if requested
        if [ "$DELETE_OUTPUT" = true ]; then
            local output_dir="${OUTPUT_DIR}/${OUTPUT_DIR_NAME}/${treatment}/${sample_name}"

            if [ -d "$output_dir" ]; then
                echo "  Deleting output directory: $output_dir"
                rm -rfv "$output_dir"
                echo "  ✓ Deleted output directory"
            else
                echo "  ℹ Output directory not found: $output_dir"
            fi
        fi
    else
        echo "ERROR: Could not find treatment for sample '$sample_name' in SAMPLE_INFO_FILE: $SAMPLE_INFO_FILE"
    fi

    echo ""
}

# Function to clear single execution checkpoint (like final report)
clear_single_checkpoint() {
    local checkpoint_file="${CHECKPOINT_DIR}/${CHECKPOINT_NAME}_complete"

    echo "Clearing checkpoint for stage $STAGE (single execution)"

    if [ -f "$checkpoint_file" ]; then
        rm -v "$checkpoint_file"
        echo "  ✓ Deleted checkpoint: $checkpoint_file"
    else
        echo "  ℹ Checkpoint not found: $checkpoint_file"
    fi

    # Delete output directory if requested
    if [ "$DELETE_OUTPUT" = true ]; then
        local output_dir="${OUTPUT_DIR}/${OUTPUT_DIR_NAME}"
        if [ -d "$output_dir" ]; then
            echo "  Deleting output directory: $output_dir"
            rm -rfv "$output_dir"
            echo "  ✓ Deleted output directory"
        else
            echo "  ℹ Output directory not found: $output_dir"
        fi
    fi

    echo ""
}

# Main logic
if [ "$STAGE" = "10" ]; then
    # Stage 10 is single execution
    clear_single_checkpoint

elif is_treatment_level_stage "$STAGE"; then
    # Treatment-level stage
    if [ "$CLEAR_ALL" = true ]; then
        # Clear all treatments
        echo "Clearing checkpoints for ALL treatments..."
        echo ""

        if [ -f "$TREATMENTS_FILE" ]; then
            while IFS= read -r treatment; do
                treatment=$(echo "$treatment" | tr -d '\r\n' | xargs)
                [ -z "$treatment" ] && continue
                clear_treatment_checkpoint "$treatment"
            done < "$TREATMENTS_FILE"
        else
            echo "ERROR: Treatments file not found: $TREATMENTS_FILE"
            echo ""
            echo "The treatments file is created by run_pipeline.sh when you first run the pipeline."
            echo "If you haven't run the pipeline yet, or if you want to clear a specific treatment,"
            echo "use: ./clear_checkpoints.sh -s $STAGE -t <treatment_name>"
            exit 1
        fi
    else
        # Clear specific treatment
        if [ -z "$TREATMENT" ]; then
            echo "ERROR: Treatment name required for treatment-level stage $STAGE"
            echo ""
            echo "Usage:"
            echo "  ./clear_checkpoints.sh -s $STAGE -t <treatment_name>"
            echo "  ./clear_checkpoints.sh -s $STAGE -a               # Clear all treatments"
            echo ""
            if [ -f "$TREATMENTS_FILE" ]; then
                echo "Available treatments (from $TREATMENTS_FILE):"
                cat "$TREATMENTS_FILE" | head -20 | sed 's/^/  - /'
            fi
            exit 1
        fi

        clear_treatment_checkpoint "$TREATMENT"
    fi

else
    # Sample-level stage
    if [ "$CLEAR_ALL" = true ]; then
        # Clear all samples
        echo "Clearing checkpoints for ALL samples..."
        echo ""

        if [ -f "$SAMPLE_INFO_FILE" ]; then
            while IFS=',' read -r sample_name rest; do
                # Skip header and empty lines
                [ "$sample_name" = "sample_name" ] && continue
                [ -z "$sample_name" ] && continue

                sample_name=$(echo "$sample_name" | tr -d '\r\n' | xargs)
                clear_sample_checkpoint "$sample_name"
            done < "$SAMPLE_INFO_FILE"
        else
            echo "ERROR: Sample info file not found: $SAMPLE_INFO_FILE"
            exit 1
        fi
    else
        # Clear specific sample
        if [ -z "$SAMPLE_NAME" ]; then
            echo "ERROR: Sample name required for sample-level stage $STAGE"
            echo ""
            echo "Usage:"
            echo "  ./clear_checkpoints.sh -s $STAGE -n <sample_name>"
            echo "  ./clear_checkpoints.sh -s $STAGE -a               # Clear all samples"
            echo ""
            if [ -f "$SAMPLE_INFO_FILE" ]; then
                echo "Available samples (from $SAMPLE_INFO_FILE):"
                awk -F',' 'NR>1 {print $1}' "$SAMPLE_INFO_FILE" | head -20 | sed 's/^/  - /'
            fi
            exit 1
        fi

        clear_sample_checkpoint "$SAMPLE_NAME"
    fi
fi

echo "========================================="
echo "Checkpoint clearing complete!"
echo "========================================="
echo ""
echo "You can now re-run stage $STAGE using:"
if [ "$STAGE" = "10" ]; then
    echo "  ./run_pipeline.sh -s $STAGE -e $STAGE"
elif is_treatment_level_stage "$STAGE"; then
    if [ -n "$TREATMENT" ]; then
        echo "  ./run_pipeline.sh -t $TREATMENT -s $STAGE -e $STAGE"
    else
        echo "  ./run_pipeline.sh -s $STAGE -e $STAGE"
    fi
else
    if [ -n "$SAMPLE_NAME" ]; then
        echo "  ./run_pipeline.sh -n $SAMPLE_NAME -s $STAGE -e $STAGE"
    else
        echo "  ./run_pipeline.sh -s $STAGE -e $STAGE"
    fi
fi
echo ""
