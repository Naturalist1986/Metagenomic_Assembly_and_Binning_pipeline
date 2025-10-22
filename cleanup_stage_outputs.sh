#!/bin/bash

# cleanup_stage_outputs.sh - Clean stage output directories for restarting pipeline stages
# This script removes output directories and checkpoints for specified stages and samples

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Try to source configuration utilities
CONFIG_SOURCED=false
for config_file in "${SCRIPT_DIR}/00_config_utilities.sh" "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"; do
    if [ -f "$config_file" ]; then
        echo "Sourcing configuration from: $config_file"
        source "$config_file" 2>/dev/null && CONFIG_SOURCED=true && break
    fi
done

if [ "$CONFIG_SOURCED" = false ]; then
    echo "ERROR: Could not source config utilities."
    echo "Please ensure 00_config_utilities.sh is available in one of these locations:"
    echo "  - ${SCRIPT_DIR}/00_config_utilities.sh"
    echo "  - ${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
    echo ""
    echo "Alternatively, set these environment variables manually:"
    echo "  - OUTPUT_DIR"
    echo "  - CHECKPOINT_DIR" 
    echo "  - SAMPLE_INFO_FILE"
    exit 1
fi

# Validate and set up required variables
if [ -z "$OUTPUT_DIR" ] || [ "$OUTPUT_DIR" = "/path/to/output/directory" ]; then
    echo "ERROR: OUTPUT_DIR is not properly configured."
    echo "Current value: ${OUTPUT_DIR:-[NOT SET]}"
    echo ""
    echo "Please set OUTPUT_DIR to your actual pipeline output directory:"
    echo "  export OUTPUT_DIR=/sci/backup/aerez/aerez/moshea/Pipeline_Test"
    echo "  $0 [your options]"
    echo ""
    echo "Or update your 00_config_utilities.sh file to set the correct path."
    exit 1
fi

# Set CHECKPOINT_DIR if not explicitly set (default to OUTPUT_DIR/checkpoints)
if [ -z "$CHECKPOINT_DIR" ]; then
    export CHECKPOINT_DIR="${OUTPUT_DIR}/checkpoints"
    echo "CHECKPOINT_DIR not set, using: $CHECKPOINT_DIR"
fi

# Available stages
AVAILABLE_STAGES=(
    "read_qc"
    "assembly"
    "binning" 
    "bin_refinement"
    "bin_reassembly"
    "bin_classification"
    "bin_annotation"
    "bin_quantification"
)

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Clean pipeline stage outputs for restarting specific stages.

OPTIONS:
    -s, --stage STAGE           Stage(s) to clean (can be specified multiple times)
    -a, --all-stages           Clean all pipeline stages
    -S, --sample SAMPLE        Specific sample(s) to clean (can be specified multiple times)  
    -T, --treatment TREATMENT  Specific treatment(s) to clean (can be specified multiple times)
    -A, --all-samples          Clean all samples
    -c, --checkpoints-only     Only remove checkpoints, keep output files
    -f, --force                Force removal without confirmation
    -n, --dry-run              Show what would be removed without actually removing
    -v, --verbose              Verbose output
    --show-config              Show current configuration and exit
    -h, --help                 Display this help message

STAGES:
$(printf '    %s\n' "${AVAILABLE_STAGES[@]}")

EXAMPLES:
    # Show current configuration
    $0 --show-config

    # Clean bin_refinement stage for all samples
    $0 -s bin_refinement -A

    # Clean multiple stages for specific samples
    $0 -s binning -s bin_refinement -S sample1 -S sample2

    # Clean all stages for specific treatment
    $0 -a -T treatment_A

    # Dry run to see what would be cleaned
    $0 -s assembly -A --dry-run --verbose

    # Only remove checkpoints without deleting output files
    $0 -s bin_refinement -A --checkpoints-only

    # Force cleanup without confirmation
    $0 -s binning -A --force

CONFIGURATION:
    This script requires proper configuration variables. If you see errors about
    missing variables, ensure your 00_config_utilities.sh file defines:
    
    export OUTPUT_DIR="/path/to/pipeline/output"
    export CHECKPOINT_DIR="/path/to/pipeline/checkpoints"
    
    Or set them manually before running this script:
    export OUTPUT_DIR="/your/output/path"
    export CHECKPOINT_DIR="/your/checkpoint/path"
    ./cleanup_stage_outputs.sh -s bin_refinement -A

EOF
}

# Initialize variables
STAGES=()
SAMPLES=()
TREATMENTS=()
ALL_STAGES=false
ALL_SAMPLES=false
CHECKPOINTS_ONLY=false
FORCE=false
DRY_RUN=false
VERBOSE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--stage)
            if [[ -n "$2" ]] && [[ "${AVAILABLE_STAGES[*]}" =~ $2 ]]; then
                STAGES+=("$2")
                shift 2
            else
                echo "ERROR: Invalid or missing stage: $2"
                echo "Available stages: ${AVAILABLE_STAGES[*]}"
                exit 1
            fi
            ;;
        -a|--all-stages)
            ALL_STAGES=true
            shift
            ;;
        -S|--sample)
            if [[ -n "$2" ]]; then
                SAMPLES+=("$2")
                shift 2
            else
                echo "ERROR: Missing sample name"
                exit 1
            fi
            ;;
        -T|--treatment)
            if [[ -n "$2" ]]; then
                TREATMENTS+=("$2")
                shift 2
            else
                echo "ERROR: Missing treatment name"
                exit 1
            fi
            ;;
        -A|--all-samples)
            ALL_SAMPLES=true
            shift
            ;;
        -c|--checkpoints-only)
            CHECKPOINTS_ONLY=true
            shift
            ;;
        -f|--force)
            FORCE=true
            shift
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        --show-config)
            echo "Current Configuration:"
            echo "  Script directory: $SCRIPT_DIR"
            echo "  CONFIG_SOURCED: $CONFIG_SOURCED"
            echo "  OUTPUT_DIR: ${OUTPUT_DIR:-[NOT SET]}"
            echo "  CHECKPOINT_DIR: ${CHECKPOINT_DIR:-[NOT SET]}"
            echo "  PIPELINE_SCRIPT_DIR: ${PIPELINE_SCRIPT_DIR:-[NOT SET]}"
            echo ""
            echo "Available functions:"
            echo "  get_sample_count: $(command -v get_sample_count >/dev/null && echo "Available" || echo "Not available")"
            echo "  get_sample_info_by_index: $(command -v get_sample_info_by_index >/dev/null && echo "Available" || echo "Not available")"
            echo ""
            echo "Output directory contents:"
            if [ -d "$OUTPUT_DIR" ]; then
                ls -la "$OUTPUT_DIR"
            else
                echo "  Output directory does not exist: $OUTPUT_DIR"
            fi
            exit 0
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validation
if [ ${#STAGES[@]} -eq 0 ] && [ "$ALL_STAGES" = false ]; then
    echo "ERROR: No stages specified. Use -s/--stage or -a/--all-stages"
    usage
    exit 1
fi

if [ ${#SAMPLES[@]} -eq 0 ] && [ ${#TREATMENTS[@]} -eq 0 ] && [ "$ALL_SAMPLES" = false ]; then
    echo "ERROR: No samples or treatments specified. Use -S/--sample, -T/--treatment, or -A/--all-samples"
    usage
    exit 1
fi

# Set stages to clean
if [ "$ALL_STAGES" = true ]; then
    STAGES=("${AVAILABLE_STAGES[@]}")
fi

# Verbose logging function
vlog() {
    if [ "$VERBOSE" = true ]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    fi
}

# Get all samples to process
get_samples_to_process() {
    local samples_to_process=()
    
    # Try to use config utility functions first
    if command -v get_sample_count &> /dev/null && command -v get_sample_info_by_index &> /dev/null; then
        # Use config utility functions
        vlog "Using config utility functions to get sample information"
        
        if [ "$ALL_SAMPLES" = true ]; then
            local sample_count=$(get_sample_count)
            for ((i=0; i<sample_count; i++)); do
                local sample_info=$(get_sample_info_by_index $i)
                if [ -n "$sample_info" ]; then
                    IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                    samples_to_process+=("${sample_name}|${treatment}")
                fi
            done
        else
            # Process specific samples and treatments using utility functions
            local sample_count=$(get_sample_count)
            
            if [ ${#SAMPLES[@]} -gt 0 ]; then
                for sample in "${SAMPLES[@]}"; do
                    local found=false
                    for ((i=0; i<sample_count; i++)); do
                        local sample_info=$(get_sample_info_by_index $i)
                        if [ -n "$sample_info" ]; then
                            IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                            if [ "$sample_name" = "$sample" ]; then
                                samples_to_process+=("${sample_name}|${treatment}")
                                found=true
                            fi
                        fi
                    done
                    if [ "$found" = false ]; then
                        echo "WARNING: Sample '$sample' not found in sample info"
                    fi
                done
            fi
            
            if [ ${#TREATMENTS[@]} -gt 0 ]; then
                for treatment in "${TREATMENTS[@]}"; do
                    local found=false
                    for ((i=0; i<sample_count; i++)); do
                        local sample_info=$(get_sample_info_by_index $i)
                        if [ -n "$sample_info" ]; then
                            IFS='|' read -r sample_name treat _ _ <<< "$sample_info"
                            if [ "$treat" = "$treatment" ]; then
                                samples_to_process+=("${sample_name}|${treat}")
                                found=true
                            fi
                        fi
                    done
                    if [ "$found" = false ]; then
                        echo "WARNING: Treatment '$treatment' not found in sample info"
                    fi
                done
            fi
        fi
    else
        # Fallback: scan output directory structure
        vlog "Config utility functions not available, scanning output directory structure"
        
        if [ "$ALL_SAMPLES" = true ]; then
            # Find all sample directories by scanning the output structure
            for stage_dir in "$OUTPUT_DIR"/*; do
                if [ -d "$stage_dir" ]; then
                    for treatment_dir in "$stage_dir"/*; do
                        if [ -d "$treatment_dir" ]; then
                            local treatment=$(basename "$treatment_dir")
                            for sample_dir in "$treatment_dir"/*; do
                                if [ -d "$sample_dir" ]; then
                                    local sample_name=$(basename "$sample_dir")
                                    samples_to_process+=("${sample_name}|${treatment}")
                                fi
                            done
                        fi
                    done
                fi
            done
        else
            # Process specific samples/treatments by scanning directories
            if [ ${#SAMPLES[@]} -gt 0 ]; then
                for sample in "${SAMPLES[@]}"; do
                    local found=false
                    for stage_dir in "$OUTPUT_DIR"/*; do
                        if [ -d "$stage_dir" ]; then
                            for treatment_dir in "$stage_dir"/*; do
                                if [ -d "$treatment_dir" ]; then
                                    local treatment=$(basename "$treatment_dir")
                                    if [ -d "$treatment_dir/$sample" ]; then
                                        samples_to_process+=("${sample}|${treatment}")
                                        found=true
                                    fi
                                fi
                            done
                        fi
                    done
                    if [ "$found" = false ]; then
                        echo "WARNING: Sample '$sample' not found in output directories"
                    fi
                done
            fi
            
            if [ ${#TREATMENTS[@]} -gt 0 ]; then
                for treatment in "${TREATMENTS[@]}"; do
                    local found=false
                    for stage_dir in "$OUTPUT_DIR"/*; do
                        if [ -d "$stage_dir/$treatment" ]; then
                            found=true
                            for sample_dir in "$stage_dir/$treatment"/*; do
                                if [ -d "$sample_dir" ]; then
                                    local sample_name=$(basename "$sample_dir")
                                    samples_to_process+=("${sample_name}|${treatment}")
                                fi
                            done
                        fi
                    done
                    if [ "$found" = false ]; then
                        echo "WARNING: Treatment '$treatment' not found in output directories"
                    fi
                done
            fi
        fi
    fi
    
    # Remove duplicates
    printf '%s\n' "${samples_to_process[@]}" | sort -u
}

# Remove checkpoint for a sample and stage
remove_checkpoint() {
    local sample_name="$1"
    local treatment="$2"
    local stage="$3"
    
    # Match the checkpoint pattern used in the config utilities
    local checkpoint_dir="${CHECKPOINT_DIR}/${treatment}"
    local checkpoint_file="${checkpoint_dir}/${sample_name}_${stage}_complete"
    
    if [ -f "$checkpoint_file" ]; then
        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY RUN] Would remove checkpoint: $checkpoint_file"
        else
            vlog "Removing checkpoint: $checkpoint_file"
            rm -f "$checkpoint_file"
            echo "  Removed checkpoint: $checkpoint_file"
        fi
        return 0
    else
        vlog "Checkpoint not found: $checkpoint_file"
        return 1
    fi
}

# Remove output directory for a stage
remove_stage_output() {
    local sample_name="$1"
    local treatment="$2"
    local stage="$3"
    
    local output_dir="${OUTPUT_DIR}/${stage}/${treatment}/${sample_name}"
    
    if [ -d "$output_dir" ]; then
        if [ "$DRY_RUN" = true ]; then
            echo "  [DRY RUN] Would remove directory: $output_dir"
            if [ "$VERBOSE" = true ]; then
                echo "    Directory size: $(du -sh "$output_dir" 2>/dev/null | cut -f1)"
                echo "    Files: $(find "$output_dir" -type f | wc -l)"
            fi
        else
            vlog "Removing directory: $output_dir"
            rm -rf "$output_dir"
            echo "  Removed directory: $output_dir"
        fi
        return 0
    else
        vlog "Output directory not found: $output_dir"
        return 1
    fi
}

# Clean a specific stage for a sample
clean_sample_stage() {
    local sample_name="$1"
    local treatment="$2"
    local stage="$3"
    
    echo "Cleaning stage '$stage' for sample '$sample_name' (treatment: $treatment)"
    
    local cleaned_something=false
    
    # Always remove checkpoint
    if remove_checkpoint "$sample_name" "$stage"; then
        cleaned_something=true
    fi
    
    # Remove output directory unless checkpoints-only mode
    if [ "$CHECKPOINTS_ONLY" = false ]; then
        if remove_stage_output "$sample_name" "$treatment" "$stage"; then
            cleaned_something=true
        fi
    fi
    
    if [ "$cleaned_something" = false ]; then
        echo "  No files found to clean for this sample/stage combination"
    fi
}

# Main cleaning function
perform_cleanup() {
    local samples_processed=0
    local stages_processed=0
    
    echo "Starting cleanup process..."
    echo "Stages to clean: ${STAGES[*]}"
    
    # Get samples to process
    local samples_to_process
    readarray -t samples_to_process < <(get_samples_to_process)
    
    if [ ${#samples_to_process[@]} -eq 0 ]; then
        echo "ERROR: No samples found to process"
        echo ""
        echo "This could mean:"
        echo "  1. No matching samples/treatments found"
        echo "  2. The output directory structure is different than expected"
        echo "  3. Configuration issues with sample information"
        echo ""
        echo "Debug information:"
        echo "  Output directory: $OUTPUT_DIR"
        echo "  Available stage directories:"
        if [ -d "$OUTPUT_DIR" ]; then
            ls -la "$OUTPUT_DIR" | head -10
        else
            echo "    Output directory does not exist!"
        fi
        exit 1
    fi
    
    echo "Samples to process: ${#samples_to_process[@]}"
    
    if [ "$VERBOSE" = true ]; then
        echo "Sample list:"
        printf '  %s\n' "${samples_to_process[@]}"
        echo ""
    fi
    
    # Process each sample
    for sample_treatment in "${samples_to_process[@]}"; do
        if [ -z "$sample_treatment" ]; then
            echo "WARNING: Empty sample_treatment entry, skipping..."
            continue
        fi
        
        IFS='|' read -r sample_name treatment <<< "$sample_treatment"
        
        if [ -z "$sample_name" ]; then
            echo "WARNING: Empty sample name in '$sample_treatment', skipping..."
            continue
        fi
        
        echo ""
        echo "Processing sample: $sample_name (treatment: $treatment)"
        
        # Process each stage for this sample
        for stage in "${STAGES[@]}"; do
            clean_sample_stage "$sample_name" "$treatment" "$stage"
            ((stages_processed++))
        done
        
        ((samples_processed++))
    done
    
    echo ""
    echo "Cleanup completed:"
    echo "  Samples processed: $samples_processed"
    echo "  Stage-sample combinations processed: $stages_processed"
    
    if [ "$DRY_RUN" = true ]; then
        echo "  (This was a dry run - no files were actually removed)"
    fi
}

# Confirmation function
confirm_cleanup() {
    if [ "$FORCE" = true ] || [ "$DRY_RUN" = true ]; then
        return 0
    fi
    
    echo ""
    echo "WARNING: This will permanently delete the following:"
    echo "  Stages: ${STAGES[*]}"
    
    if [ "$ALL_SAMPLES" = true ]; then
        echo "  Scope: ALL samples"
    else
        if [ ${#SAMPLES[@]} -gt 0 ]; then
            echo "  Samples: ${SAMPLES[*]}"
        fi
        if [ ${#TREATMENTS[@]} -gt 0 ]; then
            echo "  Treatments: ${TREATMENTS[*]}"
        fi
    fi
    
    if [ "$CHECKPOINTS_ONLY" = true ]; then
        echo "  Action: Remove checkpoints only (keep output files)"
    else
        echo "  Action: Remove checkpoints AND output directories"
    fi
    
    echo ""
    read -p "Are you sure you want to continue? [y/N]: " -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cleanup cancelled"
        exit 0
    fi
}

# Main execution
main() {
    # Show configuration for debugging
    if [ "$VERBOSE" = true ]; then
        echo "Configuration:"
        echo "  OUTPUT_DIR: ${OUTPUT_DIR:-[NOT SET]}"
        echo "  CHECKPOINT_DIR: ${CHECKPOINT_DIR:-[NOT SET]}"
        echo "  PIPELINE_SCRIPT_DIR: ${PIPELINE_SCRIPT_DIR:-[NOT SET]}"
        echo ""
    fi
    
    # Check if required directories exist
    if [ ! -d "$OUTPUT_DIR" ]; then
        echo "ERROR: Output directory not found: $OUTPUT_DIR"
        echo "Please check your configuration."
        exit 1
    fi
    
    if [ ! -d "$CHECKPOINT_DIR" ]; then
        echo "WARNING: Checkpoint directory not found: $CHECKPOINT_DIR"
        echo "Creating checkpoint directory..."
        mkdir -p "$CHECKPOINT_DIR" || {
            echo "ERROR: Could not create checkpoint directory: $CHECKPOINT_DIR"
            exit 1
        }
    fi
    
    # Show what will be cleaned and ask for confirmation
    confirm_cleanup
    
    # Perform the cleanup
    perform_cleanup
}

# Run main function
main "$@"