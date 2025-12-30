#!/bin/bash
#SBATCH --job-name=assembly_success
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# calculate_assembly_success_rates.sh - Standalone script to calculate assembly success rates
# for already-completed assemblies (both individual and coassembly)
#
# Can run in two modes:
# 1. Array job mode (recommended for many treatments): Each array task processes one treatment
# 2. Sequential mode: Process all treatments in one job

# Usage:
#   ./calculate_assembly_success_rates.sh [OPTIONS]
#
# Options:
#   -m, --mode MODE          Assembly mode: individual, coassembly, or both (default: both)
#   -s, --sample NAME        Process specific sample (for individual mode)
#   -t, --treatment NAME     Process specific treatment
#   -o, --output FILE        Output summary report file (default: assembly_success_rates_summary.tsv)
#   -h, --help               Show this help message
#
# Required Environment Variables:
#   OUTPUT_DIR               Path to pipeline output directory
#   PIPELINE_DIR             Path to pipeline scripts directory (where 00_config_utilities.sh is located)
#                            If not set, will try to detect automatically

# Default values
MODE="both"
SAMPLE_NAME=""
TREATMENT_NAME=""
OUTPUT_SUMMARY="assembly_success_rates_summary.tsv"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"
            shift 2
            ;;
        -s|--sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -t|--treatment)
            TREATMENT_NAME="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_SUMMARY="$2"
            shift 2
            ;;
        -h|--help)
            grep '^#' "$0" | grep -v '#!/bin/bash' | sed 's/^# //'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Determine pipeline directory
if [ -n "$PIPELINE_DIR" ]; then
    # User explicitly set PIPELINE_DIR
    CONFIG_FILE="${PIPELINE_DIR}/00_config_utilities.sh"
elif [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    # Legacy variable
    CONFIG_FILE="${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    # Try to auto-detect (works when running script directly, not via SLURM)
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
    if [ -f "${SCRIPT_DIR}/00_config_utilities.sh" ]; then
        CONFIG_FILE="${SCRIPT_DIR}/00_config_utilities.sh"
    else
        echo "ERROR: Cannot find 00_config_utilities.sh"
        echo ""
        echo "Please set PIPELINE_DIR environment variable to the directory containing"
        echo "the pipeline scripts (where 00_config_utilities.sh is located)."
        echo ""
        echo "Example:"
        echo "  export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline"
        echo "  sbatch calculate_assembly_success_rates.sh"
        echo ""
        echo "Or set it inline:"
        echo "  PIPELINE_DIR=/path/to/scripts OUTPUT_DIR=/path/to/output sbatch calculate_assembly_success_rates.sh"
        exit 1
    fi
fi

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file not found: $CONFIG_FILE"
    echo ""
    echo "Please ensure PIPELINE_DIR is set correctly:"
    echo "  export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline"
    exit 1
fi

# Source configuration and utilities
source "$CONFIG_FILE"

# Initialize
init_conda
TEMP_DIR=$(setup_temp_dir)

log "====== Starting Assembly Success Rate Calculation ======"
log "Mode: $MODE"
log "Output summary: $OUTPUT_SUMMARY"

# Detect if running as array job
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
    ARRAY_MODE=true
    ARRAY_INDEX=$SLURM_ARRAY_TASK_ID
    log "Running in ARRAY mode (task $ARRAY_INDEX)"
else
    ARRAY_MODE=false
    log "Running in SEQUENTIAL mode"
fi

# Create summary file header (append mode for array jobs)
if [ "$ARRAY_MODE" = false ] || [ "$ARRAY_INDEX" -eq 0 ]; then
    echo -e "Type\tTreatment\tSample\tTotal_Reads\tAssembly_Success_Rate(%)\tContigs_File\tStatus" > "$OUTPUT_SUMMARY"
fi

# Function to get all unique treatments
get_all_treatments() {
    local treatments=()

    # Get treatments from assembly directory
    if [ "$MODE" = "individual" ] || [ "$MODE" = "both" ]; then
        if [ -d "${OUTPUT_DIR}/assembly" ]; then
            for treatment_dir in "${OUTPUT_DIR}/assembly"/*; do
                if [ -d "$treatment_dir" ]; then
                    local treatment=$(basename "$treatment_dir")
                    treatments+=("$treatment")
                fi
            done
        fi
    fi

    # Get treatments from coassembly directory
    if [ "$MODE" = "coassembly" ] || [ "$MODE" = "both" ]; then
        if [ -d "${OUTPUT_DIR}/coassembly" ]; then
            for treatment_dir in "${OUTPUT_DIR}/coassembly"/*; do
                if [ -d "$treatment_dir" ]; then
                    local treatment=$(basename "$treatment_dir")
                    # Add only if not already in array
                    if [[ ! " ${treatments[@]} " =~ " ${treatment} " ]]; then
                        treatments+=("$treatment")
                    fi
                fi
            done
        fi
    fi

    # Return unique treatments
    printf '%s\n' "${treatments[@]}" | sort -u
}

# Function to calculate success rate for individual assembly
calculate_individual_assembly_rate() {
    local sample_name="$1"
    local treatment="$2"

    log "Processing individual assembly: $sample_name ($treatment)"

    local assembly_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    local contigs_file="${assembly_dir}/contigs.fasta"

    # Check if assembly exists
    if [ ! -f "$contigs_file" ] || [ ! -s "$contigs_file" ]; then
        log "  WARNING: Assembly not found or empty: $contigs_file"
        echo -e "individual\t${treatment}\t${sample_name}\t0\t0.00\t${contigs_file}\tno_assembly" >> "$OUTPUT_SUMMARY"
        return 1
    fi

    # Check if already calculated
    if [ -f "${assembly_dir}/assembly_success_rate.txt" ]; then
        local existing_rate=$(cat "${assembly_dir}/assembly_success_rate.txt")
        log "  Assembly success rate already calculated: ${existing_rate}%"

        # Get read count for summary
        local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
        local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"

        local input_r1=""
        if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
            input_r1="${validated_dir}/validated_1.fastq.gz"
        elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
            input_r1="${quality_dir}/filtered_1.fastq.gz"
        fi

        local total_reads=0
        if [ -n "$input_r1" ] && [ -f "$input_r1" ]; then
            total_reads=$(count_reads "$input_r1")
        fi

        echo -e "individual\t${treatment}\t${sample_name}\t${total_reads}\t${existing_rate}\t${contigs_file}\talready_calculated" >> "$OUTPUT_SUMMARY"
        return 0
    fi

    # Find input reads (prefer validated, fall back to quality_filtering)
    local validated_dir="${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"

    local input_r1=""
    local input_r2=""

    if [ -d "$validated_dir" ] && [ -f "${validated_dir}/validated_1.fastq.gz" ]; then
        input_r1="${validated_dir}/validated_1.fastq.gz"
        input_r2="${validated_dir}/validated_2.fastq.gz"
        log "  Using validated reads"
    elif [ -d "$quality_dir" ] && [ -f "${quality_dir}/filtered_1.fastq.gz" ]; then
        input_r1="${quality_dir}/filtered_1.fastq.gz"
        input_r2="${quality_dir}/filtered_2.fastq.gz"
        log "  Using quality-filtered reads"
    else
        log "  ERROR: No input reads found for $sample_name"
        echo -e "individual\t${treatment}\t${sample_name}\t0\t0.00\t${contigs_file}\tno_reads" >> "$OUTPUT_SUMMARY"
        return 1
    fi

    # Validate input files exist
    if [ ! -f "$input_r1" ] || [ ! -f "$input_r2" ]; then
        log "  ERROR: Input read files not found"
        echo -e "individual\t${treatment}\t${sample_name}\t0\t0.00\t${contigs_file}\tno_reads" >> "$OUTPUT_SUMMARY"
        return 1
    fi

    # Count total reads
    local total_reads=$(count_reads "$input_r1")
    log "  Total reads: $total_reads"

    # Calculate success rate
    local assembly_success_rate=$(calculate_assembly_success_rate \
        "$contigs_file" \
        "$input_r1" \
        "$input_r2" \
        "$sample_name" \
        "${assembly_dir}/assembly_mapping")

    if [ $? -eq 0 ]; then
        log "  Assembly success rate: ${assembly_success_rate}%"

        # Save to file
        echo "$assembly_success_rate" > "${assembly_dir}/assembly_success_rate.txt"

        # Update assembly statistics file if it exists
        local stats_file="${assembly_dir}/assembly_statistics.txt"
        if [ -f "$stats_file" ]; then
            # Check if success rate already in file
            if ! grep -q "Assembly success rate:" "$stats_file"; then
                echo "" >> "$stats_file"
                echo "Assembly Quality:" >> "$stats_file"
                echo "  Assembly success rate: ${assembly_success_rate}%" >> "$stats_file"
                echo "  (Percentage of input reads that mapped to assembled contigs)" >> "$stats_file"
                log "  Updated assembly statistics file"
            fi
        fi

        echo -e "individual\t${treatment}\t${sample_name}\t${total_reads}\t${assembly_success_rate}\t${contigs_file}\tsuccess" >> "$OUTPUT_SUMMARY"
        return 0
    else
        log "  ERROR: Failed to calculate assembly success rate"
        echo -e "individual\t${treatment}\t${sample_name}\t${total_reads}\t0.00\t${contigs_file}\tfailed" >> "$OUTPUT_SUMMARY"
        return 1
    fi
}

# Function to calculate success rate for coassembly
calculate_coassembly_rate() {
    local treatment="$1"

    log "Processing coassembly: $treatment"

    local coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"
    local contigs_file="${coassembly_dir}/contigs.fasta"

    # Check if coassembly exists
    if [ ! -f "$contigs_file" ] || [ ! -s "$contigs_file" ]; then
        log "  WARNING: Coassembly not found or empty: $contigs_file"
        echo -e "coassembly\t${treatment}\tALL\t0\t0.00\t${contigs_file}\tno_assembly" >> "$OUTPUT_SUMMARY"
        return 1
    fi

    # Check if already calculated
    if [ -f "${coassembly_dir}/assembly_success_rate.txt" ]; then
        local existing_rate=$(cat "${coassembly_dir}/assembly_success_rate.txt")
        log "  Assembly success rate already calculated: ${existing_rate}%"

        # Get read count for summary
        local merged_reads_dir="${coassembly_dir}/merged_reads"
        local total_reads=0
        if [ -f "${merged_reads_dir}/merged_R1.fastq.gz" ]; then
            total_reads=$(count_reads "${merged_reads_dir}/merged_R1.fastq.gz")
        fi

        echo -e "coassembly\t${treatment}\tALL\t${total_reads}\t${existing_rate}\t${contigs_file}\talready_calculated" >> "$OUTPUT_SUMMARY"
        return 0
    fi

    # Find merged reads
    local merged_reads_dir="${coassembly_dir}/merged_reads"
    local merged_r1="${merged_reads_dir}/merged_R1.fastq.gz"
    local merged_r2="${merged_reads_dir}/merged_R2.fastq.gz"

    if [ ! -f "$merged_r1" ] || [ ! -f "$merged_r2" ]; then
        log "  ERROR: Merged reads not found for coassembly"
        log "    Expected: $merged_r1 and $merged_r2"
        echo -e "coassembly\t${treatment}\tALL\t0\t0.00\t${contigs_file}\tno_reads" >> "$OUTPUT_SUMMARY"
        return 1
    fi

    log "  Using merged reads from coassembly"

    # Count total reads
    local total_reads=$(count_reads "$merged_r1")
    log "  Total reads: $total_reads"

    # Calculate success rate
    local assembly_success_rate=$(calculate_coassembly_success_rate \
        "$contigs_file" \
        "$merged_r1" \
        "$merged_r2" \
        "$treatment" \
        "${coassembly_dir}/coassembly_mapping")

    if [ $? -eq 0 ]; then
        log "  Assembly success rate: ${assembly_success_rate}%"

        # Save to file
        echo "$assembly_success_rate" > "${coassembly_dir}/assembly_success_rate.txt"

        echo -e "coassembly\t${treatment}\tALL\t${total_reads}\t${assembly_success_rate}\t${contigs_file}\tsuccess" >> "$OUTPUT_SUMMARY"
        return 0
    else
        log "  ERROR: Failed to calculate assembly success rate"
        echo -e "coassembly\t${treatment}\tALL\t${total_reads}\t0.00\t${contigs_file}\tfailed" >> "$OUTPUT_SUMMARY"
        return 1
    fi
}

# Determine which treatments to process
if [ "$ARRAY_MODE" = true ]; then
    # Array mode: Get treatment for this array index
    mapfile -t ALL_TREATMENTS < <(get_all_treatments)

    if [ "$ARRAY_INDEX" -ge "${#ALL_TREATMENTS[@]}" ]; then
        log "No treatment for array index $ARRAY_INDEX (only ${#ALL_TREATMENTS[@]} treatments found)"
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    CURRENT_TREATMENT="${ALL_TREATMENTS[$ARRAY_INDEX]}"
    log "Processing treatment: $CURRENT_TREATMENT (array task $ARRAY_INDEX)"
    TREATMENTS_TO_PROCESS=("$CURRENT_TREATMENT")
else
    # Sequential mode: Process all treatments (or specific treatment if requested)
    if [ -n "$TREATMENT_NAME" ]; then
        TREATMENTS_TO_PROCESS=("$TREATMENT_NAME")
        log "Processing specific treatment: $TREATMENT_NAME"
    else
        mapfile -t TREATMENTS_TO_PROCESS < <(get_all_treatments)
        log "Processing all treatments: ${TREATMENTS_TO_PROCESS[*]}"
    fi
fi

# Process individual assemblies
if [ "$MODE" = "individual" ] || [ "$MODE" = "both" ]; then
    log "====== Processing Individual Assemblies ======"

    if [ -n "$SAMPLE_NAME" ] && [ -n "$TREATMENT_NAME" ]; then
        # Process specific sample
        calculate_individual_assembly_rate "$SAMPLE_NAME" "$TREATMENT_NAME"
    else
        # Process samples for selected treatments
        if [ -d "${OUTPUT_DIR}/assembly" ]; then
            for treatment in "${TREATMENTS_TO_PROCESS[@]}"; do
                treatment_dir="${OUTPUT_DIR}/assembly/${treatment}"

                if [ ! -d "$treatment_dir" ]; then
                    log "  No individual assemblies found for treatment: $treatment"
                    continue
                fi

                log "  Processing individual assemblies for treatment: $treatment"

                for sample_dir in "$treatment_dir"/*; do
                    if [ ! -d "$sample_dir" ]; then
                        continue
                    fi

                    sample=$(basename "$sample_dir")
                    calculate_individual_assembly_rate "$sample" "$treatment"
                done
            done
        else
            log "WARNING: No individual assemblies found at ${OUTPUT_DIR}/assembly"
        fi
    fi
fi

# Process coassemblies
if [ "$MODE" = "coassembly" ] || [ "$MODE" = "both" ]; then
    log "====== Processing Coassemblies ======"

    if [ -n "$SAMPLE_NAME" ]; then
        log "WARNING: --sample option ignored for coassembly mode"
    fi

    # Process coassemblies for selected treatments
    if [ -d "${OUTPUT_DIR}/coassembly" ]; then
        for treatment in "${TREATMENTS_TO_PROCESS[@]}"; do
            coassembly_dir="${OUTPUT_DIR}/coassembly/${treatment}"

            if [ ! -d "$coassembly_dir" ]; then
                log "  No coassembly found for treatment: $treatment"
                continue
            fi

            log "  Processing coassembly for treatment: $treatment"
            calculate_coassembly_rate "$treatment"
        done
    else
        log "WARNING: No coassemblies found at ${OUTPUT_DIR}/coassembly"
    fi
fi

# Generate summary statistics
log "====== Summary ======"
log "Results saved to: $OUTPUT_SUMMARY"

# Count successes and failures
total_processed=$(tail -n +2 "$OUTPUT_SUMMARY" | wc -l)
successful=$(tail -n +2 "$OUTPUT_SUMMARY" | grep -c "success")
already_calculated=$(tail -n +2 "$OUTPUT_SUMMARY" | grep -c "already_calculated")
failed=$(tail -n +2 "$OUTPUT_SUMMARY" | grep -c "failed")
no_assembly=$(tail -n +2 "$OUTPUT_SUMMARY" | grep -c "no_assembly")
no_reads=$(tail -n +2 "$OUTPUT_SUMMARY" | grep -c "no_reads")

log "Total processed: $total_processed"
log "  Newly calculated: $successful"
log "  Already calculated: $already_calculated"
log "  Failed: $failed"
log "  No assembly found: $no_assembly"
log "  No reads found: $no_reads"

# Calculate average success rate for newly calculated
if [ $successful -gt 0 ]; then
    avg_rate=$(tail -n +2 "$OUTPUT_SUMMARY" | grep "success" | awk -F'\t' '{sum+=$5; count++} END {if(count>0) printf "%.2f", sum/count; else print "0.00"}')
    log "Average assembly success rate (newly calculated): ${avg_rate}%"
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"

log "====== Assembly Success Rate Calculation Complete ======"
