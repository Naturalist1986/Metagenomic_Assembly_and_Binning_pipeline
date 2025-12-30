#!/bin/bash

# submit_assembly_success_rates.sh - Helper script to submit array job with correct size
#
# This script counts the number of treatments and submits the assembly success
# rate calculation with the correct array size.

# Check required environment variables
if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR not set"
    echo "Please set OUTPUT_DIR to your pipeline output directory:"
    echo "  export OUTPUT_DIR=/path/to/your/output"
    exit 1
fi

if [ -z "$PIPELINE_DIR" ]; then
    echo "ERROR: PIPELINE_DIR not set"
    echo "Please set PIPELINE_DIR to your pipeline scripts directory:"
    echo "  export PIPELINE_DIR=/path/to/Metagenomic_Assembly_and_Binning_pipeline"
    exit 1
fi

# Default values
MODE="both"
MAX_PARALLEL=10
EXTRA_ARGS=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--mode)
            MODE="$2"
            EXTRA_ARGS="$EXTRA_ARGS --mode $2"
            shift 2
            ;;
        -p|--max-parallel)
            MAX_PARALLEL="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -m, --mode MODE           Assembly mode: individual, coassembly, or both (default: both)"
            echo "  -p, --max-parallel NUM    Maximum parallel jobs (default: 10)"
            echo "  -h, --help                Show this help message"
            echo ""
            echo "Required Environment Variables:"
            echo "  OUTPUT_DIR                Path to pipeline output directory"
            echo "  PIPELINE_DIR              Path to pipeline scripts directory"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Count treatments
echo "Counting treatments in $OUTPUT_DIR..."

declare -A TREATMENTS
TREATMENT_COUNT=0

# Count from assembly directory
if [ "$MODE" = "individual" ] || [ "$MODE" = "both" ]; then
    if [ -d "${OUTPUT_DIR}/assembly" ]; then
        for treatment_dir in "${OUTPUT_DIR}/assembly"/*; do
            if [ -d "$treatment_dir" ]; then
                treatment=$(basename "$treatment_dir")
                TREATMENTS["$treatment"]=1
            fi
        done
    fi
fi

# Count from coassembly directory
if [ "$MODE" = "coassembly" ] || [ "$MODE" = "both" ]; then
    if [ -d "${OUTPUT_DIR}/coassembly" ]; then
        for treatment_dir in "${OUTPUT_DIR}/coassembly"/*; do
            if [ -d "$treatment_dir" ]; then
                treatment=$(basename "$treatment_dir")
                TREATMENTS["$treatment"]=1
            fi
        done
    fi
fi

# Get unique treatment count
TREATMENT_COUNT=${#TREATMENTS[@]}

if [ $TREATMENT_COUNT -eq 0 ]; then
    echo "ERROR: No treatments found!"
    echo "  Checked assembly: ${OUTPUT_DIR}/assembly"
    echo "  Checked coassembly: ${OUTPUT_DIR}/coassembly"
    exit 1
fi

echo "Found $TREATMENT_COUNT unique treatment(s):"
for treatment in "${!TREATMENTS[@]}"; do
    echo "  - $treatment"
done
echo ""

# Calculate array range
MAX_INDEX=$((TREATMENT_COUNT - 1))
ARRAY_SPEC="0-${MAX_INDEX}%${MAX_PARALLEL}"

echo "Submitting SLURM array job..."
echo "  Array specification: $ARRAY_SPEC"
echo "  Mode: $MODE"
echo "  Output directory: $OUTPUT_DIR"
echo "  Pipeline directory: $PIPELINE_DIR"
echo ""

# Submit job
JOB_ID=$(sbatch \
    --array="$ARRAY_SPEC" \
    --export=OUTPUT_DIR,PIPELINE_DIR \
    "${PIPELINE_DIR}/calculate_assembly_success_rates.sh" \
    $EXTRA_ARGS | awk '{print $NF}')

if [ $? -eq 0 ]; then
    echo "✓ Job submitted successfully!"
    echo "  Job ID: $JOB_ID"
    echo ""
    echo "Monitor job status:"
    echo "  squeue -j $JOB_ID"
    echo ""
    echo "View array task details:"
    echo "  squeue -j $JOB_ID -r"
    echo ""
    echo "Output files will be in:"
    echo "  Log files: slurm-${JOB_ID}_*.out (one per treatment)"
    echo "  Summary: assembly_success_rates_summary.tsv"
    echo "  Per-assembly: \$OUTPUT_DIR/assembly/{treatment}/{sample}/assembly_success_rate.txt"
    echo "  Per-coassembly: \$OUTPUT_DIR/coassembly/{treatment}/assembly_success_rate.txt"
else
    echo "✗ Job submission failed!"
    exit 1
fi
