#!/bin/bash

# Helper script to submit bin annotation jobs

set -euo pipefail

# Default values
TREATMENT=""
BINNING_BASE="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/binning"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRY_RUN=false

usage() {
    cat << EOF
Usage: $0 -t TREATMENT [OPTIONS]

Submit Prodigal + MicrobeAnnotator annotation jobs for bins.

REQUIRED:
    -t, --treatment NAME    Treatment/group name (e.g., carR)

OPTIONS:
    -b, --binning-base DIR  Base binning directory [default: $BINNING_BASE]
    -d, --dry-run           Show what would be submitted without actually submitting
    -h, --help              Show this help message

EXAMPLE:
    $0 -t carR

This will:
1. Find all COMEBin and SemiBin bin directories for the treatment
2. Submit an array job that runs Prodigal + MicrobeAnnotator on each bin directory
3. Output annotations to: binning/TREATMENT/SAMPLE/annotations/BINNER/

EOF
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -t|--treatment)
            TREATMENT="$2"
            shift 2
            ;;
        -b|--binning-base)
            BINNING_BASE="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [ -z "$TREATMENT" ]; then
    echo "ERROR: Treatment name is required (-t)"
    usage
fi

# Check if treatment directory exists
if [ ! -d "${BINNING_BASE}/${TREATMENT}" ]; then
    echo "ERROR: Treatment directory not found: ${BINNING_BASE}/${TREATMENT}"
    exit 1
fi

# Count bin directories to process
echo "Scanning for bin directories in treatment: $TREATMENT"
echo ""

BIN_DIRS=()

for sample_dir in "${BINNING_BASE}/${TREATMENT}"/*; do
    if [ ! -d "$sample_dir" ]; then
        continue
    fi

    sample_name=$(basename "$sample_dir")

    # Check for COMEBin bins
    comebin_dir="${sample_dir}/comebin/comebin_res/comebin_res_bins"
    if [ ! -d "$comebin_dir" ]; then
        comebin_dir="${sample_dir}/comebin/comebin_res_bins"
    fi

    if [ -d "$comebin_dir" ] && [ "$(ls -A ${comebin_dir}/*.fa 2>/dev/null)" ]; then
        bin_count=$(ls -1 ${comebin_dir}/*.fa 2>/dev/null | wc -l)
        echo "  ✓ Found COMEBin bins: ${sample_name} ($bin_count bins)"
        BIN_DIRS+=("${sample_name}|comebin|${comebin_dir}")
    fi

    # Check for SemiBin bins
    semibin_dir="${sample_dir}/semibin/output_bins"
    if [ -d "$semibin_dir" ] && [ "$(ls -A ${semibin_dir}/*.fa 2>/dev/null)" ]; then
        bin_count=$(ls -1 ${semibin_dir}/*.fa 2>/dev/null | wc -l)
        echo "  ✓ Found SemiBin bins: ${sample_name} ($bin_count bins)"
        BIN_DIRS+=("${sample_name}|semibin|${semibin_dir}")
    fi
done

echo ""
echo "Total bin directories to annotate: ${#BIN_DIRS[@]}"

if [ ${#BIN_DIRS[@]} -eq 0 ]; then
    echo "ERROR: No bin directories found for treatment $TREATMENT"
    exit 1
fi

# Create logs directory
mkdir -p "${BINNING_BASE}/../logs/slurm"

# Build sbatch command
ARRAY_SIZE=${#BIN_DIRS[@]}
ARRAY_MAX=$((ARRAY_SIZE - 1))

CMD="sbatch"
CMD+=" --export=ALL,TREATMENT=${TREATMENT}"
CMD+=" --array=0-${ARRAY_MAX}"
CMD+=" ${SCRIPT_DIR}/annotate_bins.sh"

echo ""
echo "Submitting annotation job array..."
echo "Command: $CMD"
echo ""

if [ "$DRY_RUN" = true ]; then
    echo "[DRY RUN] Would execute: $CMD"
    exit 0
fi

# Submit the job
JOB_ID=$(eval "$CMD" | awk '{print $NF}')

if [ -n "$JOB_ID" ]; then
    echo "✅ Job submitted successfully: $JOB_ID"
    echo ""
    echo "Monitor with: squeue -j $JOB_ID"
    echo "Cancel with: scancel $JOB_ID"
    echo ""
    echo "Output will be in:"
    echo "  Proteins: binning/${TREATMENT}/*/annotations/*/proteins/"
    echo "  Annotations: binning/${TREATMENT}/*/annotations/*/microbeannotator/"
else
    echo "❌ Failed to submit job"
    exit 1
fi
