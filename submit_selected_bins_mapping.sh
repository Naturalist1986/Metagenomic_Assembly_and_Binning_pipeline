#!/bin/bash
# submit_selected_bins_mapping.sh - Submit selected bins (MAG) mapping workflow

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Selected Bins (MAG) Mapping Workflow"
echo "====================================================================="
echo ""
echo "This workflow maps reads to the final selected/refined bins (MAGs)"
echo "to measure what fraction of reads map to high-quality MAGs."
echo ""

# Step 1: Check selected bins exist
echo "Step 1: Checking selected bins..."
echo "-------------------------------------------------------------------"

SELECTED_BINS_DIR="${OUTPUT_DIR}/selected_bins"

if [ ! -d "$SELECTED_BINS_DIR" ]; then
    echo "ERROR: Selected bins directory not found: $SELECTED_BINS_DIR"
    echo "Please run bin selection first (step 08_bin_selection.sh)"
    exit 1
fi

# Count treatments
NUM_TREATMENTS=$(ls -1d "${SELECTED_BINS_DIR}"/*/ 2>/dev/null | wc -l)

if [ "$NUM_TREATMENTS" -eq 0 ]; then
    echo "ERROR: No treatments found in selected bins directory"
    echo "Please run bin selection first"
    exit 1
fi

echo "Found $NUM_TREATMENTS treatments with selected bins"
echo ""

# Show bin counts per treatment
echo "Selected bins per treatment:"
for treatment_dir in "${SELECTED_BINS_DIR}"/*/; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    treatment=$(basename "$treatment_dir")
    bin_count=$(ls -1 "${treatment_dir}"/*.fa 2>/dev/null | wc -l)
    echo "  $treatment: $bin_count bins"
done
echo ""

# Step 2: Submit mapping jobs
echo "Step 2: Submitting SLURM array job for MAG mapping..."
echo "-------------------------------------------------------------------"

MAX_ARRAY_INDEX=$((NUM_TREATMENTS - 1))

if [ "$MAX_ARRAY_INDEX" -lt 0 ]; then
    MAX_ARRAY_INDEX=0
fi

echo "Array indices: 0-${MAX_ARRAY_INDEX}"
echo ""

# Check if we're in a SLURM environment
if command -v sbatch &> /dev/null; then
    # Create a temporary submit script with the correct array size
    SUBMIT_SCRIPT=$(mktemp)

    cat > "$SUBMIT_SCRIPT" << EOF
#!/bin/bash
#SBATCH --job-name=mag_map
#SBATCH --array=0-${MAX_ARRAY_INDEX}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel
#SBATCH --output=${LOG_DIR}/selected_bins_mapping/map_%A_%a.out
#SBATCH --error=${LOG_DIR}/selected_bins_mapping/map_%A_%a.err

# Set environment variable for script location
export PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}"

# Run the mapping script
bash "${SCRIPT_DIR}/map_to_selected_bins.sh"
EOF

    # Make log directory
    mkdir -p "${LOG_DIR}/selected_bins_mapping"

    # Submit the job
    JOB_ID=$(sbatch --parsable "$SUBMIT_SCRIPT")
    EXIT_CODE=$?

    rm -f "$SUBMIT_SCRIPT"

    if [ $EXIT_CODE -eq 0 ]; then
        echo "Job submitted successfully!"
        echo "Job ID: $JOB_ID"
        echo ""
        echo "Monitor job status with:"
        echo "  squeue -j $JOB_ID"
        echo ""
        echo "View logs at:"
        echo "  ${LOG_DIR}/selected_bins_mapping/map_${JOB_ID}_*.out"
        echo "  ${LOG_DIR}/selected_bins_mapping/map_${JOB_ID}_*.err"
        echo ""
        echo "After jobs complete, generate combined summary with:"
        echo "  ./summarize_with_mags.sh"
        echo ""
    else
        echo "ERROR: Failed to submit job"
        exit 1
    fi
else
    echo "WARNING: sbatch not found. Not in a SLURM environment."
    echo ""
    echo "To run mapping manually, use:"
    echo "  bash ${SCRIPT_DIR}/map_to_selected_bins.sh"
    echo ""
fi

echo "====================================================================="
echo "Submission complete"
echo "====================================================================="
