#!/bin/bash
# submit_eukfinder.sh - Identify largest bins and submit EukFinder jobs

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Parameters
NUM_LARGEST="${1:-2}"  # Number of largest bins per binner (default: 2)

echo "====================================================================="
echo "EukFinder Submission Workflow"
echo "====================================================================="
echo ""

# Step 1: Identify largest bins
echo "Step 1: Identifying largest bins..."
echo "-------------------------------------------------------------------"

LARGEST_BINS_FILE="${OUTPUT_DIR}/largest_bins_list.txt"

"${SCRIPT_DIR}/identify_largest_bins.sh" "$NUM_LARGEST" "$LARGEST_BINS_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to identify largest bins"
    exit 1
fi

# Check if any bins were found
TOTAL_BINS=$(grep -v "^#" "$LARGEST_BINS_FILE" | wc -l)

if [ "$TOTAL_BINS" -eq 0 ]; then
    echo "ERROR: No bins found to process"
    echo "Please check that binning has completed successfully"
    exit 1
fi

echo ""
echo "Found $TOTAL_BINS bins to process with EukFinder"
echo ""

# Check if taxonomy database has been updated
TAXA_DB="${HOME}/.etetoolkit/taxa.sqlite"
if [ ! -f "$TAXA_DB" ]; then
    echo "WARNING: NCBI taxonomy database not found at $TAXA_DB"
    echo ""
    echo "To avoid database locking issues with parallel jobs, you should"
    echo "update the taxonomy database once before running:"
    echo ""
    echo "  ./update_eukfinder_taxonomy.sh"
    echo ""
    read -p "Do you want to continue anyway? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted. Please run ./update_eukfinder_taxonomy.sh first"
        exit 1
    fi
elif [ "$TOTAL_BINS" -gt 1 ]; then
    echo "NOTE: Running $TOTAL_BINS parallel jobs"
    echo "If you encounter 'database is locked' errors, run this first:"
    echo "  ./update_eukfinder_taxonomy.sh"
    echo ""
fi

# Step 2: Calculate array size
MAX_ARRAY_INDEX=$((TOTAL_BINS - 1))

if [ "$MAX_ARRAY_INDEX" -lt 0 ]; then
    MAX_ARRAY_INDEX=0
fi

echo "Step 2: Submitting SLURM array job..."
echo "-------------------------------------------------------------------"
echo "Array indices: 0-${MAX_ARRAY_INDEX}"
echo ""

# Step 3: Submit the job
# Check if we're in a SLURM environment
if command -v sbatch &> /dev/null; then
    # Create a temporary submit script with the correct array size
    SUBMIT_SCRIPT=$(mktemp)

    cat > "$SUBMIT_SCRIPT" << EOF
#!/bin/bash
#SBATCH --job-name=eukfinder
#SBATCH --array=0-${MAX_ARRAY_INDEX}
#SBATCH --account=ofinkel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=${LOG_DIR}/eukfinder_output/eukfinder_%A_%a.out
#SBATCH --error=${LOG_DIR}/eukfinder_output/eukfinder_%A_%a.err

# Set environment variable for script location
export PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}"

# Run the main EukFinder script
bash "${SCRIPT_DIR}/10_eukfinder.sh"
EOF

    # Make log directory
    mkdir -p "${LOG_DIR}/eukfinder_output"

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
        echo "  ${LOG_DIR}/eukfinder_output/eukfinder_${JOB_ID}_*.out"
        echo "  ${LOG_DIR}/eukfinder_output/eukfinder_${JOB_ID}_*.err"
        echo ""
        echo "Results will be saved to:"
        echo "  ${OUTPUT_DIR}/eukfinder_output/"
        echo ""
    else
        echo "ERROR: Failed to submit job"
        exit 1
    fi
else
    echo "WARNING: sbatch not found. Not in a SLURM environment."
    echo ""
    echo "To run EukFinder manually, use:"
    echo "  bash ${SCRIPT_DIR}/10_eukfinder.sh"
    echo ""
    echo "Or submit with sbatch after modifying the array size in 10_eukfinder.sh:"
    echo "  #SBATCH --array=0-${MAX_ARRAY_INDEX}%10"
    echo ""
fi

echo "====================================================================="
echo "EukFinder submission complete"
echo "====================================================================="
