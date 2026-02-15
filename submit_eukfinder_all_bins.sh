#!/bin/bash
# submit_eukfinder_all_bins.sh - Identify all bins and submit EukFinder jobs

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "EukFinder Submission Workflow - ALL BINS"
echo "====================================================================="
echo ""

# Step 1: Identify all bins
echo "Step 1: Identifying ALL bins..."
echo "-------------------------------------------------------------------"

ALL_BINS_FILE="${OUTPUT_DIR}/all_bins_list.txt"

"${SCRIPT_DIR}/identify_all_bins.sh" "$ALL_BINS_FILE"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to identify bins"
    exit 1
fi

# Check if any bins were found
TOTAL_BINS=$(grep -v "^#" "$ALL_BINS_FILE" | wc -l)

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
#SBATCH --job-name=eukfinder_all
#SBATCH --array=0-${MAX_ARRAY_INDEX}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --account=ofinkel
#SBATCH --output=${LOG_DIR}/eukfinder_output/eukfinder_all_%A_%a.out
#SBATCH --error=${LOG_DIR}/eukfinder_output/eukfinder_all_%A_%a.err

# Set environment variable for script location
export PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}"

# Use all bins list instead of largest bins list
export BINS_LIST_FILE="${ALL_BINS_FILE}"

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
        echo "  ${LOG_DIR}/eukfinder_output/eukfinder_all_${JOB_ID}_*.out"
        echo "  ${LOG_DIR}/eukfinder_output/eukfinder_all_${JOB_ID}_*.err"
        echo ""
        echo "Results will be saved to:"
        echo "  ${OUTPUT_DIR}/eukfinder_output/"
        echo ""
        echo "After jobs complete:"
        echo "  1. Generate EukFinder summary:"
        echo "     ./create_eukfinder_summary_table.sh"
        echo ""
        echo "  2. Map reads to non-prokaryotic sequences:"
        echo "     ./submit_nonprokaryotic_mapping.sh"
        echo ""
    else
        echo "ERROR: Failed to submit job"
        exit 1
    fi
else
    echo "WARNING: sbatch not found. Not in a SLURM environment."
    echo ""
    echo "To run EukFinder manually, use:"
    echo "  export BINS_LIST_FILE='${ALL_BINS_FILE}'"
    echo "  bash ${SCRIPT_DIR}/10_eukfinder.sh"
    echo ""
    echo "Or submit with sbatch after modifying the array size in 10_eukfinder.sh:"
    echo "  #SBATCH --array=0-${MAX_ARRAY_INDEX}"
    echo ""
fi

echo "====================================================================="
echo "EukFinder submission complete"
echo "====================================================================="
