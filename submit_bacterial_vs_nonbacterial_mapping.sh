#!/bin/bash
# submit_bacterial_vs_nonbacterial_mapping.sh - Submit bacterial vs non-bacterial mapping workflow

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Bacterial vs Non-Bacterial Mapping Workflow"
echo "====================================================================="
echo ""

# Step 1: Collect bacterial sequences
echo "Step 1: Collecting bacterial sequences from EukFinder results..."
echo "-------------------------------------------------------------------"

"${SCRIPT_DIR}/collect_bacterial_sequences.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to collect bacterial sequences"
    exit 1
fi

# Step 2: Collect non-bacterial sequences
echo ""
echo "Step 2: Collecting non-bacterial sequences from EukFinder results..."
echo "-------------------------------------------------------------------"

"${SCRIPT_DIR}/collect_nonbacterial_sequences.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to collect non-bacterial sequences"
    exit 1
fi

# Check how many treatments were processed
BACTERIAL_DIR="${OUTPUT_DIR}/bacterial_sequences"
NUM_TREATMENTS=$(ls -1 "${BACTERIAL_DIR}"/*_bacterial.fasta 2>/dev/null | wc -l)

if [ "$NUM_TREATMENTS" -eq 0 ]; then
    echo "ERROR: No sequence files found"
    exit 1
fi

echo ""
echo "Found $NUM_TREATMENTS treatments with sequences"
echo ""

# Step 3: Submit mapping jobs
echo "Step 3: Submitting SLURM array job for read mapping..."
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
#SBATCH --job-name=bact_map
#SBATCH --array=0-${MAX_ARRAY_INDEX}%5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel
#SBATCH --output=${LOG_DIR}/bacterial_vs_nonbacterial_mapping/map_%A_%a.out
#SBATCH --error=${LOG_DIR}/bacterial_vs_nonbacterial_mapping/map_%A_%a.err

# Set environment variable for script location
export PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}"

# Run the mapping script
bash "${SCRIPT_DIR}/map_bacterial_vs_nonbacterial.sh"
EOF

    # Make log directory
    mkdir -p "${LOG_DIR}/bacterial_vs_nonbacterial_mapping"

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
        echo "  ${LOG_DIR}/bacterial_vs_nonbacterial_mapping/map_${JOB_ID}_*.out"
        echo "  ${LOG_DIR}/bacterial_vs_nonbacterial_mapping/map_${JOB_ID}_*.err"
        echo ""
        echo "After jobs complete, generate comparison table with:"
        echo "  ./summarize_bacterial_vs_nonbacterial.sh"
        echo ""
    else
        echo "ERROR: Failed to submit job"
        exit 1
    fi
else
    echo "WARNING: sbatch not found. Not in a SLURM environment."
    echo ""
    echo "To run mapping manually, use:"
    echo "  bash ${SCRIPT_DIR}/map_bacterial_vs_nonbacterial.sh"
    echo ""
fi

echo "====================================================================="
echo "Submission complete"
echo "====================================================================="
