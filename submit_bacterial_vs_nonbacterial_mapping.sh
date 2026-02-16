#!/bin/bash
# submit_bacterial_vs_nonbacterial_mapping.sh - Submit bacterial vs non-bacterial mapping workflow

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Bacterial vs Non-Bacterial Mapping Workflow (Per-Binner)"
echo "====================================================================="
echo ""
echo "This workflow maps reads to sequences from EACH BINNER separately,"
echo "then averages results to avoid duplicate contig counting."
echo ""

# Step 1: Check EukFinder results
echo "Step 1: Checking EukFinder results..."
echo "-------------------------------------------------------------------"

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder_output"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder directory not found: $EUKFINDER_DIR"
    echo "Please run EukFinder first:"
    echo "  ./submit_eukfinder_all_bins.sh"
    exit 1
fi

# Count treatments
NUM_TREATMENTS=$(ls -1d "${EUKFINDER_DIR}"/*/ 2>/dev/null | wc -l)

if [ "$NUM_TREATMENTS" -eq 0 ]; then
    echo "ERROR: No treatments found in EukFinder results"
    echo "Please run EukFinder first:"
    echo "  ./submit_eukfinder_all_bins.sh"
    exit 1
fi

echo "Found $NUM_TREATMENTS treatments in EukFinder results"
echo ""

# Show which binners are available
echo "Checking available binners..."
for treatment_dir in "${EUKFINDER_DIR}"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    treatment=$(basename "$treatment_dir")
    binners=$(ls -1 "$treatment_dir" | grep -oE '_[^_]+_' | tr -d '_' | sort -u | tr '\n' ' ')
    echo "  $treatment: $binners"
done
echo ""

# Step 2: Submit mapping jobs
echo "Step 2: Submitting SLURM array job for per-binner mapping..."
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
#SBATCH --array=0-${MAX_ARRAY_INDEX}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
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
