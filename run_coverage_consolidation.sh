#!/bin/bash
# Wrapper script to run coverage consolidation with proper configuration

# Set the actual paths for your pipeline
export OUTPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly"
export INPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene"

# Set script directory for sourcing config
export PIPELINE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "=========================================="
echo "Coverage Consolidation Wrapper"
echo "=========================================="
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "INPUT_DIR: $INPUT_DIR"
echo "SCRIPT_DIR: $PIPELINE_SCRIPT_DIR"
echo ""
echo "This script will consolidate:"
echo "  - All bins from all treatments"
echo "  - Coverage data with unmapped reads"
echo "  - GTDB-Tk taxonomy"
echo "  - Binette quality metrics"
echo ""
echo "Output will be saved to:"
echo "  - $OUTPUT_DIR/consolidated_bins/"
echo "  - $OUTPUT_DIR/consolidated_coverage_data/"
echo ""

# Check if this should be submitted to SLURM or run directly
if [ "$1" == "--submit" ]; then
    echo "Submitting to SLURM..."
    sbatch "$PIPELINE_SCRIPT_DIR/12_coverage_consolidation.sh"
elif [ "$1" == "--interactive" ] || [ "$1" == "-i" ]; then
    echo "Running interactively..."
    echo ""
    bash "$PIPELINE_SCRIPT_DIR/12_coverage_consolidation.sh"
else
    echo "Usage:"
    echo "  $0 --submit       # Submit to SLURM queue"
    echo "  $0 --interactive  # Run directly (interactive mode)"
    echo ""
    echo "Example:"
    echo "  $0 --interactive"
    exit 1
fi
