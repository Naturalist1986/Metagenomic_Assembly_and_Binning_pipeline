#!/bin/bash
# update_eukfinder_taxonomy.sh - Update EukFinder NCBI taxonomy database once before parallel runs

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "EukFinder Taxonomy Database Update"
echo "====================================================================="
echo ""
echo "This script updates the NCBI taxonomy database used by EukFinder."
echo "Run this ONCE before submitting parallel EukFinder jobs to avoid"
echo "database locking issues."
echo ""

# Initialize conda
init_conda

# Activate eukfinder environment
echo "Activating eukfinder conda environment..."
activate_env eukfinder

# Check if eukfinder is available
if ! command -v eukfinder &> /dev/null; then
    echo "ERROR: eukfinder command not found"
    echo "Please ensure the eukfinder conda environment is properly installed"
    exit 1
fi

echo "EukFinder version:"
eukfinder --help | head -5
echo ""

# Create a test directory
TEST_DIR=$(mktemp -d -p "${WORK_DIR}" eukfinder_taxonomy_update_XXXXXX)
echo "Created test directory: $TEST_DIR"
cd "$TEST_DIR"

# Create a minimal test FASTA file
cat > test.fasta << 'EOF'
>test_seq_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
EOF

echo "Running EukFinder with taxonomy update enabled..."
echo "This will download and update the NCBI taxonomy database."
echo "This may take several minutes..."
echo ""

# Run eukfinder with a minimal file and taxonomy update enabled
# This will trigger the taxonomy database update
eukfinder long_seqs \
    -l test.fasta \
    -o test_output \
    --mhlen 100 \
    --cdb "$EUKFINDER_CENTRIFUGE_DB" \
    -n 2 \
    -z 1 \
    -t True \
    -p "$EUKFINDER_PLAST_DB" \
    -m "$EUKFINDER_PLAST_ID_MAP" \
    -e 0.01 \
    --pid 60 \
    --cov 30 \
    2>&1 | tee taxonomy_update.log

EXIT_CODE=$?

# Cleanup
cd ..
rm -rf "$TEST_DIR"

deactivate_env

if [ $EXIT_CODE -eq 0 ] || grep -q "generating entries" taxonomy_update.log; then
    echo ""
    echo "====================================================================="
    echo "Taxonomy database update completed successfully!"
    echo "====================================================================="
    echo ""
    echo "The NCBI taxonomy database has been updated in:"
    echo "  ~/.etetoolkit/taxa.sqlite"
    echo ""
    echo "You can now run parallel EukFinder jobs without database locking issues:"
    echo "  ./submit_eukfinder.sh"
    echo ""
else
    echo ""
    echo "====================================================================="
    echo "WARNING: Taxonomy database update may have failed"
    echo "====================================================================="
    echo ""
    echo "Exit code: $EXIT_CODE"
    echo "Check the log for details"
    echo ""
fi
