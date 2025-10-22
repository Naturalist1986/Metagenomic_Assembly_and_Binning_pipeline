#!/bin/bash
# debug_sample_parsing.sh - Debug the sample parsing step by step

export INPUT_DIR="/sci/backup/aerez/aerez/moshea/Pipeline_Test/raw_fastqs/"
export OUTPUT_DIR="/sci/backup/aerez/aerez/moshea/Pipeline_Test"
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"
export SAMPLE_SHEET="/sci/backup/aerez/aerez/moshea/Pipeline_Test/sample_sheet.xlsx"

echo "=== DEBUG SAMPLE PARSING ==="
echo "INPUT_DIR: $INPUT_DIR"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "WORK_DIR: $WORK_DIR"
echo "SAMPLE_SHEET: $SAMPLE_SHEET"

# Check if sample sheet exists
echo "Sample sheet exists: $([ -f "$SAMPLE_SHEET" ] && echo "YES" || echo "NO")"
if [ -f "$SAMPLE_SHEET" ]; then
    echo "Sample sheet size: $(ls -la "$SAMPLE_SHEET")"
fi

# Source the utilities
echo "Sourcing config utilities..."
source ./00_config_utilities.sh

# Create work directory manually
echo "Creating work directory..."
mkdir -p "$WORK_DIR"
echo "Work dir created: $([ -d "$WORK_DIR" ] && echo "YES" || echo "NO")"

# Test the init_sample_storage function
echo "Testing init_sample_storage..."
init_sample_storage
echo "Sample info file after init: $([ -f "${WORK_DIR}/sample_info.txt" ] && echo "EXISTS" || echo "MISSING")"

# Test parsing manually
echo "Testing parse_sample_sheet..."
parse_sample_sheet "$SAMPLE_SHEET"

# Check results
echo "=== RESULTS ==="
echo "Sample info file exists: $([ -f "${WORK_DIR}/sample_info.txt" ] && echo "YES" || echo "NO")"
if [ -f "${WORK_DIR}/sample_info.txt" ]; then
    echo "Sample info file content:"
    cat "${WORK_DIR}/sample_info.txt"
    echo "Total samples: $(get_total_samples)"
else
    echo "Sample info file was not created!"
    echo "Files in work dir:"
    ls -la "$WORK_DIR/"
fi

echo "=== TESTING SAMPLE RETRIEVAL ==="
for i in {0..3}; do
    echo "Testing index $i:"
    sample_info=$(get_sample_info_by_index $i 2>/dev/null)
    if [ -n "$sample_info" ]; then
        echo "  Found: $sample_info"
    else
        echo "  Not found"
    fi
done
