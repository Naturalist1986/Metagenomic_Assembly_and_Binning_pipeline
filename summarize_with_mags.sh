#!/bin/bash
# summarize_with_mags.sh - Merge MAG mapping stats with bacterial/non-bacterial summary
# Creates a combined TSV with MAGs as a subset of bacterial reads

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Creating Combined Summary with MAGs"
echo "====================================================================="
echo ""

BACT_MAPPING_DIR="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping"
MAG_MAPPING_DIR="${OUTPUT_DIR}/selected_bins_mapping"
INPUT_TSV="${BACT_MAPPING_DIR}/bacterial_vs_nonbacterial_summary.tsv"
OUTPUT_TSV="${BACT_MAPPING_DIR}/bacterial_vs_nonbacterial_with_mags_summary.tsv"

# Check prerequisites
if [ ! -f "$INPUT_TSV" ]; then
    echo "ERROR: Bacterial vs non-bacterial summary not found: $INPUT_TSV"
    echo "Please run summarize_bacterial_vs_nonbacterial.sh first"
    exit 1
fi

if [ ! -d "$MAG_MAPPING_DIR" ]; then
    echo "ERROR: MAG mapping directory not found: $MAG_MAPPING_DIR"
    echo "Please run submit_selected_bins_mapping.sh first"
    exit 1
fi

echo "Reading bacterial/non-bacterial summary from: $INPUT_TSV"
echo "Reading MAG mapping results from: $MAG_MAPPING_DIR"
echo ""

# Create output TSV header
cat > "$OUTPUT_TSV" << 'EOF'
Treatment	Total_Reads	MAG_Reads	MAG_%	Other_Bacterial_Reads	Other_Bacterial_%	Non-Bacterial_Reads	Non-Bacterial_%	Unmapped_Reads	Unmapped_%
EOF

# Process each treatment from the existing summary
TREATMENT_COUNT=0

# Read the input TSV line by line (skip header)
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r treatment total_reads bacterial_reads bacterial_pct nonbacterial_reads nonbacterial_pct unmapped_reads unmapped_pct; do

    echo "Processing: $treatment"

    # Find MAG mapping summary for this treatment
    mag_summary="${MAG_MAPPING_DIR}/${treatment}/mag_mapping_summary.txt"

    if [ ! -f "$mag_summary" ]; then
        echo "  WARNING: No MAG mapping summary found for $treatment"
        echo "  Using MAG_Reads=0"
        mag_mapped=0
    else
        # Extract MAG-mapped reads from summary file
        mag_mapped=$(grep "Mapped reads:" "$mag_summary" | awk '{print $3}')
        mag_mapped=${mag_mapped:-0}
    fi

    # Calculate Other Bacterial = Bacterial - MAGs
    other_bacterial=$((bacterial_reads - mag_mapped))

    # Guard against negative values (shouldn't happen, but be safe)
    if [ "$other_bacterial" -lt 0 ]; then
        echo "  WARNING: MAG reads ($mag_mapped) exceed bacterial reads ($bacterial_reads)"
        echo "  Setting Other_Bacterial to 0"
        other_bacterial=0
    fi

    # Calculate percentages
    if [ "$total_reads" -gt 0 ]; then
        mag_pct=$(awk "BEGIN {printf \"%.2f\", ($mag_mapped / $total_reads) * 100}")
        other_bact_pct=$(awk "BEGIN {printf \"%.2f\", ($other_bacterial / $total_reads) * 100}")
    else
        mag_pct="0.00"
        other_bact_pct="0.00"
    fi

    # Write to output TSV
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$treatment" "$total_reads" "$mag_mapped" "$mag_pct" \
        "$other_bacterial" "$other_bact_pct" \
        "$nonbacterial_reads" "$nonbacterial_pct" \
        "$unmapped_reads" "$unmapped_pct" >> "$OUTPUT_TSV"

    echo "  MAG: $mag_mapped (${mag_pct}%), Other Bacterial: $other_bacterial (${other_bact_pct}%)"
done

echo ""
echo "====================================================================="
echo "Summary Complete"
echo "====================================================================="
echo ""
echo "Output file:"
echo "  $OUTPUT_TSV"
echo ""

# Display the output
echo "Contents:"
column -t -s $'\t' "$OUTPUT_TSV" 2>/dev/null || cat "$OUTPUT_TSV"
echo ""
echo "To plot, copy this TSV to your local machine and run:"
echo "  python plot_bacterial_vs_nonbacterial.py"
echo ""
