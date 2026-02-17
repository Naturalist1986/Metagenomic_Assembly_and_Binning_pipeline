#!/bin/bash
# summarize_with_mags.sh - Create summary TSV from competitive MAG mapping results
# Reads category_counts.txt from each treatment's mapping directory

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Creating Competitive MAG Mapping Summary"
echo "====================================================================="
echo ""

MAG_MAPPING_DIR="${OUTPUT_DIR}/selected_bins_mapping"
OUTPUT_TSV="${MAG_MAPPING_DIR}/bacterial_vs_nonbacterial_with_mags_summary.tsv"

if [ ! -d "$MAG_MAPPING_DIR" ]; then
    echo "ERROR: MAG mapping directory not found: $MAG_MAPPING_DIR"
    echo "Please run submit_selected_bins_mapping.sh first"
    exit 1
fi

echo "Reading mapping results from: $MAG_MAPPING_DIR"
echo ""

# Create TSV header
cat > "$OUTPUT_TSV" << 'EOF'
Treatment	Total_Reads	MAG_Reads	MAG_%	Other_Bacterial_Reads	Other_Bacterial_%	Non-Bacterial_Reads	Non-Bacterial_%	Unmapped_Reads	Unmapped_%
EOF

TREATMENT_COUNT=0

for treatment_dir in "$MAG_MAPPING_DIR"/*/; do
    [ -d "$treatment_dir" ] || continue

    treatment=$(basename "$treatment_dir")
    counts_file="${treatment_dir}/category_counts.txt"

    if [ ! -f "$counts_file" ]; then
        echo "WARNING: No category_counts.txt for $treatment - skipping"
        continue
    fi

    # Source the counts (total_reads, mag_reads, bact_reads, nonbact_reads, unmapped_reads)
    source "$counts_file"

    # Calculate percentages
    if [ "$total_reads" -gt 0 ]; then
        mag_pct=$(awk "BEGIN {printf \"%.2f\", ($mag_reads / $total_reads) * 100}")
        bact_pct=$(awk "BEGIN {printf \"%.2f\", ($bact_reads / $total_reads) * 100}")
        nonbact_pct=$(awk "BEGIN {printf \"%.2f\", ($nonbact_reads / $total_reads) * 100}")
        unmapped_pct=$(awk "BEGIN {printf \"%.2f\", ($unmapped_reads / $total_reads) * 100}")
    else
        mag_pct="0.00"
        bact_pct="0.00"
        nonbact_pct="0.00"
        unmapped_pct="0.00"
    fi

    # Write to TSV
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$treatment" "$total_reads" "$mag_reads" "$mag_pct" \
        "$bact_reads" "$bact_pct" \
        "$nonbact_reads" "$nonbact_pct" \
        "$unmapped_reads" "$unmapped_pct" >> "$OUTPUT_TSV"

    TREATMENT_COUNT=$((TREATMENT_COUNT + 1))
    echo "  $treatment: MAG=${mag_pct}%  OtherBact=${bact_pct}%  NonBact=${nonbact_pct}%  Unmapped=${unmapped_pct}%"
done

if [ $TREATMENT_COUNT -eq 0 ]; then
    echo "ERROR: No completed mapping results found"
    echo "Make sure mapping jobs have completed successfully"
    exit 1
fi

echo ""
echo "====================================================================="
echo "Summary Complete ($TREATMENT_COUNT treatments)"
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
