#!/bin/bash
# summarize_bacterial_vs_nonbacterial.sh - Create comparison table for bacterial vs non-bacterial mapping

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Creating Bacterial vs Non-Bacterial Comparison Table"
echo "====================================================================="
echo ""

MAPPING_DIR="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping"
OUTPUT_TSV="${MAPPING_DIR}/bacterial_vs_nonbacterial_summary.tsv"
OUTPUT_TXT="${MAPPING_DIR}/bacterial_vs_nonbacterial_summary.txt"

if [ ! -d "$MAPPING_DIR" ]; then
    echo "ERROR: Mapping directory not found: $MAPPING_DIR"
    echo "Please run submit_bacterial_vs_nonbacterial_mapping.sh first"
    exit 1
fi

echo "Processing mapping results from: $MAPPING_DIR"
echo ""

# Create TSV header
cat > "$OUTPUT_TSV" << 'EOF'
Treatment	Total_Reads	Bacterial_Reads	Bacterial_%	Non-Bacterial_Reads	Non-Bacterial_%	Unmapped_Reads	Unmapped_%
EOF

# Process each treatment
TREATMENT_COUNT=0

for treatment_dir in "$MAPPING_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    summary_file="${treatment_dir}/mapping_summary.txt"

    if [ ! -f "$summary_file" ]; then
        echo "WARNING: No summary file found for $treatment"
        continue
    fi

    echo "Processing: $treatment"

    # Extract values from summary file
    total_reads=$(grep "Total Reads:" "$summary_file" | awk '{print $3}')
    bacterial_reads=$(grep -A1 "Bacterial (Bact):" "$summary_file" | grep "Mapped reads:" | awk '{print $3}')
    bacterial_pct=$(grep -A2 "Bacterial (Bact):" "$summary_file" | grep "Percentage:" | awk '{print $2}' | tr -d '%')
    nonbacterial_reads=$(grep -A1 "Non-Bacterial" "$summary_file" | grep "Mapped reads:" | awk '{print $3}')
    nonbacterial_pct=$(grep -A2 "Non-Bacterial" "$summary_file" | grep "Percentage:" | awk '{print $2}' | tr -d '%')
    unmapped_reads=$(grep -A1 "Unmapped:" "$summary_file" | grep "Unmapped reads:" | awk '{print $3}')
    unmapped_pct=$(grep -A2 "Unmapped:" "$summary_file" | grep "Percentage:" | awk '{print $2}' | tr -d '%')

    # Add to TSV
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$treatment" "$total_reads" "$bacterial_reads" "$bacterial_pct" \
        "$nonbacterial_reads" "$nonbacterial_pct" "$unmapped_reads" "$unmapped_pct" >> "$OUTPUT_TSV"

    ((TREATMENT_COUNT++))
done

if [ $TREATMENT_COUNT -eq 0 ]; then
    echo "ERROR: No completed mapping results found"
    echo "Make sure mapping jobs have completed successfully"
    exit 1
fi

echo ""
echo "Processed $TREATMENT_COUNT treatments"
echo ""

# Create formatted text summary
cat > "$OUTPUT_TXT" << EOF
Bacterial vs Non-Bacterial Mapping Summary
==========================================
Generated: $(date)
Total treatments: $TREATMENT_COUNT

EOF

# Add formatted table
awk -F'\t' 'NR==1 {
    printf "%-20s %15s %15s %12s %18s %15s %15s %12s\n", $1, $2, $3, $4, $5, $6, $7, $8
    printf "%-20s %15s %15s %12s %18s %15s %15s %12s\n",
        "--------------------", "---------------", "---------------", "------------",
        "------------------", "---------------", "---------------", "------------"
    next
}
{
    printf "%-20s %15s %15s %11s%% %18s %14s%% %15s %11s%%\n",
        $1, $2, $3, $4, $5, $6, $7, $8
}' "$OUTPUT_TSV" >> "$OUTPUT_TXT"

# Add overall statistics
echo "" >> "$OUTPUT_TXT"
echo "Overall Statistics:" >> "$OUTPUT_TXT"
echo "------------------" >> "$OUTPUT_TXT"

awk -F'\t' 'NR>1 {
    total_reads += $2
    bacterial_reads += $3
    nonbacterial_reads += $5
    unmapped_reads += $7
}
END {
    bacterial_pct = (bacterial_reads / total_reads) * 100
    nonbacterial_pct = (nonbacterial_reads / total_reads) * 100
    unmapped_pct = (unmapped_reads / total_reads) * 100

    printf "Total reads across all treatments: %d\n", total_reads
    printf "Bacterial reads: %d (%.2f%%)\n", bacterial_reads, bacterial_pct
    printf "Non-bacterial reads: %d (%.2f%%)\n", nonbacterial_reads, nonbacterial_pct
    printf "Unmapped reads: %d (%.2f%%)\n", unmapped_reads, unmapped_pct
}' "$OUTPUT_TSV" >> "$OUTPUT_TXT"

echo ""
echo "====================================================================="
echo "Summary Complete"
echo "====================================================================="
echo ""
echo "Output files:"
echo "  TSV: $OUTPUT_TSV"
echo "  TXT: $OUTPUT_TXT"
echo ""

# Display the summary
cat "$OUTPUT_TXT"

echo ""
echo "To view in spreadsheet format:"
echo "  column -t -s $'\\t' $OUTPUT_TSV | less -S"
echo ""
echo "To import into Excel/Google Sheets:"
echo "  Open: $OUTPUT_TSV"
echo ""
