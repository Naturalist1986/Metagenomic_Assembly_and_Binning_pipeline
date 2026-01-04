#!/bin/bash
# summarize_nonprokaryotic_mapping.sh - Create summary of non-prokaryotic read mapping

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Creating Non-Prokaryotic Mapping Summary"
echo "====================================================================="
echo ""

MAPPING_DIR="${OUTPUT_DIR}/nonprokaryotic_mapping"

if [ ! -d "$MAPPING_DIR" ]; then
    echo "ERROR: Mapping directory not found: $MAPPING_DIR"
    echo "Please run submit_nonprokaryotic_mapping.sh first"
    exit 1
fi

# Output files
SUMMARY_TSV="${MAPPING_DIR}/nonprokaryotic_mapping_summary.tsv"
SUMMARY_REPORT="${MAPPING_DIR}/nonprokaryotic_mapping_report.txt"

echo "Collecting mapping results from: $MAPPING_DIR"
echo ""

# Create TSV header
cat > "$SUMMARY_TSV" << 'EOF'
Treatment	Total_Read_Pairs	Mapped_Read_Pairs	Percent_Mapped	Nonprok_Sequences	Nonprok_Size_bp
EOF

# Collect results from each treatment
TOTAL_TREATMENTS=0
TOTAL_READS=0
TOTAL_MAPPED=0

for treatment_dir in "$MAPPING_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    summary_file="${treatment_dir}/mapping_summary.tsv"

    if [ ! -f "$summary_file" ]; then
        echo "WARNING: No mapping summary found for $treatment"
        continue
    fi

    echo "Processing: $treatment"

    # Read the summary TSV
    while IFS=$'\t' read -r treat total_reads mapped_reads percent_mapped; do
        # Get non-prokaryotic sequence info
        nonprok_fasta="${OUTPUT_DIR}/nonprokaryotic_sequences/${treatment}_nonprokaryotic.fasta"

        if [ -f "$nonprok_fasta" ]; then
            num_seqs=$(grep -c "^>" "$nonprok_fasta" 2>/dev/null || echo "0")
            size_bp=$(grep -v "^>" "$nonprok_fasta" | tr -d '\n' | wc -c 2>/dev/null || echo "0")
        else
            num_seqs=0
            size_bp=0
        fi

        # Write to summary TSV
        printf "%s\t%d\t%d\t%.4f\t%d\t%d\n" \
            "$treatment" "$total_reads" "$mapped_reads" "$percent_mapped" "$num_seqs" "$size_bp" >> "$SUMMARY_TSV"

        # Accumulate totals
        ((TOTAL_TREATMENTS++))
        TOTAL_READS=$((TOTAL_READS + total_reads))
        TOTAL_MAPPED=$((TOTAL_MAPPED + mapped_reads))

    done < "$summary_file"
done

# Calculate overall percentage
if [ $TOTAL_READS -gt 0 ]; then
    OVERALL_PERCENT=$(awk "BEGIN {printf \"%.4f\", ($TOTAL_MAPPED / $TOTAL_READS) * 100}")
else
    OVERALL_PERCENT=0.0000
fi

# Create report
cat > "$SUMMARY_REPORT" << EOF
Non-Prokaryotic Read Mapping Summary Report
============================================
Generated: $(date)

Total treatments analyzed: $TOTAL_TREATMENTS

Overall Statistics:
-------------------
Total read pairs across all treatments: $TOTAL_READS
Total read pairs mapped to non-prokaryotic sequences: $TOTAL_MAPPED
Overall percentage mapped: $OVERALL_PERCENT%

Per-Treatment Results:
----------------------
EOF

# Add formatted table to report
awk -F'\t' 'NR>1 {
    printf "%-15s %15s %15s %12s %12s %15s\n", $1, $2, $3, $4"%", $5, $6
}' "$SUMMARY_TSV" | sort -t$'\t' -k4 -nr | \
awk 'BEGIN {
    printf "%-15s %15s %15s %12s %12s %15s\n", "Treatment", "Total Reads", "Mapped Reads", "% Mapped", "Nonprok Seqs", "Nonprok Size"
    printf "%-15s %15s %15s %12s %12s %15s\n", "---------", "-----------", "------------", "--------", "------------", "------------"
}
{print}' >> "$SUMMARY_REPORT"

echo "" >> "$SUMMARY_REPORT"
echo "Top 5 Treatments by Mapping Percentage:" >> "$SUMMARY_REPORT"
echo "----------------------------------------" >> "$SUMMARY_REPORT"

awk -F'\t' 'NR>1 {print $1"\t"$4}' "$SUMMARY_TSV" | \
    sort -t$'\t' -k2 -nr | \
    head -5 | \
    awk -F'\t' 'BEGIN {printf "%-15s %12s\n", "Treatment", "% Mapped"; printf "%-15s %12s\n", "---------", "--------"}
                {printf "%-15s %11.4f%%\n", $1, $2}' >> "$SUMMARY_REPORT"

echo "" >> "$SUMMARY_REPORT"
echo "Detailed Results by Treatment:" >> "$SUMMARY_REPORT"
echo "------------------------------" >> "$SUMMARY_REPORT"

for treatment_dir in "$MAPPING_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    stats_file="${treatment_dir}/mapping_stats.txt"

    if [ -f "$stats_file" ]; then
        echo "" >> "$SUMMARY_REPORT"
        echo "=== $treatment ===" >> "$SUMMARY_REPORT"
        tail -n +3 "$stats_file" | head -n 20 >> "$SUMMARY_REPORT"
    fi
done

# Display the report
echo ""
echo "====================================================================="
echo "Summary Complete"
echo "====================================================================="
echo ""
cat "$SUMMARY_REPORT"
echo ""
echo "Files created:"
echo "  Summary table (TSV): $SUMMARY_TSV"
echo "  Detailed report: $SUMMARY_REPORT"
echo ""
echo "To view the TSV in formatted columns:"
echo "  column -t -s $'\\t' $SUMMARY_TSV | less -S"
echo ""
echo "To import into Excel/Google Sheets:"
echo "  Open: $SUMMARY_TSV"
echo ""
