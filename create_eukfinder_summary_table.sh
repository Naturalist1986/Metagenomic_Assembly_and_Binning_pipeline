#!/bin/bash
# create_eukfinder_summary_table.sh - Create consolidated EukFinder results table with percentages

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Output file
OUTPUT_FILE="${1:-${OUTPUT_DIR}/eukfinder_consolidated_results.tsv}"

echo "====================================================================="
echo "Creating Consolidated EukFinder Results Table"
echo "====================================================================="
echo ""

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder_output"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder results directory not found: $EUKFINDER_DIR"
    exit 1
fi

echo "Searching for EukFinder results in: $EUKFINDER_DIR"
echo "Output will be saved to: $OUTPUT_FILE"
echo ""

# Create header for TSV file
cat > "$OUTPUT_FILE" << 'EOF'
Treatment	Sample	Binner	Bin	Total_Sequences	Total_Size_bp	Euk_Seqs	Euk_bp	Euk_%	Bact_Seqs	Bact_bp	Bact_%	Arch_Seqs	Arch_bp	Arch_%	Unk_Seqs	Unk_bp	Unk_%	EUnk_Seqs	EUnk_bp	EUnk_%	Misc_Seqs	Misc_bp	Misc_%
EOF

# Find all summary_table.txt files and process them
RESULTS_COUNT=0

for treatment_dir in "$EUKFINDER_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")

    for bin_dir in "$treatment_dir"/*; do
        if [ ! -d "$bin_dir" ]; then
            continue
        fi

        bin_dirname=$(basename "$bin_dir")
        results_dir="${bin_dir}/Eukfinder_results"
        summary_table="${results_dir}/summary_table.txt"

        if [ ! -f "$summary_table" ]; then
            continue
        fi

        echo "Processing: $treatment/$bin_dirname"

        # Parse bin information from directory name
        # Format: <sample>_<binner>_<binname>
        # Example: RH_metabat2_bin.4
        if [[ "$bin_dirname" =~ ^(.+)_([^_]+)_(.+)$ ]]; then
            sample="${BASH_REMATCH[1]}"
            binner="${BASH_REMATCH[2]}"
            bin_name="${BASH_REMATCH[3]}"
        else
            # Fallback if name doesn't match expected pattern
            sample="$bin_dirname"
            binner="unknown"
            bin_name="$bin_dirname"
        fi

        # Initialize variables for all categories
        declare -A seqs
        declare -A sizes
        seqs=(["Euk"]=0 ["Bact"]=0 ["Arch"]=0 ["Unk"]=0 ["EUnk"]=0 ["Misc"]=0)
        sizes=(["Euk"]=0 ["Bact"]=0 ["Arch"]=0 ["Unk"]=0 ["EUnk"]=0 ["Misc"]=0)

        # Parse the summary_table.txt
        while IFS=$'\t' read -r group num_seq total_size; do
            # Skip the header line
            if [[ "$group" == "Group" ]]; then
                continue
            fi

            # Store the values
            if [[ -n "$group" ]] && [[ -n "$num_seq" ]] && [[ -n "$total_size" ]]; then
                seqs["$group"]=$num_seq
                sizes["$group"]=$total_size
            fi
        done < "$summary_table"

        # Calculate totals
        total_seqs=0
        total_size=0
        for category in Euk Bact Arch Unk EUnk Misc; do
            total_seqs=$((total_seqs + ${seqs[$category]}))
            total_size=$((total_size + ${sizes[$category]}))
        done

        # Calculate percentages (avoid division by zero)
        declare -A pct_seqs
        if [ $total_seqs -gt 0 ]; then
            for category in Euk Bact Arch Unk EUnk Misc; do
                pct_seqs[$category]=$(awk "BEGIN {printf \"%.2f\", (${seqs[$category]} / $total_seqs) * 100}")
            done
        else
            for category in Euk Bact Arch Unk EUnk Misc; do
                pct_seqs[$category]="0.00"
            done
        fi

        # Output row to TSV
        # Columns: Treatment, Sample, Binner, Bin, Total_Sequences, Total_Size_bp, ...
        printf "%s\t%s\t%s\t%s\t%d\t%d" "$treatment" "$sample" "$binner" "$bin_name" "$total_seqs" "$total_size" >> "$OUTPUT_FILE"

        for category in Euk Bact Arch Unk EUnk Misc; do
            printf "\t%d\t%d\t%s" "${seqs[$category]}" "${sizes[$category]}" "${pct_seqs[$category]}" >> "$OUTPUT_FILE"
        done

        printf "\n" >> "$OUTPUT_FILE"

        ((RESULTS_COUNT++))
    done
done

echo ""
echo "====================================================================="
echo "Summary Table Creation Complete"
echo "====================================================================="
echo ""
echo "Processed $RESULTS_COUNT result sets"
echo "Results saved to: $OUTPUT_FILE"
echo ""

if [ $RESULTS_COUNT -eq 0 ]; then
    echo "No EukFinder results found to process"
    echo "Expected structure:"
    echo "  ${EUKFINDER_DIR}/<treatment>/<sample>/Eukfinder_results/summary_table.txt"
    exit 1
fi

# Create a human-readable summary report
REPORT_FILE="${OUTPUT_FILE%.tsv}_report.txt"

cat > "$REPORT_FILE" << EOF
EukFinder Results Summary Report
=================================
Generated: $(date)
Total bins analyzed: $RESULTS_COUNT

EOF

echo "Summary by Category:" >> "$REPORT_FILE"
echo "-------------------" >> "$REPORT_FILE"

# Calculate overall statistics
awk -F'\t' 'NR>1 {
    total_seqs += $5
    total_size += $6
    euk_seqs += $7
    euk_bp += $8
    bact_seqs += $10
    bact_bp += $11
    arch_seqs += $13
    arch_bp += $14
    unk_seqs += $16
    unk_bp += $17
    eunk_seqs += $19
    eunk_bp += $20
}
END {
    printf "Total sequences analyzed: %d (%d bp)\n\n", total_seqs, total_size
    printf "%-15s %12s %15s %10s\n", "Category", "Sequences", "Size (bp)", "Percent"
    printf "%-15s %12s %15s %10s\n", "--------", "---------", "--------", "-------"
    if (total_seqs > 0) {
        printf "%-15s %12d %15d %9.2f%%\n", "Eukaryotic", euk_seqs, euk_bp, (euk_seqs/total_seqs)*100
        printf "%-15s %12d %15d %9.2f%%\n", "Bacterial", bact_seqs, bact_bp, (bact_seqs/total_seqs)*100
        printf "%-15s %12d %15d %9.2f%%\n", "Archaeal", arch_seqs, arch_bp, (arch_seqs/total_seqs)*100
        printf "%-15s %12d %15d %9.2f%%\n", "Unknown", unk_seqs, unk_bp, (unk_seqs/total_seqs)*100
        printf "%-15s %12d %15d %9.2f%%\n", "Euk+Unknown", eunk_seqs, eunk_bp, (eunk_seqs/total_seqs)*100
    }
}' "$OUTPUT_FILE" >> "$REPORT_FILE"

echo "" >> "$REPORT_FILE"
echo "Top 10 Bins by Eukaryotic Content (by sequence count):" >> "$REPORT_FILE"
echo "-------------------------------------------------------" >> "$REPORT_FILE"

awk -F'\t' 'NR>1 {print $1"\t"$3"\t"$4"\t"$7"\t"$9"%"}' "$OUTPUT_FILE" | \
    sort -t$'\t' -k4 -nr | \
    head -10 | \
    awk -F'\t' 'BEGIN {printf "%-15s %-15s %-20s %12s %10s\n", "Treatment", "Binner", "Bin", "Euk Seqs", "Euk %"}
                {printf "%-15s %-15s %-20s %12s %10s\n", $1, $2, $3, $4, $5}' >> "$REPORT_FILE"

echo "" >> "$REPORT_FILE"
echo "Results by Treatment:" >> "$REPORT_FILE"
echo "--------------------" >> "$REPORT_FILE"

awk -F'\t' 'NR>1 {
    treatment[$1]++
    euk[$1] += $7
    bact[$1] += $10
    total[$1] += $5
}
END {
    printf "%-15s %10s %12s %12s %12s %10s\n", "Treatment", "Bins", "Total Seqs", "Euk Seqs", "Bact Seqs", "Euk %"
    for (t in treatment) {
        euk_pct = (total[t] > 0) ? (euk[t]/total[t])*100 : 0
        printf "%-15s %10d %12d %12d %12d %9.2f%%\n", t, treatment[t], total[t], euk[t], bact[t], euk_pct
    }
}' "$OUTPUT_FILE" >> "$REPORT_FILE"

# Display the report
echo ""
cat "$REPORT_FILE"
echo ""
echo "Full results table: $OUTPUT_FILE"
echo "Summary report: $REPORT_FILE"
echo ""
echo "To view the full table:"
echo "  column -t -s $'\\t' $OUTPUT_FILE | less -S"
echo ""
echo "To import into Excel/Google Sheets:"
echo "  Open the TSV file: $OUTPUT_FILE"
echo ""
