#!/bin/bash
# 11_collect_bins_for_gtdbtk.sh - Collect all non-Soil bins for GTDB-Tk classification

BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
COVERM_DIR="${BASE_DIR}/coverm_treatment_level"
OUTPUT_DIR="${BASE_DIR}/all_bins_for_gtdbtk"
SUMMARY_FILE="${OUTPUT_DIR}/bin_collection_summary.txt"

echo "========================================="
echo "Collecting Bins for GTDB-Tk"
echo "========================================="
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Initialize summary
cat > "$SUMMARY_FILE" << EOF
Bin Collection Summary for GTDB-Tk Classification
==================================================

Date: $(date)
Source: ${COVERM_DIR}
Output: ${OUTPUT_DIR}

Bins collected (excluding Soil bins):
======================================

EOF

# Counters
total_bins=0
soil_bins_excluded=0
treatments_processed=0

# Process each treatment
for treatment_dir in "$COVERM_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    
    treatment=$(basename "$treatment_dir")
    bins_dir="${treatment_dir}/all_bins"
    
    if [ ! -d "$bins_dir" ]; then
        echo "WARNING: No all_bins directory found for $treatment"
        continue
    fi
    
    ((treatments_processed++))
    echo "Processing treatment: $treatment"
    
    local_count=0
    local_excluded=0
    
    # Copy bins, excluding Soil bins
    for bin_file in "$bins_dir"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        bin_name=$(basename "$bin_file")
        
        # Check if it's a Soil bin
        if [[ "$bin_name" =~ _Soil\. ]]; then
            ((local_excluded++))
            ((soil_bins_excluded++))
            continue
        fi
        
        # Copy bin to output directory
        cp "$bin_file" "$OUTPUT_DIR/"
        ((local_count++))
        ((total_bins++))
    done
    
    echo "  ✓ Collected $local_count bins (excluded $local_excluded Soil bins)"
    echo "$treatment: $local_count bins" >> "$SUMMARY_FILE"
done

echo ""
echo "========================================="
echo "Collection Complete"
echo "========================================="
echo ""
echo "Summary:"
echo "  Treatments processed: $treatments_processed"
echo "  Total bins collected: $total_bins"
echo "  Soil bins excluded: $soil_bins_excluded"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""

# Add summary statistics
cat >> "$SUMMARY_FILE" << EOF

Summary Statistics:
===================
Total treatments: $treatments_processed
Total bins collected: $total_bins
Soil bins excluded: $soil_bins_excluded

Output directory: $OUTPUT_DIR

Bins by treatment:
==================
EOF

# Count bins per treatment in output
for treatment in carK carR ces hok mtz RH; do
    count=$(ls -1 "$OUTPUT_DIR"/${treatment}_*.fa 2>/dev/null | wc -l)
    if [ $count -gt 0 ]; then
        echo "  $treatment: $count bins" >> "$SUMMARY_FILE"
    fi
done

echo "Summary file created: $SUMMARY_FILE"
echo ""

# List some examples
echo "Example bins collected:"
ls "$OUTPUT_DIR"/*.fa | head -10 | while read bin; do
    bin_name=$(basename "$bin")
    bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
    contigs=$(grep -c "^>" "$bin")
    printf "  %-50s %10s bp, %4d contigs\n" "$bin_name" "$bin_size" "$contigs"
done
echo ""

# Check bin quality
echo "Bin size distribution:"
echo "  Total bins: $total_bins"
echo "  Large bins (>500kb): $(find "$OUTPUT_DIR" -name "*.fa" -exec sh -c 'size=$(grep -v "^>" "$1" | tr -d "\n" | wc -c); [ $size -gt 500000 ] && echo 1' _ {} \; | wc -l)"
echo "  Medium bins (100-500kb): $(find "$OUTPUT_DIR" -name "*.fa" -exec sh -c 'size=$(grep -v "^>" "$1" | tr -d "\n" | wc -c); [ $size -ge 100000 ] && [ $size -le 500000 ] && echo 1' _ {} \; | wc -l)"
echo "  Small bins (<100kb): $(find "$OUTPUT_DIR" -name "*.fa" -exec sh -c 'size=$(grep -v "^>" "$1" | tr -d "\n" | wc -c); [ $size -lt 100000 ] && echo 1' _ {} \; | wc -l)"
echo ""

echo "========================================="
echo "Ready for GTDB-Tk Classification"
echo "========================================="
echo ""
echo "Next steps:"
echo "  1. Review collected bins: ls -lh $OUTPUT_DIR"
echo "  2. Run GTDB-Tk: ./12_run_gtdbtk.sh"
echo ""
