#!/bin/bash
# identify_largest_bins.sh - Identify the largest bins from each binning tool
# Supports both individual assembly and coassembly folder structures

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Parameters
NUM_LARGEST="${1:-2}"  # Number of largest bins to select per binner (default: 2)
OUTPUT_FILE="${2:-${OUTPUT_DIR}/largest_bins_list.txt}"

echo "Identifying the ${NUM_LARGEST} largest bins from each binning tool..."

# Create output file header
cat > "$OUTPUT_FILE" << 'EOF'
# Largest bins identified for EukFinder analysis
# Format: BIN_PATH|TREATMENT|SAMPLE|BINNER|BIN_NAME|SIZE_BP|NUM_CONTIGS
EOF

# Function to process bins in a binner directory
process_binner_dir() {
    local binner_dir="$1"
    local treatment="$2"
    local sample="$3"
    local binner="$4"

    # Find all bins and calculate sizes
    local tmp_file=$(mktemp)

    for bin in "$binner_dir"/*.fa; do
        if [ ! -f "$bin" ]; then
            continue
        fi

        bin_name=$(basename "$bin")
        # Calculate bin size (total bp)
        bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
        # Count contigs
        num_contigs=$(grep -c "^>" "$bin")

        echo "${bin}|${bin_size}|${num_contigs}" >> "$tmp_file"
    done

    # Sort by size (descending) and take the largest N bins
    if [ -s "$tmp_file" ]; then
        binner_name=$(echo "$binner" | sed 's/_bins$//')

        sort -t'|' -k2 -nr "$tmp_file" | head -n "$NUM_LARGEST" | while IFS='|' read -r bin_path bin_size num_contigs; do
            bin_name=$(basename "$bin_path")
            echo "${bin_path}|${treatment}|${sample}|${binner_name}|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
            echo "    ${binner_name}: ${bin_name} - ${bin_size} bp, ${num_contigs} contigs"
        done
    fi

    rm -f "$tmp_file"
}

# Scan all binning directories
for treatment_dir in "${OUTPUT_DIR}/binning"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"

    # Check if this is a coassembly structure (bins directly under treatment)
    # by checking if metabat2_bins, maxbin2_bins, or concoct_bins exist at treatment level
    is_coassembly=false
    for binner in metabat2_bins maxbin2_bins concoct_bins; do
        if [ -d "${treatment_dir}/${binner}" ]; then
            is_coassembly=true
            break
        fi
    done

    if [ "$is_coassembly" = true ]; then
        # Coassembly structure: bins are directly under treatment
        echo "  Detected coassembly structure for treatment: $treatment"

        for binner in metabat2_bins maxbin2_bins concoct_bins; do
            binner_dir="${treatment_dir}/${binner}"

            if [ ! -d "$binner_dir" ]; then
                continue
            fi

            # For coassembly, use treatment name as sample
            process_binner_dir "$binner_dir" "$treatment" "$treatment" "$binner"
        done
    else
        # Individual assembly structure: bins are under treatment/sample
        echo "  Detected individual assembly structure for treatment: $treatment"

        for sample_dir in "$treatment_dir"/*; do
            if [ ! -d "$sample_dir" ]; then
                continue
            fi

            sample=$(basename "$sample_dir")

            # Skip non-sample directories
            if [[ "$sample" == "work_files" ]] || [[ "$sample" == *"_bins" ]]; then
                continue
            fi

            echo "    Sample: $sample"

            # Check each binner directory
            for binner in metabat2_bins maxbin2_bins concoct_bins; do
                binner_dir="${sample_dir}/${binner}"

                if [ ! -d "$binner_dir" ]; then
                    continue
                fi

                process_binner_dir "$binner_dir" "$treatment" "$sample" "$binner"
            done
        done
    fi
done

# Count total bins identified
total_bins=$(grep -v "^#" "$OUTPUT_FILE" | wc -l)
echo ""
echo "Summary:"
echo "--------"
echo "Total bins identified: $total_bins"
echo "Output saved to: $OUTPUT_FILE"
echo ""

# Create a summary by treatment
echo "Bins per treatment:"
for treatment_dir in "${OUTPUT_DIR}/binning"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    treatment=$(basename "$treatment_dir")
    count=$(grep -v "^#" "$OUTPUT_FILE" | grep "|${treatment}|" | wc -l)
    echo "  $treatment: $count bins"
done
