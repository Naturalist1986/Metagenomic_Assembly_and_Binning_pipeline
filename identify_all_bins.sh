#!/bin/bash
# identify_all_bins.sh - Identify ALL bins from each binning tool for EukFinder analysis

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# Parameters
OUTPUT_FILE="${1:-${OUTPUT_DIR}/all_bins_list.txt}"

echo "Identifying ALL bins from each binning tool for EukFinder analysis..."

# Create output file header
cat > "$OUTPUT_FILE" << 'EOF'
# All bins identified for EukFinder analysis
# Format: BIN_PATH|TREATMENT|SAMPLE|BINNER|BIN_NAME|SIZE_BP|NUM_CONTIGS
EOF

# Scan all binning directories
for treatment_dir in "${OUTPUT_DIR}/binning"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"

    # Check if this is a coassembly structure (bins directly under treatment)
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

            binner_name=$(echo "$binner" | sed 's/_bins$//')
            bin_count=0

            # Process ALL bins in this binner directory
            for bin in "$binner_dir"/*.fa; do
                if [ ! -f "$bin" ]; then
                    continue
                fi

                bin_name=$(basename "$bin")
                # Calculate bin size (total bp)
                bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                # Count contigs
                num_contigs=$(grep -c "^>" "$bin")

                # Add to output file
                echo "${bin}|${treatment}|${treatment}|${binner_name}|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                ((bin_count++))
            done

            if [ $bin_count -gt 0 ]; then
                echo "    ${binner_name}: $bin_count bins"
            fi
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

                binner_name=$(echo "$binner" | sed 's/_bins$//')
                bin_count=0

                # Process ALL bins in this binner directory
                for bin in "$binner_dir"/*.fa; do
                    if [ ! -f "$bin" ]; then
                        continue
                    fi

                    bin_name=$(basename "$bin")
                    # Calculate bin size (total bp)
                    bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    # Count contigs
                    num_contigs=$(grep -c "^>" "$bin")

                    # Add to output file
                    echo "${bin}|${treatment}|${sample}|${binner_name}|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                    ((bin_count++))
                done

                if [ $bin_count -gt 0 ]; then
                    echo "      ${binner_name}: $bin_count bins"
                fi
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

echo ""
echo "Summary by binner:"
for binner in metabat2 maxbin2 concoct; do
    count=$(grep -v "^#" "$OUTPUT_FILE" | grep "|${binner}|" | wc -l)
    echo "  $binner: $count bins"
done
