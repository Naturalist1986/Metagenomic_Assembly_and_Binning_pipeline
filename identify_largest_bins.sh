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

# Function to process COMEBin output (nested directory structure)
process_comebin_dir() {
    local treatment_dir="$1"
    local treatment="$2"
    local sample="$3"

    local comebin_dir="${treatment_dir}/comebin"
    if [ ! -d "$comebin_dir" ]; then
        return
    fi

    # COMEBin has nested structure: comebin/comebin_res/comebin_res_bins/
    local bin_dir=""
    if [ -d "${comebin_dir}/comebin_res/comebin_res_bins" ]; then
        bin_dir="${comebin_dir}/comebin_res/comebin_res_bins"
    elif [ -d "${comebin_dir}/comebin_res_bins" ]; then
        bin_dir="${comebin_dir}/comebin_res_bins"
    fi

    if [ -n "$bin_dir" ] && [ -d "$bin_dir" ]; then
        process_binner_dir "$bin_dir" "$treatment" "$sample" "comebin"
    fi
}

# Function to process SemiBin output (output_bins subdirectory)
process_semibin_dir() {
    local treatment_dir="$1"
    local treatment="$2"
    local sample="$3"

    local semibin_dir="${treatment_dir}/semibin"
    if [ ! -d "$semibin_dir" ]; then
        return
    fi

    # SemiBin outputs to output_bins subdirectory
    local bin_dir="${semibin_dir}/output_bins"
    if [ -d "$bin_dir" ]; then
        process_binner_dir "$bin_dir" "$treatment" "$sample" "semibin"
    fi
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

        # Process standard binners
        for binner in metabat2_bins maxbin2_bins concoct_bins; do
            binner_dir="${treatment_dir}/${binner}"

            if [ ! -d "$binner_dir" ]; then
                continue
            fi

            # For coassembly, use treatment name as sample
            process_binner_dir "$binner_dir" "$treatment" "$treatment" "$binner"
        done

        # Process COMEBin (nested structure)
        process_comebin_dir "$treatment_dir" "$treatment" "$treatment"

        # Process SemiBin (output_bins subdirectory)
        process_semibin_dir "$treatment_dir" "$treatment" "$treatment"
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

            # Check standard binner directories
            for binner in metabat2_bins maxbin2_bins concoct_bins; do
                binner_dir="${sample_dir}/${binner}"

                if [ ! -d "$binner_dir" ]; then
                    continue
                fi

                process_binner_dir "$binner_dir" "$treatment" "$sample" "$binner"
            done

            # Process COMEBin (nested structure)
            process_comebin_dir "$sample_dir" "$treatment" "$sample"

            # Process SemiBin (output_bins subdirectory)
            process_semibin_dir "$sample_dir" "$treatment" "$sample"
        done
    fi

    # Process bin refinement outputs (DAS Tool and Binette)
    # These exist at the bin_refinement level
    refinement_dir="${OUTPUT_DIR}/bin_refinement/${treatment}"

    if [ -d "$refinement_dir" ]; then
        echo "  Checking bin refinement results for treatment: $treatment"

        # Check for treatment-level refinement (coassembly)
        if [ -d "${refinement_dir}/dastool/dastool_DASTool_bins" ]; then
            echo "    Found treatment-level DAS Tool bins"
            process_binner_dir "${refinement_dir}/dastool/dastool_DASTool_bins" "$treatment" "$treatment" "dastool"
        fi

        if [ -d "${refinement_dir}/binette/final_bins" ]; then
            echo "    Found treatment-level Binette bins"
            process_binner_dir "${refinement_dir}/binette/final_bins" "$treatment" "$treatment" "binette"
        fi

        # Check for sample-level refinement (individual assembly)
        for sample_ref_dir in "$refinement_dir"/*; do
            if [ ! -d "$sample_ref_dir" ]; then
                continue
            fi

            sample=$(basename "$sample_ref_dir")

            # Skip non-sample directories
            if [[ "$sample" == "dastool" ]] || [[ "$sample" == "binette" ]]; then
                continue
            fi

            echo "    Checking refinement for sample: $sample"

            if [ -d "${sample_ref_dir}/dastool/dastool_DASTool_bins" ]; then
                process_binner_dir "${sample_ref_dir}/dastool/dastool_DASTool_bins" "$treatment" "$sample" "dastool"
            fi

            if [ -d "${sample_ref_dir}/binette/final_bins" ]; then
                process_binner_dir "${sample_ref_dir}/binette/final_bins" "$treatment" "$sample" "binette"
            fi
        done
    fi

    # Process selected bins (from 08_bin_selection.sh)
    selected_bins_dir="${OUTPUT_DIR}/selected_bins/${treatment}"

    if [ -d "$selected_bins_dir" ]; then
        echo "  Checking selected bins for treatment: $treatment"

        # Check for treatment-level selected bins (coassembly)
        if ls "${selected_bins_dir}"/*.fa &>/dev/null; then
            echo "    Found treatment-level selected bins"
            process_binner_dir "$selected_bins_dir" "$treatment" "$treatment" "selected"
        fi

        # Check for sample-level selected bins (individual assembly)
        for sample_sel_dir in "$selected_bins_dir"/*; do
            if [ ! -d "$sample_sel_dir" ]; then
                continue
            fi

            sample=$(basename "$sample_sel_dir")
            echo "    Checking selected bins for sample: $sample"

            if ls "${sample_sel_dir}"/*.fa &>/dev/null; then
                process_binner_dir "$sample_sel_dir" "$treatment" "$sample" "selected"
            fi
        done
    fi

    # Process final bin collection (from 09_bin_collection.sh)
    collection_dir="${OUTPUT_DIR}/bin_collection/${treatment}/bins"

    if [ -d "$collection_dir" ] && ls "${collection_dir}"/*.fa &>/dev/null; then
        echo "  Found final bin collection for treatment: $treatment"
        process_binner_dir "$collection_dir" "$treatment" "$treatment" "collection"
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
