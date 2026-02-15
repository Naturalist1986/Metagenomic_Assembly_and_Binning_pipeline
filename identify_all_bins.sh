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

# Function to process COMEBin bins
process_comebin_bins() {
    local treatment_dir="$1"
    local treatment="$2"
    local sample="$3"
    local comebin_dir="${treatment_dir}/comebin"

    if [ ! -d "$comebin_dir" ]; then
        return 0
    fi

    # COMEBin has nested structure
    local bin_dir=""
    if [ -d "${comebin_dir}/comebin_res/comebin_res_bins" ]; then
        bin_dir="${comebin_dir}/comebin_res/comebin_res_bins"
    elif [ -d "${comebin_dir}/comebin_res_bins" ]; then
        bin_dir="${comebin_dir}/comebin_res_bins"
    fi

    if [ -z "$bin_dir" ] || [ ! -d "$bin_dir" ]; then
        return 0
    fi

    local bin_count=0
    for bin in "$bin_dir"/*.fa; do
        if [ ! -f "$bin" ]; then
            continue
        fi

        bin_name=$(basename "$bin")
        bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
        num_contigs=$(grep -c "^>" "$bin")

        echo "${bin}|${treatment}|${sample}|comebin|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
        ((bin_count++))
    done

    if [ $bin_count -gt 0 ]; then
        echo "      comebin: $bin_count bins"
    fi
}

# Function to process SemiBin bins
process_semibin_bins() {
    local treatment_dir="$1"
    local treatment="$2"
    local sample="$3"
    local semibin_dir="${treatment_dir}/semibin"

    if [ ! -d "$semibin_dir" ]; then
        return 0
    fi

    # SemiBin outputs to output_bins subdirectory
    local bin_dir="${semibin_dir}/output_bins"
    if [ ! -d "$bin_dir" ]; then
        return 0
    fi

    local bin_count=0
    for bin in "$bin_dir"/*.fa; do
        if [ ! -f "$bin" ]; then
            continue
        fi

        bin_name=$(basename "$bin")
        bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
        num_contigs=$(grep -c "^>" "$bin")

        echo "${bin}|${treatment}|${sample}|semibin|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
        ((bin_count++))
    done

    if [ $bin_count -gt 0 ]; then
        echo "      semibin: $bin_count bins"
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

        # Process COMEBin and SemiBin
        process_comebin_bins "$treatment_dir" "$treatment" "$treatment"
        process_semibin_bins "$treatment_dir" "$treatment" "$treatment"
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

            # Process COMEBin and SemiBin
            process_comebin_bins "$sample_dir" "$treatment" "$sample"
            process_semibin_bins "$sample_dir" "$treatment" "$sample"
        done
    fi

    # Process bin refinement outputs (DAS Tool and Binette)
    refinement_dir="${OUTPUT_DIR}/bin_refinement/${treatment}"

    if [ -d "$refinement_dir" ]; then
        echo "  Checking bin refinement results for treatment: $treatment"

        # Check for treatment-level refinement (coassembly)
        if [ -d "${refinement_dir}/dastool/dastool_DASTool_bins" ]; then
            echo "    Processing treatment-level DAS Tool bins"
            bin_count=0
            for bin in "${refinement_dir}/dastool/dastool_DASTool_bins"/*.fa; do
                if [ ! -f "$bin" ]; then
                    continue
                fi

                bin_name=$(basename "$bin")
                bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                num_contigs=$(grep -c "^>" "$bin")

                echo "${bin}|${treatment}|${treatment}|dastool|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                ((bin_count++))
            done
            if [ $bin_count -gt 0 ]; then
                echo "      dastool: $bin_count bins"
            fi
        fi

        if [ -d "${refinement_dir}/binette/final_bins" ]; then
            echo "    Processing treatment-level Binette bins"
            bin_count=0
            for bin in "${refinement_dir}/binette/final_bins"/*.fa; do
                if [ ! -f "$bin" ]; then
                    continue
                fi

                bin_name=$(basename "$bin")
                bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                num_contigs=$(grep -c "^>" "$bin")

                echo "${bin}|${treatment}|${treatment}|binette|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                ((bin_count++))
            done
            if [ $bin_count -gt 0 ]; then
                echo "      binette: $bin_count bins"
            fi
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
                bin_count=0
                for bin in "${sample_ref_dir}/dastool/dastool_DASTool_bins"/*.fa; do
                    if [ ! -f "$bin" ]; then
                        continue
                    fi

                    bin_name=$(basename "$bin")
                    bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    num_contigs=$(grep -c "^>" "$bin")

                    echo "${bin}|${treatment}|${sample}|dastool|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                    ((bin_count++))
                done
                if [ $bin_count -gt 0 ]; then
                    echo "        dastool: $bin_count bins"
                fi
            fi

            if [ -d "${sample_ref_dir}/binette/final_bins" ]; then
                bin_count=0
                for bin in "${sample_ref_dir}/binette/final_bins"/*.fa; do
                    if [ ! -f "$bin" ]; then
                        continue
                    fi

                    bin_name=$(basename "$bin")
                    bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    num_contigs=$(grep -c "^>" "$bin")

                    echo "${bin}|${treatment}|${sample}|binette|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                    ((bin_count++))
                done
                if [ $bin_count -gt 0 ]; then
                    echo "        binette: $bin_count bins"
                fi
            fi
        done
    fi

    # Process selected bins (from 08_bin_selection.sh)
    selected_bins_dir="${OUTPUT_DIR}/selected_bins/${treatment}"

    if [ -d "$selected_bins_dir" ]; then
        echo "  Checking selected bins for treatment: $treatment"

        # Check for treatment-level selected bins (coassembly)
        if ls "${selected_bins_dir}"/*.fa &>/dev/null; then
            echo "    Processing treatment-level selected bins"
            bin_count=0
            for bin in "${selected_bins_dir}"/*.fa; do
                if [ ! -f "$bin" ]; then
                    continue
                fi

                bin_name=$(basename "$bin")
                bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                num_contigs=$(grep -c "^>" "$bin")

                echo "${bin}|${treatment}|${treatment}|selected|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                ((bin_count++))
            done
            if [ $bin_count -gt 0 ]; then
                echo "      selected: $bin_count bins"
            fi
        fi

        # Check for sample-level selected bins (individual assembly)
        for sample_sel_dir in "$selected_bins_dir"/*; do
            if [ ! -d "$sample_sel_dir" ]; then
                continue
            fi

            sample=$(basename "$sample_sel_dir")
            echo "    Checking selected bins for sample: $sample"

            if ls "${sample_sel_dir}"/*.fa &>/dev/null; then
                bin_count=0
                for bin in "${sample_sel_dir}"/*.fa; do
                    if [ ! -f "$bin" ]; then
                        continue
                    fi

                    bin_name=$(basename "$bin")
                    bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
                    num_contigs=$(grep -c "^>" "$bin")

                    echo "${bin}|${treatment}|${sample}|selected|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
                    ((bin_count++))
                done
                if [ $bin_count -gt 0 ]; then
                    echo "        selected: $bin_count bins"
                fi
            fi
        done
    fi

    # Process final bin collection (from 09_bin_collection.sh)
    collection_dir="${OUTPUT_DIR}/bin_collection/${treatment}/bins"

    if [ -d "$collection_dir" ] && ls "${collection_dir}"/*.fa &>/dev/null; then
        echo "  Processing final bin collection for treatment: $treatment"
        bin_count=0
        for bin in "${collection_dir}"/*.fa; do
            if [ ! -f "$bin" ]; then
                continue
            fi

            bin_name=$(basename "$bin")
            bin_size=$(grep -v "^>" "$bin" | tr -d '\n' | wc -c)
            num_contigs=$(grep -c "^>" "$bin")

            echo "${bin}|${treatment}|${treatment}|collection|${bin_name}|${bin_size}|${num_contigs}" >> "$OUTPUT_FILE"
            ((bin_count++))
        done
        if [ $bin_count -gt 0 ]; then
            echo "    collection: $bin_count bins"
        fi
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
for binner in metabat2 maxbin2 concoct comebin semibin dastool binette selected collection; do
    count=$(grep -v "^#" "$OUTPUT_FILE" | grep "|${binner}|" | wc -l)
    if [ $count -gt 0 ]; then
        echo "  $binner: $count bins"
    fi
done
