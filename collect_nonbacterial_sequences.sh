#!/bin/bash
# collect_nonbacterial_sequences.sh - Collect non-bacterial sequences from EukFinder results
# Includes: Archaeal, Eukaryotic, Unknown, Euk+Unknown, and Misc

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Collecting Non-Bacterial Sequences from EukFinder Results"
echo "====================================================================="
echo ""

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder_output"
NONBACTERIAL_DIR="${OUTPUT_DIR}/nonbacterial_sequences"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder results directory not found: $EUKFINDER_DIR"
    exit 1
fi

mkdir -p "$NONBACTERIAL_DIR"

echo "Processing EukFinder results from: $EUKFINDER_DIR"
echo "Output directory: $NONBACTERIAL_DIR"
echo ""

# Categories to collect (everything except Bact)
NONBACT_CATEGORIES="Arch Euk Unk EUnk Misc"

# Process each treatment
for treatment_dir in "$EUKFINDER_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"

    # Create output file for this treatment
    output_fasta="${NONBACTERIAL_DIR}/${treatment}_nonbacterial.fasta"
    > "$output_fasta"  # Clear/create file

    declare -A category_counts
    for cat in $NONBACT_CATEGORIES; do
        category_counts[$cat]=0
    done

    total_seqs=0
    total_bp=0

    # Find all bin directories for this treatment
    for bin_dir in "$treatment_dir"/*; do
        if [ ! -d "$bin_dir" ]; then
            continue
        fi

        results_dir="${bin_dir}/Eukfinder_results"
        if [ ! -d "$results_dir" ]; then
            continue
        fi

        bin_name=$(basename "$bin_dir")
        echo "  Processing bin: $bin_name"

        # Collect sequences from each non-bacterial category
        for category in $NONBACT_CATEGORIES; do
            cat_file="${results_dir}/${category}.fasta"

            if [ -f "$cat_file" ] && [ -s "$cat_file" ]; then
                # Count sequences before adding
                seq_count=$(grep -c "^>" "$cat_file" 2>/dev/null || echo "0")

                if [ $seq_count -gt 0 ]; then
                    # Add bin name and category to headers for tracking
                    awk -v bin="$bin_name" -v cat="$category" '/^>/ {print $0"|bin:"bin"|cat:"cat; next} {print}' "$cat_file" >> "$output_fasta"

                    category_counts[$category]=$((${category_counts[$category]} + seq_count))
                    total_seqs=$((total_seqs + seq_count))
                    echo "    $category: $seq_count sequences"
                fi
            fi
        done
    done

    # Calculate total size
    if [ -s "$output_fasta" ]; then
        total_bp=$(grep -v "^>" "$output_fasta" | tr -d '\n' | wc -c)
        echo ""
        echo "  Treatment $treatment summary:"
        for category in $NONBACT_CATEGORIES; do
            if [ ${category_counts[$category]} -gt 0 ]; then
                echo "    $category: ${category_counts[$category]} sequences"
            fi
        done
        echo "    Total non-bacterial sequences: $total_seqs"
        echo "    Total size: $total_bp bp"
        echo "    Output: $output_fasta"
        echo ""
    else
        echo "  No non-bacterial sequences found for treatment $treatment"
        rm -f "$output_fasta"
    fi
done

# Summary
echo "====================================================================="
echo "Non-Bacterial Sequence Collection Complete"
echo "====================================================================="
echo ""

num_treatments=$(ls -1 "${NONBACTERIAL_DIR}"/*_nonbacterial.fasta 2>/dev/null | wc -l)
echo "Processed $num_treatments treatments"
echo "Output directory: $NONBACTERIAL_DIR"
echo ""

if [ $num_treatments -eq 0 ]; then
    echo "WARNING: No non-bacterial sequences found"
    echo "Make sure EukFinder has completed successfully"
fi
