#!/bin/bash
# collect_bacterial_sequences.sh - Collect bacterial sequences from EukFinder results

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Collecting Bacterial Sequences from EukFinder Results"
echo "====================================================================="
echo ""

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder_output"
BACTERIAL_DIR="${OUTPUT_DIR}/bacterial_sequences"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder results directory not found: $EUKFINDER_DIR"
    exit 1
fi

mkdir -p "$BACTERIAL_DIR"

echo "Processing EukFinder results from: $EUKFINDER_DIR"
echo "Output directory: $BACTERIAL_DIR"
echo ""

# Process each treatment
for treatment_dir in "$EUKFINDER_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"

    # Create output file for this treatment
    output_fasta="${BACTERIAL_DIR}/${treatment}_bacterial.fasta"
    > "$output_fasta"  # Clear/create file

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

        # Collect bacterial sequences (Bact only)
        bact_file="${results_dir}/Bact.fasta"

        if [ -f "$bact_file" ] && [ -s "$bact_file" ]; then
            # Count sequences before adding
            seq_count=$(grep -c "^>" "$bact_file" 2>/dev/null || echo "0")

            if [ $seq_count -gt 0 ]; then
                # Add bin name to headers for tracking
                awk -v bin="$bin_name" '/^>/ {print $0"|bin:"bin; next} {print}' "$bact_file" >> "$output_fasta"

                total_seqs=$((total_seqs + seq_count))
                echo "    Bact: $seq_count sequences"
            fi
        fi
    done

    # Calculate total size
    if [ -s "$output_fasta" ]; then
        total_bp=$(grep -v "^>" "$output_fasta" | tr -d '\n' | wc -c)
        echo ""
        echo "  Treatment $treatment summary:"
        echo "    Total bacterial sequences: $total_seqs"
        echo "    Total size: $total_bp bp"
        echo "    Output: $output_fasta"
        echo ""
    else
        echo "  No bacterial sequences found for treatment $treatment"
        rm -f "$output_fasta"
    fi
done

# Summary
echo "====================================================================="
echo "Bacterial Sequence Collection Complete"
echo "====================================================================="
echo ""

num_treatments=$(ls -1 "${BACTERIAL_DIR}"/*_bacterial.fasta 2>/dev/null | wc -l)
echo "Processed $num_treatments treatments"
echo "Output directory: $BACTERIAL_DIR"
echo ""

if [ $num_treatments -eq 0 ]; then
    echo "WARNING: No bacterial sequences found"
    echo "Make sure EukFinder has completed successfully"
fi
