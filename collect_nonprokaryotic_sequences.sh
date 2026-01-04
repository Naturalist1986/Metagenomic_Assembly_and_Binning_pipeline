#!/bin/bash
# collect_nonprokaryotic_sequences.sh - Collect non-prokaryotic sequences per treatment

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Collecting Non-Prokaryotic Sequences from EukFinder Results"
echo "====================================================================="
echo ""

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder"
NONPROK_DIR="${OUTPUT_DIR}/nonprokaryotic_sequences"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder results directory not found: $EUKFINDER_DIR"
    exit 1
fi

mkdir -p "$NONPROK_DIR"

echo "Processing EukFinder results from: $EUKFINDER_DIR"
echo "Output directory: $NONPROK_DIR"
echo ""

# Categories to collect (everything except Bact and Arch)
NONPROK_CATEGORIES="Euk Unk EUnk Misc"

# Process each treatment
for treatment_dir in "$EUKFINDER_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")
    echo "Processing treatment: $treatment"

    # Create output file for this treatment
    output_fasta="${NONPROK_DIR}/${treatment}_nonprokaryotic.fasta"
    > "$output_fasta"  # Clear/create file

    total_seqs=0
    total_bp=0

    # Find all bin directories for this treatment
    for bin_dir in "$treatment_dir"/*; do
        if [ ! -d "$bin_dir" ]; then
            continue
        fi

        bin_name=$(basename "$bin_dir")
        results_dir="${bin_dir}/Eukfinder_results"

        if [ ! -d "$results_dir" ]; then
            continue
        fi

        echo "  Processing bin: $bin_name"

        # Collect sequences from non-prokaryotic categories
        for category in $NONPROK_CATEGORIES; do
            fasta_file="${results_dir}/${bin_name}.${category}.fasta"

            if [ ! -f "$fasta_file" ]; then
                continue
            fi

            # Check if file is not empty
            if [ ! -s "$fasta_file" ]; then
                continue
            fi

            # Count sequences in this file
            num_seqs=$(grep -c "^>" "$fasta_file" 2>/dev/null || echo "0")

            if [ "$num_seqs" -gt 0 ]; then
                echo "    Found $num_seqs sequences in ${category}.fasta"

                # Append to combined file with modified headers to include source info
                awk -v bin="$bin_name" -v cat="$category" '
                    /^>/ {
                        # Modify header to include bin and category info
                        print $0 "|bin=" bin "|category=" cat
                        next
                    }
                    {print}
                ' "$fasta_file" >> "$output_fasta"

                total_seqs=$((total_seqs + num_seqs))
            fi
        done
    done

    # Calculate total size
    if [ -s "$output_fasta" ]; then
        total_bp=$(grep -v "^>" "$output_fasta" | tr -d '\n' | wc -c)
        echo "  Total non-prokaryotic sequences for $treatment: $total_seqs sequences, $total_bp bp"
        echo ""
    else
        echo "  No non-prokaryotic sequences found for $treatment"
        echo ""
        rm -f "$output_fasta"
    fi
done

# Create summary file
summary_file="${NONPROK_DIR}/collection_summary.txt"
cat > "$summary_file" << EOF
Non-Prokaryotic Sequence Collection Summary
============================================
Date: $(date)

Sequences collected from categories: $NONPROK_CATEGORIES
(Excludes: Bact, Arch)

Per Treatment Summary:
EOF

for fasta in "${NONPROK_DIR}"/*_nonprokaryotic.fasta; do
    if [ -f "$fasta" ]; then
        treatment=$(basename "$fasta" _nonprokaryotic.fasta)
        num_seqs=$(grep -c "^>" "$fasta" 2>/dev/null || echo "0")
        size_bp=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c 2>/dev/null || echo "0")

        printf "  %-15s %10d sequences, %15d bp\n" "$treatment:" "$num_seqs" "$size_bp" >> "$summary_file"
    fi
done

echo ""
echo "====================================================================="
echo "Collection Complete"
echo "====================================================================="
echo ""
cat "$summary_file"
echo ""
echo "Combined FASTA files saved to: $NONPROK_DIR"
echo "Ready for read mapping with: ./map_to_nonprokaryotic.sh"
echo ""
