#!/bin/bash
# regenerate_eukfinder_summaries.sh - Regenerate EukFinder summary files from existing results

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

echo "====================================================================="
echo "Regenerating EukFinder Summary Files"
echo "====================================================================="
echo ""

EUKFINDER_DIR="${OUTPUT_DIR}/eukfinder"

if [ ! -d "$EUKFINDER_DIR" ]; then
    echo "ERROR: EukFinder results directory not found: $EUKFINDER_DIR"
    exit 1
fi

echo "Searching for EukFinder results in: $EUKFINDER_DIR"
echo ""

# Find all directories with Eukfinder_results subdirectories
RESULTS_COUNT=0

for treatment_dir in "$EUKFINDER_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi

    treatment=$(basename "$treatment_dir")

    for sample_dir in "$treatment_dir"/*; do
        if [ ! -d "$sample_dir" ]; then
            continue
        fi

        sample=$(basename "$sample_dir")
        results_dir="${sample_dir}/Eukfinder_results"

        if [ ! -d "$results_dir" ]; then
            continue
        fi

        summary_table="${results_dir}/summary_table.txt"
        summary_file="${sample_dir}/eukfinder_summary.txt"

        if [ ! -f "$summary_table" ]; then
            echo "  Skipping $treatment/$sample - no summary_table.txt found"
            continue
        fi

        echo "Processing: $treatment/$sample"

        # Create new summary file
        cat > "$summary_file" << EOF
# EukFinder Results Summary
# Auto-generated from summary_table.txt

Treatment: $treatment
Sample: $sample
Date regenerated: $(date)
EOF
        echo "" >> "$summary_file"
        echo "Results Summary:" >> "$summary_file"
        echo "----------------" >> "$summary_file"

        # Parse the summary_table.txt
        while IFS=$'\t' read -r group num_seq total_size; do
            # Skip the header line
            if [[ "$group" == "Group" ]]; then
                continue
            fi
            # Format and add to summary
            printf "  %-6s %8s sequences, %12s bp\n" "$group:" "$num_seq" "$total_size" >> "$summary_file"
        done < "$summary_table"

        echo "" >> "$summary_file"
        echo "Summary updated: $summary_file"

        ((RESULTS_COUNT++))
    done
done

echo ""
echo "====================================================================="
echo "Summary Regeneration Complete"
echo "====================================================================="
echo ""
echo "Processed $RESULTS_COUNT result sets"
echo ""

if [ $RESULTS_COUNT -eq 0 ]; then
    echo "No EukFinder results found to process"
    echo "Expected structure:"
    echo "  ${EUKFINDER_DIR}/<treatment>/<sample>/Eukfinder_results/summary_table.txt"
fi
