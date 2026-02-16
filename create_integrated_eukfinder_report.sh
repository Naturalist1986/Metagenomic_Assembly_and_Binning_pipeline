#!/bin/bash
# create_integrated_eukfinder_report.sh
# User-friendly wrapper for bacterial vs non-bacterial mapping summary
#
# This script provides a single-command interface to create an integrated
# summary table showing bacterial, non-bacterial, and unmapped read percentages
# per treatment based on EukFinder classification results.

set -euo pipefail

# Source configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config_utilities.sh"

# ===================================================================
# Functions
# ===================================================================

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

print_banner() {
    echo "==================================================================="
    echo "Creating Integrated EukFinder Mapping Report"
    echo "==================================================================="
    echo ""
    echo "This creates a table with one row per treatment showing:"
    echo "  • Bacterial %: Reads mapping to bacterial sequences (Bact)"
    echo "  • Non-Bacterial %: Reads mapping to Arch, Euk, Unk, EUnk, Misc"
    echo "  • Unmapped %: Reads not mapping to any binned sequences"
    echo ""
}

check_eukfinder_results() {
    # Check if EukFinder results exist and count treatments
    # Returns: number of treatments found (0 = error)

    echo "Step 1: Checking EukFinder results..."
    echo "-------------------------------------------------------------------"

    local eukfinder_dir="${OUTPUT_DIR}/eukfinder_output"

    if [ ! -d "$eukfinder_dir" ]; then
        echo "ERROR: EukFinder directory not found: $eukfinder_dir"
        echo ""
        echo "Please run EukFinder first:"
        echo "  ./submit_eukfinder_all_bins.sh"
        echo ""
        echo "This will classify sequences in all bins into:"
        echo "  - Bacterial (Bact)"
        echo "  - Archaeal (Arch)"
        echo "  - Eukaryotic (Euk)"
        echo "  - Unknown (Unk)"
        echo "  - Eukaryotic+Unknown (EUnk)"
        echo "  - Miscellaneous (Misc)"
        echo ""
        return 0
    fi

    # Count treatments
    local treatment_count=0
    local treatments=()

    for treatment_dir in "${eukfinder_dir}"/*; do
        if [ -d "$treatment_dir" ]; then
            local treatment=$(basename "$treatment_dir")

            # Check if this treatment has EukFinder results
            local has_results=0
            for result_dir in "${treatment_dir}"/*; do
                if [ -d "$result_dir/Eukfinder_results" ] && \
                   [ -f "$result_dir/Eukfinder_results/summary_table.txt" ]; then
                    has_results=1
                    break
                fi
            done

            if [ "$has_results" -eq 1 ]; then
                treatments+=("$treatment")
                ((treatment_count++))
            fi
        fi
    done

    if [ "$treatment_count" -eq 0 ]; then
        echo "ERROR: No treatments with EukFinder results found in $eukfinder_dir"
        echo ""
        echo "Please run EukFinder first:"
        echo "  ./submit_eukfinder_all_bins.sh"
        echo ""
        return 0
    fi

    echo "✓ Found EukFinder results for $treatment_count treatments:"
    for treatment in "${treatments[@]}"; do
        echo "  - $treatment"
    done
    echo ""

    return "$treatment_count"
}

check_mapping_status() {
    # Check if mapping jobs are complete, in progress, or not started
    # Returns: status string ("complete", "partial", or "not_started")

    echo "Step 2: Checking mapping status..."
    echo "-------------------------------------------------------------------"

    local mapping_dir="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping"
    local eukfinder_dir="${OUTPUT_DIR}/eukfinder_output"

    # Check if mapping directory exists
    if [ ! -d "$mapping_dir" ]; then
        echo "Mapping not started yet"
        echo ""
        echo "not_started"
        return 0
    fi

    # Count completed treatments
    local total_treatments=0
    local completed_treatments=0
    local missing_treatments=()

    for treatment_dir in "${eukfinder_dir}"/*; do
        if [ ! -d "$treatment_dir" ]; then
            continue
        fi

        local treatment=$(basename "$treatment_dir")

        # Check if this treatment has EukFinder results
        local has_results=0
        for result_dir in "${treatment_dir}"/*; do
            if [ -d "$result_dir/Eukfinder_results" ]; then
                has_results=1
                break
            fi
        done

        if [ "$has_results" -eq 0 ]; then
            continue
        fi

        ((total_treatments++))

        # Check if mapping summary exists for this treatment
        local summary_file="${mapping_dir}/${treatment}/mapping_summary.txt"
        if [ -f "$summary_file" ]; then
            ((completed_treatments++))
        else
            missing_treatments+=("$treatment")
        fi
    done

    if [ "$completed_treatments" -eq 0 ]; then
        echo "Mapping directory exists but no completed mappings found"
        echo ""
        echo "not_started"
        return 0
    elif [ "$completed_treatments" -eq "$total_treatments" ]; then
        echo "✓ Mapping completed for $completed_treatments/$total_treatments treatments"
        echo ""
        echo "complete"
        return 0
    else
        echo "⚠ Mapping partially completed: $completed_treatments/$total_treatments treatments"
        echo ""
        echo "Missing results for:"
        for treatment in "${missing_treatments[@]}"; do
            echo "  - $treatment"
        done
        echo ""
        echo "partial"
        return 0
    fi
}

check_slurm_jobs() {
    # Check if mapping jobs are currently running in SLURM
    # Returns: 0 if jobs found, 1 if no jobs found

    if ! command -v squeue &> /dev/null; then
        return 1
    fi

    # Check for jobs with name containing "bact_map"
    local job_count=$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -c "bact_map" || true)

    if [ "$job_count" -gt 0 ]; then
        echo "Found $job_count running mapping jobs in queue"
        return 0
    fi

    return 1
}

submit_mapping_jobs() {
    # Submit mapping jobs via existing infrastructure

    echo "Step 3: Submitting mapping jobs..."
    echo "-------------------------------------------------------------------"
    echo ""

    if [ ! -f "${SCRIPT_DIR}/submit_bacterial_vs_nonbacterial_mapping.sh" ]; then
        echo "ERROR: Cannot find submit_bacterial_vs_nonbacterial_mapping.sh"
        echo "Expected location: ${SCRIPT_DIR}/submit_bacterial_vs_nonbacterial_mapping.sh"
        return 1
    fi

    # Run the submission script
    echo "Running submission script..."
    echo ""
    "${SCRIPT_DIR}/submit_bacterial_vs_nonbacterial_mapping.sh"

    echo ""
    echo "==================================================================="
    echo "Mapping jobs submitted!"
    echo "==================================================================="
    echo ""
    echo "Monitor job status:"
    echo "  squeue -u \$USER | grep bact_map"
    echo ""
    echo "Monitor logs:"
    echo "  tail -f \${LOG_DIR}/bacterial_vs_nonbacterial_mapping/bact_map_*.out"
    echo ""
    echo "When jobs complete, run this script again to generate the summary:"
    echo "  ./create_integrated_eukfinder_report.sh"
    echo ""

    return 0
}

generate_summary_table() {
    # Generate final integrated summary using existing script
    # Returns: 0 on success, 1 on failure

    echo "Step 3: Generating integrated summary..."
    echo "-------------------------------------------------------------------"

    if [ ! -f "${SCRIPT_DIR}/summarize_bacterial_vs_nonbacterial.sh" ]; then
        echo "ERROR: Cannot find summarize_bacterial_vs_nonbacterial.sh"
        echo "Expected location: ${SCRIPT_DIR}/summarize_bacterial_vs_nonbacterial.sh"
        return 1
    fi

    # Run the summary script
    if ! bash "${SCRIPT_DIR}/summarize_bacterial_vs_nonbacterial.sh"; then
        echo "ERROR: Failed to generate summary"
        return 1
    fi

    # Verify output files were created
    local tsv_file="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.tsv"
    local txt_file="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.txt"

    if [ ! -f "$tsv_file" ]; then
        echo "ERROR: Summary TSV file not created: $tsv_file"
        return 1
    fi

    echo "✓ Summary table created"
    echo ""

    return 0
}

display_results() {
    # Display results with interpretation guidance

    local tsv_file="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.tsv"
    local txt_file="${OUTPUT_DIR}/bacterial_vs_nonbacterial_mapping/bacterial_vs_nonbacterial_summary.txt"

    echo "==================================================================="
    echo "Results Ready!"
    echo "==================================================================="
    echo ""
    echo "Output files:"
    echo "  TSV:    $tsv_file"
    echo "  Report: $txt_file"
    echo ""
    echo "Table columns:"
    echo "  Treatment           - Treatment name"
    echo "  Total_Reads         - Total read pairs analyzed"
    echo "  Bacterial_Reads     - Reads mapping to bacteria"
    echo "  Bacterial_%         - Percent mapping to bacteria"
    echo "  Non-Bacterial_Reads - Reads mapping to non-bacteria"
    echo "  Non-Bacterial_%     - Percent mapping to non-bacteria"
    echo "  Unmapped_Reads      - Reads unmapped"
    echo "  Unmapped_%          - Percent unmapped"
    echo ""
    echo "View commands:"
    echo "  # Formatted text report"
    echo "  cat $txt_file"
    echo ""
    echo "  # TSV in columns"
    echo "  column -t -s \$'\\t' $tsv_file | less -S"
    echo ""
    echo "Import into Excel/R:"
    echo "  Open: $tsv_file"
    echo ""

    # Display sample of results if file exists
    if [ -f "$txt_file" ]; then
        echo "Preview of results:"
        echo "-------------------------------------------------------------------"
        head -20 "$txt_file" 2>/dev/null || true
        echo ""
    fi
}

display_partial_status() {
    # Display status for partial completion

    echo "==================================================================="
    echo "Mapping In Progress or Incomplete"
    echo "==================================================================="
    echo ""

    # Check if jobs are running
    if check_slurm_jobs; then
        echo "Status: Jobs are currently running"
        echo ""
        echo "Check job status:"
        echo "  squeue -u \$USER | grep bact_map"
        echo ""
        echo "Monitor logs:"
        echo "  tail -f \${LOG_DIR}/bacterial_vs_nonbacterial_mapping/bact_map_*.out"
        echo ""
        echo "Run this script again when jobs complete to generate the summary."
        echo ""
    else
        echo "Status: No jobs running - some mappings may have failed"
        echo ""
        echo "Check error logs:"
        echo "  less \${LOG_DIR}/bacterial_vs_nonbacterial_mapping/bact_map_*.err"
        echo ""
        echo "To rerun failed jobs:"
        echo "  ./submit_bacterial_vs_nonbacterial_mapping.sh"
        echo ""
        echo "Or generate summary with available data:"
        echo "  ./summarize_bacterial_vs_nonbacterial.sh"
        echo ""
    fi
}

# ===================================================================
# Main Workflow
# ===================================================================

main() {
    print_banner

    # Step 1: Check EukFinder results exist
    local treatment_count
    treatment_count=$(check_eukfinder_results)

    if [ "$treatment_count" -eq 0 ]; then
        exit 1
    fi

    # Step 2: Check mapping status
    local status_output
    status_output=$(check_mapping_status)
    local mapping_status=$(echo "$status_output" | tail -1)

    # Step 3: Take action based on status
    case "$mapping_status" in
        "complete")
            # Generate summary table
            if generate_summary_table; then
                display_results
                exit 0
            else
                echo "ERROR: Failed to generate summary table"
                exit 1
            fi
            ;;

        "not_started")
            # Submit mapping jobs
            if submit_mapping_jobs; then
                exit 0
            else
                echo "ERROR: Failed to submit mapping jobs"
                exit 1
            fi
            ;;

        "partial")
            # Show status and guidance
            display_partial_status
            exit 0
            ;;

        *)
            echo "ERROR: Unknown mapping status: $mapping_status"
            exit 1
            ;;
    esac
}

# Run main workflow
main
