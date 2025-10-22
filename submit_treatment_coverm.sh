#!/bin/bash
# submit_treatment_coverm.sh - Submit treatment-level CoverM jobs

BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
BIN_REFINEMENT_DIR="${BASE_DIR}/bin_refinement"
COMBINED_READS_DIR="${BASE_DIR}/combined_reads"
TREATMENT_LIST="${BASE_DIR}/treatment_list.txt"

echo "========================================="
echo "Treatment-Level CoverM Submission"
echo "========================================="
echo ""
echo "Approach: Map ALL bins in each treatment"
echo "          to each sample's reads within that treatment"
echo ""
echo "Base directory: $BASE_DIR"
echo ""

# Create treatment list
> "$TREATMENT_LIST"

echo "Discovering treatments..."
for treatment_dir in "$BIN_REFINEMENT_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    
    treatment=$(basename "$treatment_dir")
    
    # Count bins and samples for this treatment
    total_bins=0
    sample_count=0
    sample_list=""
    
    for sample_dir in "$treatment_dir"/*; do
        if [ ! -d "$sample_dir" ]; then
            continue
        fi
        
        sample=$(basename "$sample_dir")
        bins_dir="${sample_dir}/metawrap_50_10_bins"
        
        if [ -d "$bins_dir" ]; then
            bin_count=$(find "$bins_dir" -name "*.fa" -type f -exec sh -c 'grep -v "^>" "$1" | tr -d "\n" | wc -c' _ {} \; | awk '$1 >= 10000' | wc -l)
            total_bins=$((total_bins + bin_count))
        fi
        
        read1="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_1.fq.gz"
        read2="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_2.fq.gz"
        
        if [ -f "$read1" ] && [ -f "$read2" ]; then
            ((sample_count++))
            sample_list="${sample_list}${sample}, "
        fi
    done
    
    if [ $total_bins -gt 0 ] && [ $sample_count -gt 0 ]; then
        echo "$treatment" >> "$TREATMENT_LIST"
        sample_list="${sample_list%, }"  # Remove trailing comma
        echo "  ✓ $treatment: $total_bins bins across $sample_count samples"
        echo "    Samples: $sample_list"
    else
        echo "  ✗ $treatment: No valid bins or samples"
    fi
done

echo ""
treatment_count=$(wc -l < "$TREATMENT_LIST")
echo "Found $treatment_count treatments ready for analysis"
echo "Treatment list saved to: $TREATMENT_LIST"
echo ""

if [ $treatment_count -eq 0 ]; then
    echo "ERROR: No treatments found!"
    echo ""
    echo "Please check:"
    echo "  1. Refined bins exist in: $BIN_REFINEMENT_DIR"
    echo "  2. Combined reads exist in: $COMBINED_READS_DIR"
    exit 1
fi

# Show what will happen
echo "========================================="
echo "ANALYSIS PLAN"
echo "========================================="
cat "$TREATMENT_LIST" | while read treatment; do
    echo ""
    echo "Treatment: $treatment"
    echo "  1. Collect all bins from all ${treatment} samples"
    echo "  2. Map these bins to EACH sample's reads"
    echo "  3. Generate per-sample abundance profiles"
    echo "  4. Create treatment-wide summary"
done
echo ""
echo "========================================="
echo ""

# Calculate array range
max_index=$((treatment_count - 1))

# Estimate resources
echo "Resource Requirements:"
echo "  - Treatments: $treatment_count"
echo "  - CPUs per job: 32"
echo "  - Memory per job: 128 GB"
echo "  - Time per job: 24 hours"
echo ""

read -p "Proceed with submission? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Submission cancelled"
    exit 0
fi

echo "Submitting SLURM array job (0-${max_index})..."
echo ""

# Submit job
job_id=$(sbatch --array=0-${max_index}%6 --parsable 08_coverm_treatment_level.sh)

if [ $? -eq 0 ]; then
    echo "✓ Job submitted successfully!"
    echo ""
    echo "Job ID: $job_id"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -j $job_id"
    echo "  watch -n 10 'squeue -j $job_id'"
    echo ""
    echo "Check logs (as they run):"
    echo "  tail -f ${BASE_DIR}/logs/coverm_treatment/coverm_${job_id}_*.out"
    echo ""
    echo "Check for errors:"
    echo "  tail ${BASE_DIR}/logs/coverm_treatment/coverm_${job_id}_*.err"
    echo ""
    echo "View results (after completion):"
    echo "  # Overall treatment summaries"
    echo "  ls ${BASE_DIR}/coverm_treatment_level/*/TREATMENT_SUMMARY.txt"
    echo ""
    echo "  # Example: View RH treatment summary"
    echo "  cat ${BASE_DIR}/coverm_treatment_level/RH/TREATMENT_SUMMARY.txt"
    echo ""
    echo "  # Example: View RH_8 sample abundance"
    echo "  less -S ${BASE_DIR}/coverm_treatment_level/RH/RH_8_abundance.tsv"
    echo ""
    echo "  # Example: View RH_8 sample summary"
    echo "  cat ${BASE_DIR}/coverm_treatment_level/RH/RH_8_summary.txt"
else
    echo "✗ Job submission failed!"
    exit 1
fi
