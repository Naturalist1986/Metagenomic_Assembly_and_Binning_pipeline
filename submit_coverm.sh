#!/bin/bash
# submit_coverm.sh - Submit CoverM standalone jobs

BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
BIN_REFINEMENT_DIR="${BASE_DIR}/bin_refinement"
COMBINED_READS_DIR="${BASE_DIR}/combined_reads"

echo "========================================="
echo "CoverM Standalone Submission"
echo "========================================="
echo ""
echo "Base directory: $BASE_DIR"
echo ""

# Count samples
sample_count=0
echo "Discovering samples..."
for treatment_dir in "$BIN_REFINEMENT_DIR"/*; do
    if [ ! -d "$treatment_dir" ]; then
        continue
    fi
    
    treatment=$(basename "$treatment_dir")
    
    for sample_dir in "$treatment_dir"/*; do
        if [ ! -d "$sample_dir" ]; then
            continue
        fi
        
        sample=$(basename "$sample_dir")
        bins_dir="${sample_dir}/metawrap_50_10_bins"
        
        if [ -d "$bins_dir" ] && [ "$(ls -A "$bins_dir"/*.fa 2>/dev/null)" ]; then
            read1="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_1.fq.gz"
            read2="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_2.fq.gz"
            
            if [ -f "$read1" ] && [ -f "$read2" ]; then
                bin_count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
                echo "  ✓ ${treatment}/${sample} ($bin_count bins)"
                ((sample_count++))
            else
                echo "  ✗ ${treatment}/${sample} (bins found, reads missing)"
            fi
        fi
    done
done

echo ""
echo "Found $sample_count samples ready for CoverM analysis"
echo ""

if [ $sample_count -eq 0 ]; then
    echo "ERROR: No samples found!"
    echo ""
    echo "Please check:"
    echo "  1. Refined bins exist in: $BIN_REFINEMENT_DIR"
    echo "  2. Combined reads exist in: $COMBINED_READS_DIR"
    exit 1
fi

# Calculate array range
max_index=$((sample_count - 1))

echo "Submitting SLURM array job (0-${max_index})..."
echo ""

# Submit job
job_id=$(sbatch --array=0-${max_index}%10 --parsable 08_coverm_standalone.sh)

if [ $? -eq 0 ]; then
    echo "✓ Job submitted successfully!"
    echo ""
    echo "Job ID: $job_id"
    echo ""
    echo "Monitor progress:"
    echo "  squeue -j $job_id"
    echo ""
    echo "Check logs:"
    echo "  ls ${BASE_DIR}/logs/coverm/"
    echo ""
    echo "View results:"
    echo "  ls ${BASE_DIR}/coverm/*/*/abundance.tsv"
else
    echo "✗ Job submission failed!"
    exit 1
fi
