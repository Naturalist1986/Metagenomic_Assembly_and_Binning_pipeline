#!/bin/bash

# merge_lanes.sh - Merge multiple sequencing lanes with validation and repair

set -e

INPUT_DIR="${1:?Error: Please provide input directory}"
OUTPUT_DIR="${2:?Error: Please provide output directory}"
THREADS="${3:-8}"  # Default 8 threads

# Conda environment with BBMap
BBMAP_ENV="${BBMAP_ENV:-bbmap}"

echo "========================================="
echo "Merging Multiple Lanes Per Sample"
echo "========================================="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo

mkdir -p "$OUTPUT_DIR"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/reports"

# Initialize conda
if [ -f "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh" ]; then
    source "${CONDA_BASE:-$HOME/miniconda3}/etc/profile.d/conda.sh"
    conda activate "$BBMAP_ENV"
else
    echo "WARNING: Conda not found, assuming BBMap is in PATH"
fi

# Find all unique sample prefixes
declare -A samples

for file in "$INPUT_DIR"/*_L*_1.fq.gz; do
    [ -f "$file" ] || continue
    
    basename=$(basename "$file")
    # Extract sample name (everything before _DKDN or _L)
    sample_name=$(echo "$basename" | sed -E 's/(_DKDN.*|_L[0-9]+.*)//')
    
    samples["$sample_name"]=1
done

echo "Found ${#samples[@]} unique samples"
echo

# Create summary report
SUMMARY_REPORT="${OUTPUT_DIR}/reports/merge_summary.txt"
echo "Lane Merge Summary - $(date)" > "$SUMMARY_REPORT"
echo "=======================================" >> "$SUMMARY_REPORT"
echo "" >> "$SUMMARY_REPORT"

# Process each sample
TOTAL_SAMPLES=0
MERGED_SAMPLES=0
SYMLINKED_SAMPLES=0
FAILED_SAMPLES=0

for sample in "${!samples[@]}"; do
    ((TOTAL_SAMPLES++))
    echo "========================================="
    echo "Processing sample: $sample"
    echo "========================================="
    
    LOG_FILE="${OUTPUT_DIR}/logs/${sample}_merge.log"
    exec > >(tee -a "$LOG_FILE") 2>&1
    
    # Find all R1 and R2 files for this sample
    r1_files=($(find "$INPUT_DIR" -name "${sample}_*_L*_1.fq.gz" -o -name "${sample}_*_1.fq.gz" | sort))
    r2_files=($(find "$INPUT_DIR" -name "${sample}_*_L*_2.fq.gz" -o -name "${sample}_*_2.fq.gz" | sort))
    
    if [ ${#r1_files[@]} -eq 0 ]; then
        echo "  WARNING: No files found for $sample"
        echo "$sample: NO FILES FOUND" >> "$SUMMARY_REPORT"
        ((FAILED_SAMPLES++))
        continue
    fi
    
    echo "Found ${#r1_files[@]} lane(s):"
    for i in "${!r1_files[@]}"; do
        echo "  Lane $((i+1)):"
        echo "    R1: $(basename "${r1_files[$i]}")"
        echo "    R2: $(basename "${r2_files[$i]}")"
    done
    echo
    
    # Count reads in input files
    echo "Counting reads in input files..."
    declare -a r1_counts
    declare -a r2_counts
    total_r1=0
    total_r2=0
    
    for i in "${!r1_files[@]}"; do
        r1_count=$(zcat "${r1_files[$i]}" | wc -l | awk '{print $1/4}')
        r2_count=$(zcat "${r2_files[$i]}" | wc -l | awk '{print $1/4}')
        
        r1_counts[$i]=$r1_count
        r2_counts[$i]=$r2_count
        total_r1=$((total_r1 + r1_count))
        total_r2=$((total_r2 + r2_count))
        
        echo "  Lane $((i+1)): R1=$r1_count, R2=$r2_count reads"
        
        # Check if lane is synchronized
        if [ "$r1_count" -ne "$r2_count" ]; then
            echo "  ⚠️  WARNING: Lane $((i+1)) has unequal R1/R2 counts!"
        fi
    done
    
    echo "Total input reads: R1=$total_r1, R2=$total_r2"
    echo
    
    # Output files
    OUTPUT_R1="${OUTPUT_DIR}/${sample}_R1.fq.gz"
    OUTPUT_R2="${OUTPUT_DIR}/${sample}_R2.fq.gz"
    
    # Single lane - create symlinks
    if [ ${#r1_files[@]} -eq 1 ]; then
        echo "Single lane detected, creating symlinks..."
        ln -sf "$(realpath "${r1_files[0]}")" "$OUTPUT_R1"
        ln -sf "$(realpath "${r2_files[0]}")" "$OUTPUT_R2"
        
        echo "✅ Symlinked: $OUTPUT_R1"
        echo "$sample: SYMLINKED (single lane, $total_r1 reads)" >> "$SUMMARY_REPORT"
        ((SYMLINKED_SAMPLES++))
        
    # Multiple lanes - merge with BBMap reformat.sh
    else
        echo "Merging ${#r1_files[@]} lanes using BBMap reformat.sh..."
        
        # Create comma-separated file lists for BBMap
        r1_list=$(IFS=,; echo "${r1_files[*]}")
        r2_list=$(IFS=,; echo "${r2_files[*]}")
        
        # Temporary unrepaired output
        TEMP_R1="${OUTPUT_DIR}/.tmp_${sample}_R1.fq.gz"
        TEMP_R2="${OUTPUT_DIR}/.tmp_${sample}_R2.fq.gz"
        
        # Use BBMap's reformat.sh to merge and validate
        echo "Step 1: Merging lanes..."
        reformat.sh \
            in="$r1_list" \
            in2="$r2_list" \
            out="$TEMP_R1" \
            out2="$TEMP_R2" \
            threads="$THREADS" \
            ziplevel=6 \
            2>&1 | tee -a "${OUTPUT_DIR}/logs/${sample}_reformat.log"
        
        if [ $? -ne 0 ]; then
            echo "❌ ERROR: reformat.sh failed for $sample"
            echo "$sample: MERGE FAILED" >> "$SUMMARY_REPORT"
            ((FAILED_SAMPLES++))
            rm -f "$TEMP_R1" "$TEMP_R2"
            continue
        fi
        
        # Count reads in merged files
        echo
        echo "Step 2: Validating merged files..."
        merged_r1=$(zcat "$TEMP_R1" | wc -l | awk '{print $1/4}')
        merged_r2=$(zcat "$TEMP_R2" | wc -l | awk '{print $1/4}')
        
        echo "Merged read counts: R1=$merged_r1, R2=$merged_r2"
        
        # Check for read loss
        if [ "$merged_r1" -ne "$total_r1" ]; then
            echo "⚠️  WARNING: R1 read count mismatch! Input=$total_r1, Merged=$merged_r1"
            echo "  Lost reads: $((total_r1 - merged_r1))"
        fi
        
        if [ "$merged_r2" -ne "$total_r2" ]; then
            echo "⚠️  WARNING: R2 read count mismatch! Input=$total_r2, Merged=$merged_r2"
            echo "  Lost reads: $((total_r2 - merged_r2))"
        fi
        
        # Repair paired-end synchronization if needed
        if [ "$merged_r1" -ne "$merged_r2" ]; then
            echo
            echo "⚠️  Paired-end files are NOT synchronized!"
            echo "Step 3: Repairing with BBMap repair.sh..."
            
            repair.sh \
                in="$TEMP_R1" \
                in2="$TEMP_R2" \
                out="$OUTPUT_R1" \
                out2="$OUTPUT_R2" \
                outs="${OUTPUT_DIR}/${sample}_singletons.fq.gz" \
                threads="$THREADS" \
                repair=t \
                2>&1 | tee -a "${OUTPUT_DIR}/logs/${sample}_repair.log"
            
            if [ $? -ne 0 ]; then
                echo "❌ ERROR: repair.sh failed for $sample"
                echo "$sample: REPAIR FAILED" >> "$SUMMARY_REPORT"
                ((FAILED_SAMPLES++))
                rm -f "$TEMP_R1" "$TEMP_R2" "$OUTPUT_R1" "$OUTPUT_R2"
                continue
            fi
            
            # Count after repair
            repaired_r1=$(zcat "$OUTPUT_R1" | wc -l | awk '{print $1/4}')
            repaired_r2=$(zcat "$OUTPUT_R2" | wc -l | awk '{print $1/4}')
            singleton_count=$(zcat "${OUTPUT_DIR}/${sample}_singletons.fq.gz" | wc -l | awk '{print $1/4}')
            
            echo "After repair: R1=$repaired_r1, R2=$repaired_r2, Singletons=$singleton_count"
            
            if [ "$repaired_r1" -ne "$repaired_r2" ]; then
                echo "❌ ERROR: Files still not synchronized after repair!"
                echo "$sample: STILL UNSYNCHRONIZED" >> "$SUMMARY_REPORT"
                ((FAILED_SAMPLES++))
                continue
            fi
            
            # Clean up temp files
            rm -f "$TEMP_R1" "$TEMP_R2"
            
            echo "$sample: MERGED & REPAIRED (${#r1_files[@]} lanes, $repaired_r1 paired + $singleton_count singletons)" >> "$SUMMARY_REPORT"
            ((MERGED_SAMPLES++))
            
        else
            # Files are synchronized, just move them
            echo "✅ Paired-end files are synchronized!"
            mv "$TEMP_R1" "$OUTPUT_R1"
            mv "$TEMP_R2" "$OUTPUT_R2"
            
            echo "$sample: MERGED (${#r1_files[@]} lanes, $merged_r1 reads)" >> "$SUMMARY_REPORT"
            ((MERGED_SAMPLES++))
        fi
        
        echo "✅ Complete: $OUTPUT_R1"
    fi
    
    echo
done

# Final summary
echo "========================================="
echo "Lane Merging Complete!"
echo "========================================="
echo "Total samples processed: $TOTAL_SAMPLES"
echo "  Merged (multiple lanes): $MERGED_SAMPLES"
echo "  Symlinked (single lane): $SYMLINKED_SAMPLES"
echo "  Failed: $FAILED_SAMPLES"
echo
echo "Output directory: $OUTPUT_DIR"
echo "Summary report: $SUMMARY_REPORT"
echo "Individual logs: ${OUTPUT_DIR}/logs/"
echo

# Append final summary to report
echo "" >> "$SUMMARY_REPORT"
echo "=======================================" >> "$SUMMARY_REPORT"
echo "Summary:" >> "$SUMMARY_REPORT"
echo "  Total samples: $TOTAL_SAMPLES" >> "$SUMMARY_REPORT"
echo "  Merged: $MERGED_SAMPLES" >> "$SUMMARY_REPORT"
echo "  Symlinked: $SYMLINKED_SAMPLES" >> "$SUMMARY_REPORT"
echo "  Failed: $FAILED_SAMPLES" >> "$SUMMARY_REPORT"

if [ $FAILED_SAMPLES -gt 0 ]; then
    echo
    echo "⚠️  WARNING: $FAILED_SAMPLES sample(s) failed. Check logs in ${OUTPUT_DIR}/logs/"
    exit 1
fi

echo "Next steps:"
echo "1. Review the summary report: $SUMMARY_REPORT"
echo "2. Run pipeline: ./run_pipeline.sh -i $OUTPUT_DIR -o <output_dir>"
