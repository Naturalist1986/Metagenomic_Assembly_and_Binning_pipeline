#!/bin/bash
# generate_sample_sheet.sh - Create a sample sheet from existing fastq files

INPUT_DIR="${1:-/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/raw_fastqs}"
OUTPUT_FILE="${2:-${INPUT_DIR}/sample_sheet.tsv}"

echo "Generating sample sheet from files in: $INPUT_DIR"
echo "Output file: $OUTPUT_FILE"

# Create header
echo -e "Sample_Name\tTreatment\tR1_File\tR2_File" > "$OUTPUT_FILE"

# Find all R1 files and process
find "$INPUT_DIR" -name "*_1.fq.gz" -type f | sort | while read -r r1_file; do
    # Get the R2 file
    r2_file="${r1_file/_1.fq.gz/_2.fq.gz}"
    
    if [ ! -f "$r2_file" ]; then
        echo "WARNING: No R2 file found for $r1_file"
        continue
    fi
    
    # Get basenames
    r1_basename=$(basename "$r1_file")
    r2_basename=$(basename "$r2_file")
    
    # Extract sample name (everything before _DKDN or _1.fq.gz)
    sample_name=$(basename "$r1_file" | sed -E 's/_DKDN.*/_/; s/_1\.fq\.gz$//')
    
    # Extract treatment from the beginning of the filename
    treatment=""
    if [[ "$sample_name" =~ ^(carK|carR|ces|hok|mtz|RH) ]]; then
        treatment="${BASH_REMATCH[1]}"
    else
        treatment="unknown"
    fi
    
    # Use full filename (without extension) as unique sample name to avoid duplicates
    unique_sample_name=$(basename "$r1_file" .fq.gz | sed 's/_1$//')
    
    # Add to sample sheet
    echo -e "${unique_sample_name}\t${treatment}\t${r1_basename}\t${r2_basename}" >> "$OUTPUT_FILE"
done

# Display summary
echo ""
echo "Sample sheet created: $OUTPUT_FILE"
echo ""
echo "Summary:"
echo "--------"
total_samples=$(grep -v "^Sample_Name" "$OUTPUT_FILE" | wc -l)
echo "Total samples: $total_samples"
echo ""
echo "Samples per treatment:"
for treatment in carK carR ces hok mtz RH unknown; do
    count=$(grep -v "^Sample_Name" "$OUTPUT_FILE" | grep -c "	${treatment}	" || echo "0")
    if [ $count -gt 0 ]; then
        echo "  $treatment: $count"
    fi
done
echo ""
echo "First 5 entries:"
head -6 "$OUTPUT_FILE" | column -t
echo ""
echo "To review the full sample sheet:"
echo "  less $OUTPUT_FILE"
echo "  # or"
echo "  cat $OUTPUT_FILE | column -t | less -S"
