#!/bin/bash
# setup_abundance_analysis.sh - Setup script for running abundance analysis

echo "=========================================="
echo "Abundance Analysis Setup"
echo "=========================================="
echo ""

# Detect script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_SCRIPT_DIR="$SCRIPT_DIR"

echo "Script directory: $SCRIPT_DIR"
echo ""

# Check for config file
if [ ! -f "${SCRIPT_DIR}/00_config_utilities.sh" ]; then
    echo "ERROR: Cannot find 00_config_utilities.sh in $SCRIPT_DIR"
    exit 1
fi

source "${SCRIPT_DIR}/00_config_utilities.sh"

# Prompt for or detect OUTPUT_DIR
if [ -z "$OUTPUT_DIR" ]; then
    echo "OUTPUT_DIR is not set."
    echo ""
    echo "Please enter the full path to your MetaWRAP output directory"
    echo "(This should contain bin_refinement/, quality_filtering/, etc.)"
    echo ""
    read -p "OUTPUT_DIR: " OUTPUT_DIR
    export OUTPUT_DIR
fi

echo "OUTPUT_DIR: $OUTPUT_DIR"
echo ""

# Validate OUTPUT_DIR
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Directory does not exist: $OUTPUT_DIR"
    exit 1
fi

# Check for bin_refinement directory
if [ ! -d "${OUTPUT_DIR}/bin_refinement" ]; then
    echo "WARNING: bin_refinement directory not found in $OUTPUT_DIR"
    echo "Expected: ${OUTPUT_DIR}/bin_refinement"
    echo ""
    read -p "Continue anyway? (y/n): " continue_anyway
    if [ "$continue_anyway" != "y" ]; then
        exit 1
    fi
fi

# Set WORK_DIR
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"
echo "WORK_DIR: $WORK_DIR"
echo ""

# Check if sample info exists
SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"
TREATMENTS_FILE="${WORK_DIR}/treatments.txt"

if [ -f "$SAMPLE_INFO_FILE" ]; then
    echo "✓ Found existing sample info file"
    
    # Count samples
    sample_count=$(wc -l < "$SAMPLE_INFO_FILE")
    echo "  Samples: $sample_count"
    
    # List treatments
    if [ -f "$TREATMENTS_FILE" ]; then
        echo "  Treatments:"
        while read treatment; do
            sample_count=$(grep -c "|${treatment}|" "$SAMPLE_INFO_FILE")
            echo "    - $treatment ($sample_count samples)"
        done < "$TREATMENTS_FILE"
    fi
    echo ""
else
    echo "⚠ Sample info file not found"
    echo ""
    echo "The pipeline needs sample information to run."
    echo "This can be created from:"
    echo "  1. A sample sheet (Excel/TSV file)"
    echo "  2. Auto-discovery from FASTQ files"
    echo ""
    
    read -p "Do you have a sample sheet? (y/n): " has_sheet
    
    if [ "$has_sheet" = "y" ]; then
        echo ""
        read -p "Enter path to sample sheet: " SAMPLE_SHEET
        export SAMPLE_SHEET
        
        if [ ! -f "$SAMPLE_SHEET" ]; then
            echo "ERROR: Sample sheet not found: $SAMPLE_SHEET"
            exit 1
        fi
        
        # Need INPUT_DIR for file paths
        echo ""
        echo "Sample sheet contains relative paths to FASTQ files."
        read -p "Enter INPUT_DIR (directory containing FASTQ files): " INPUT_DIR
        export INPUT_DIR
        
        if [ ! -d "$INPUT_DIR" ]; then
            echo "ERROR: Directory does not exist: $INPUT_DIR"
            exit 1
        fi
    else
        echo ""
        read -p "Enter directory containing FASTQ files (for auto-discovery): " INPUT_DIR
        export INPUT_DIR
        
        if [ ! -d "$INPUT_DIR" ]; then
            echo "ERROR: Directory does not exist: $INPUT_DIR"
            exit 1
        fi
    fi
    
    echo ""
    echo "Initializing sample information..."
    init_sample_info
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to initialize sample information"
        exit 1
    fi
    
    echo "✓ Sample information initialized"
    echo ""
    
    # Display results
    sample_count=$(wc -l < "$SAMPLE_INFO_FILE")
    echo "  Samples found: $sample_count"
    
    if [ -f "$TREATMENTS_FILE" ]; then
        echo "  Treatments:"
        while read treatment; do
            sample_count=$(grep -c "|${treatment}|" "$SAMPLE_INFO_FILE")
            echo "    - $treatment ($sample_count samples)"
        done < "$TREATMENTS_FILE"
    fi
    echo ""
fi

# Verify bins exist
echo "Checking for refined bins..."
bin_count=0
treatment_count=0

if [ -d "${OUTPUT_DIR}/bin_refinement" ]; then
    for treatment_dir in "${OUTPUT_DIR}/bin_refinement"/*; do
        if [ -d "$treatment_dir" ]; then
            treatment=$(basename "$treatment_dir")
            for sample_dir in "$treatment_dir"/*; do
                if [ -d "$sample_dir" ]; then
                    bins_dir="${sample_dir}/metawrap_50_10_bins"
                    if [ -d "$bins_dir" ]; then
                        count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
                        if [ $count -gt 0 ]; then
                            ((bin_count += count))
                            ((treatment_count++))
                        fi
                    fi
                fi
            done
        fi
    done
fi

if [ $bin_count -gt 0 ]; then
    echo "✓ Found $bin_count refined bins across $treatment_count samples"
else
    echo "⚠ No refined bins found in ${OUTPUT_DIR}/bin_refinement"
    echo "  Make sure bin refinement has been completed"
fi
echo ""

# Verify quality-filtered reads exist
echo "Checking for quality-filtered reads..."
reads_count=0

if [ -d "${OUTPUT_DIR}/quality_filtering" ]; then
    for treatment_dir in "${OUTPUT_DIR}/quality_filtering"/*; do
        if [ -d "$treatment_dir" ]; then
            for sample_dir in "$treatment_dir"/*; do
                if [ -d "$sample_dir" ]; then
                    if [ -f "${sample_dir}/filtered_1.fastq.gz" ] && [ -f "${sample_dir}/filtered_2.fastq.gz" ]; then
                        ((reads_count++))
                    fi
                fi
            done
        fi
    done
fi

if [ $reads_count -gt 0 ]; then
    echo "✓ Found quality-filtered reads for $reads_count samples"
else
    echo "⚠ No quality-filtered reads found in ${OUTPUT_DIR}/quality_filtering"
    echo "  Make sure quality filtering has been completed"
fi
echo ""

# Create environment file
ENV_FILE="${SCRIPT_DIR}/abundance_env.sh"
cat > "$ENV_FILE" << EOF
# Auto-generated environment file for abundance analysis
# Source this file before running abundance scripts:
# source ${ENV_FILE}

export PIPELINE_SCRIPT_DIR="${SCRIPT_DIR}"
export OUTPUT_DIR="${OUTPUT_DIR}"
export WORK_DIR="${WORK_DIR}"
export INPUT_DIR="${INPUT_DIR:-}"
export SAMPLE_SHEET="${SAMPLE_SHEET:-}"
export CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"
EOF

echo "✓ Created environment file: $ENV_FILE"
echo ""

# Summary
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "To run abundance analysis:"
echo ""
echo "1. Load the environment:"
echo "   source $ENV_FILE"
echo ""
echo "2. Submit CoverM job:"
echo "   sbatch 08_coverm_refined_bins.sh"
echo ""
echo "3. After completion, create summary:"
echo "   sbatch 09_coverm_treatment_summary.sh"
echo ""
echo "Or for MetaWRAP quant_bins:"
echo "   sbatch 08_metawrap_quant_bins.sh"
echo "   sbatch 09_metawrap_treatment_summary.sh"
echo ""
echo "Configuration saved to: $ENV_FILE"
echo "=========================================="
