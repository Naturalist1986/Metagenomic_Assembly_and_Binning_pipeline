#!/bin/bash
#SBATCH --job-name=coverm
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00

# 08_coverm_refined_bins.sh - CoverM abundance calculation using refined bins only

# Source configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"

# Get sample info from array task ID
SAMPLE_INFO=$(get_sample_info_by_index $SLURM_ARRAY_TASK_ID)
if [ -z "$SAMPLE_INFO" ]; then
    log "No sample found for array index $SLURM_ARRAY_TASK_ID"
    exit 0
fi

# Parse sample information
IFS='|' read -r SAMPLE_NAME TREATMENT _ _ <<< "$SAMPLE_INFO"
export SAMPLE_NAME TREATMENT

# Initialize
init_conda
create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"
TEMP_DIR=$(setup_temp_dir)

log "====== Starting CoverM Abundance Calculation for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "coverm"; then
    log "CoverM already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Function to validate bins for CoverM
validate_bins_for_coverm() {
    local bins_dir="$1"
    local min_size="${2:-10000}"  # 10KB minimum
    local min_contigs="${3:-1}"   # 1 contig minimum
    
    log "Validating bins for CoverM analysis..."
    
    if [ ! -d "$bins_dir" ]; then
        log "ERROR: Bins directory not found: $bins_dir"
        return 1
    fi
    
    local valid_bins=0
    local total_bins=0
    local skipped_bins=0
    
    for bin_file in "$bins_dir"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        ((total_bins++))
        
        local bin_name=$(basename "$bin_file")
        local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
        local contig_count=$(grep -c "^>" "$bin_file")
        
        if [ $bin_size -ge $min_size ] && [ $contig_count -ge $min_contigs ]; then
            ((valid_bins++))
            log "  Valid: $bin_name ($bin_size bp, $contig_count contigs)"
        else
            ((skipped_bins++))
            log "  Skipped: $bin_name (too small: $bin_size bp, $contig_count contigs)"
        fi
    done
    
    log "Bin validation complete: $valid_bins/$total_bins valid, $skipped_bins skipped"
    
    if [ $valid_bins -eq 0 ]; then
        log "ERROR: No valid bins for CoverM analysis"
        return 1
    fi
    
    echo "$valid_bins"
    return 0
}

# Function to get refined bins (ONLY from bin_refinement)
get_refined_bins() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Looking for refined bins for $sample_name..."
    
    local bins_dir="${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}/metawrap_50_10_bins"
    
    if [ -d "$bins_dir" ] && [ "$(ls -A "$bins_dir"/*.fa 2>/dev/null)" ]; then
        local bin_count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
        log "  Found $bin_count refined bins in: $bins_dir"
        echo "$bins_dir"
        return 0
    else
        log "ERROR: No refined bins found for $sample_name in $bins_dir"
        return 1
    fi
}

# Function to run CoverM abundance calculation
run_coverm_analysis() {
    local sample_name="$1"
    local treatment="$2"
    local bins_dir="$3"
    local read1="$4"
    local read2="$5"
    local output_dir="$6"
    
    log "Running CoverM abundance calculation for $sample_name..."
    log "  Using refined bins from: $bins_dir"
    
    # Activate CoverM environment
    activate_env checkm
    
    # Check if CoverM is available
    if ! command -v coverm &> /dev/null; then
        log "ERROR: CoverM not available in environment"
        conda deactivate
        return 1
    fi
    
    # Create mapping directory
    local mapping_dir="${output_dir}/mapping"
    mkdir -p "$mapping_dir"
    
    # Count valid bins
    local valid_bins=$(validate_bins_for_coverm "$bins_dir")
    if [ $? -ne 0 ]; then
        log "ERROR: No valid bins for CoverM analysis"
        conda deactivate
        return 1
    fi
    
    # Get list of valid bin files
    local bin_files=()
    for bin_file in "$bins_dir"/*.fa; do
        if [ -f "$bin_file" ]; then
            local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
            local contig_count=$(grep -c "^>" "$bin_file")
            
            # Only include bins that meet minimum thresholds
            if [ $bin_size -ge 10000 ] && [ $contig_count -ge 1 ]; then
                bin_files+=("$bin_file")
            fi
        fi
    done
    
    if [ ${#bin_files[@]} -eq 0 ]; then
        log "ERROR: No valid bin files found"
        conda deactivate
        return 1
    fi
    
    log "Calculating abundance for ${#bin_files[@]} bins..."
    
    # Run CoverM genome mode
    coverm genome \
        --genome-fasta-files "${bin_files[@]}" \
        --coupled "$read1" "$read2" \
        --output-file "${output_dir}/abundance.tsv" \
        --output-format dense \
        --min-read-aligned-percent 0.75 \
        --min-read-percent-identity 0.95 \
        --min-covered-fraction 0 \
        --contig-end-exclusion 75 \
        --trim-min 0.05 \
        --trim-max 0.95 \
        --proper-pairs-only \
        --methods relative_abundance mean trimmed_mean covered_fraction reads_per_base rpkm tpm \
        --threads $SLURM_CPUS_PER_TASK \
        --bam-file-cache-directory "$mapping_dir" \
        2>&1 | tee "${LOG_DIR}/${treatment}/${sample_name}_coverm.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/abundance.tsv" ]; then
        log "CoverM analysis completed successfully"
        return 0
    else
        log "ERROR: CoverM analysis failed (exit code: $exit_code)"
        return 1
    fi
}

# Function to generate abundance statistics
generate_abundance_stats() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="$3"
    local abundance_file="${output_dir}/abundance.tsv"
    local stats_file="${output_dir}/abundance_summary.txt"
    
    log "Generating abundance statistics for $sample_name..."
    
    if [ ! -f "$abundance_file" ]; then
        log "ERROR: Abundance file not found: $abundance_file"
        return 1
    fi
    
    cat > "$stats_file" << EOF
CoverM Abundance Summary for $sample_name
=========================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment
Bins Source: Refined bins (bin_refinement)

EOF
    
    # Parse abundance file and create summary
    awk -F'\t' '
        BEGIN {
            print "Top 10 Most Abundant Bins (by relative abundance):"
            print "------------------------------------------------"
        }
        NR == 1 {
            # Find column indices for different methods
            for (i = 2; i <= NF; i++) {
                if ($i ~ /Relative Abundance/) rel_abund_col = i
                if ($i ~ /Mean/) mean_col = i
                if ($i ~ /TPM/) tpm_col = i
                if ($i ~ /Covered Fraction/) cov_frac_col = i
            }
            if (!rel_abund_col) {
                print "ERROR: Could not find Relative Abundance column"
                exit 1
            }
        }
        NR > 1 && $1 != "" && $1 != "unmapped" {
            # Store bin data
            bins[NR-1] = $1
            if (rel_abund_col) rel_abund[NR-1] = $rel_abund_col
            if (mean_col) mean_cov[NR-1] = $mean_col
            if (tpm_col) tpm_val[NR-1] = $tpm_col
            if (cov_frac_col) cov_frac[NR-1] = $cov_frac_col
            
            # Count bins with coverage
            if ($rel_abund_col > 0) detected++
            total_bins++
        }
        END {
            # Sort by relative abundance (simple bubble sort for top 10)
            for (i = 1; i <= total_bins; i++) {
                for (j = i+1; j <= total_bins; j++) {
                    if (rel_abund[j] > rel_abund[i]) {
                        # Swap all arrays
                        tmp = bins[i]; bins[i] = bins[j]; bins[j] = tmp
                        tmp = rel_abund[i]; rel_abund[i] = rel_abund[j]; rel_abund[j] = tmp
                        tmp = mean_cov[i]; mean_cov[i] = mean_cov[j]; mean_cov[j] = tmp
                        tmp = tpm_val[i]; tpm_val[i] = tpm_val[j]; tpm_val[j] = tmp
                        tmp = cov_frac[i]; cov_frac[i] = cov_frac[j]; cov_frac[j] = tmp
                    }
                }
            }
            
            # Print top 10
            count = 0
            for (i = 1; i <= total_bins && count < 10; i++) {
                if (rel_abund[i] > 0) {
                    # Convert to percentage if needed
                    rel_percent = rel_abund[i] * 100
                    if (rel_percent > 100) {
                        # Raw value, likely already percentage
                        rel_percent = rel_abund[i]
                    }
                    printf "%2d. %-40s %6.2f%% (Mean: %.2f, TPM: %.2f, Covered: %.2f%%)\n", 
                        ++count, bins[i], rel_percent, mean_cov[i], tpm_val[i], cov_frac[i]*100
                }
            }
            
            # Print summary
            print ""
            print "Summary Statistics:"
            print "-----------------"
            printf "Total bins: %d\n", total_bins
            printf "Bins with coverage (>0%%): %d\n", detected
            if (total_bins > 0) {
                printf "Detection rate: %.1f%%\n", (detected/total_bins)*100
            }
            
            # Calculate total relative abundance
            total_rel = 0
            for (i = 1; i <= total_bins; i++) {
                total_rel += rel_abund[i]
            }
            
            # Handle percentage vs fraction
            if (total_rel > 10) {
                printf "Total relative abundance: %.2f%%\n", total_rel
            } else {
                printf "Total relative abundance: %.2f (%.2f%%)\n", total_rel, total_rel * 100
            }
        }
    ' "$abundance_file" >> "$stats_file"
    
    log "Abundance statistics created: $stats_file"
}

# Function to extract relative abundance for treatment summary
extract_relative_abundance() {
    local abundance_file="$1"
    local sample_name="$2"
    local output_file="$3"
    
    awk -F'\t' -v sample="$sample_name" '
        NR == 1 {
            # Find relative abundance column
            for (i = 2; i <= NF; i++) {
                if ($i ~ /Relative Abundance/) {
                    rel_col = i
                    break
                }
            }
            if (!rel_col) {
                print "ERROR: Could not find Relative Abundance column" > "/dev/stderr"
                exit 1
            }
            print "Bin_Name\t" sample
        }
        NR > 1 && $1 != "" && $1 != "unmapped" {
            # Extract bin name (remove path and .fa extension)
            bin_name = $1
            gsub(/^.*\//, "", bin_name)
            gsub(/\.fa$/, "", bin_name)
            
            # Get relative abundance as percentage
            rel_abund = $rel_col
            if (rel_abund < 10) {
                # Likely a fraction, convert to percentage
                rel_abund = rel_abund * 100
            }
            
            printf "%s\t%.4f\n", bin_name, rel_abund
        }
    ' "$abundance_file" > "$output_file"
}

# Main processing function
stage_coverm() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Running CoverM abundance calculation for $sample_name ($treatment)"
    
    local output_dir="${OUTPUT_DIR}/coverm/${treatment}/${sample_name}"
    local quality_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    
    mkdir -p "$output_dir"
    
    # Check if already processed
    if [ -f "${output_dir}/abundance.tsv" ] && [ -f "${output_dir}/coverm_complete.flag" ]; then
        log "Sample $sample_name already processed, skipping..."
        return 0
    fi
    
    # Get refined bins (ONLY from bin_refinement)
    local bins_dir=$(get_refined_bins "$sample_name" "$treatment")
    if [ $? -ne 0 ]; then
        log "ERROR: No refined bins found for $sample_name"
        return 1
    fi
    
    # Check for quality-filtered reads
    local read1="${quality_dir}/filtered_1.fastq.gz"
    local read2="${quality_dir}/filtered_2.fastq.gz"
    
    if [ ! -f "$read1" ] || [ ! -f "$read2" ]; then
        log "ERROR: Missing quality-filtered reads for $sample_name"
        log "  Expected: $read1 and $read2"
        return 1
    fi
    
    # Run CoverM analysis
    if run_coverm_analysis "$sample_name" "$treatment" "$bins_dir" "$read1" "$read2" "$output_dir"; then
        log "CoverM analysis completed successfully"
    else
        log "ERROR: CoverM analysis failed"
        return 1
    fi
    
    # Generate statistics
    generate_abundance_stats "$sample_name" "$treatment" "$output_dir"
    
    # Extract relative abundance for treatment summary
    local temp_abundance="${output_dir}/relative_abundance_${sample_name}.tsv"
    extract_relative_abundance "${output_dir}/abundance.tsv" "$sample_name" "$temp_abundance"
    
    # Create completion flag
    touch "${output_dir}/coverm_complete.flag"
    
    log "CoverM abundance calculation completed for $sample_name"
    return 0
}

# Validation function
validate_coverm() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/coverm/${treatment}/${sample_name}"
    
    # Check completion flag
    if [ ! -f "${output_dir}/coverm_complete.flag" ]; then
        return 1
    fi
    
    # Check abundance file exists and is not empty
    if [ ! -f "${output_dir}/abundance.tsv" ] || [ ! -s "${output_dir}/abundance.tsv" ]; then
        return 1
    fi
    
    # Check that abundance file has more than just header
    local lines=$(wc -l < "${output_dir}/abundance.tsv")
    if [ $lines -le 1 ]; then
        return 1
    fi
    
    # Check that summary file exists
    if [ ! -f "${output_dir}/abundance_summary.txt" ]; then
        return 1
    fi
    
    return 0
}

# Run the CoverM stage
if stage_coverm "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_coverm "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "coverm"
        log "====== CoverM abundance calculation completed for $SAMPLE_NAME ======"
    else
        log "ERROR: CoverM validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: CoverM stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"