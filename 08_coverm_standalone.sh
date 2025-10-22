#!/bin/bash
#SBATCH --job-name=coverm_standalone
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00

# 08_coverm_standalone.sh - CoverM for existing refined bins and combined reads

# Configuration
BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
BIN_REFINEMENT_DIR="${BASE_DIR}/bin_refinement"
COMBINED_READS_DIR="${BASE_DIR}/combined_reads"
OUTPUT_DIR="${BASE_DIR}/coverm"
LOG_DIR="${BASE_DIR}/logs/coverm"
CHECKPOINT_DIR="${BASE_DIR}/checkpoints/coverm"

# Conda setup
CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"

# Create directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$CHECKPOINT_DIR"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_DIR}/coverm_${SLURM_ARRAY_TASK_ID}.log"
}

# Initialize conda
init_conda() {
    if command -v module >/dev/null 2>&1; then
        module purge
    fi
    source ${CONDA_BASE}/etc/profile.d/conda.sh
}

# Activate conda environment
activate_env() {
    local env_name="$1"
    conda activate "$env_name" || {
        log "ERROR: Failed to activate $env_name"
        exit 1
    }
    log "Activated conda env: $env_name"
}

# Discover all samples with refined bins
discover_samples() {
    local sample_file="${BASE_DIR}/coverm_samples.txt"
    > "$sample_file"
    
    log "Discovering samples with refined bins..."
    
    # Find all treatment directories
    for treatment_dir in "$BIN_REFINEMENT_DIR"/*; do
        if [ ! -d "$treatment_dir" ]; then
            continue
        fi
        
        local treatment=$(basename "$treatment_dir")
        
        # Find all sample directories within treatment
        for sample_dir in "$treatment_dir"/*; do
            if [ ! -d "$sample_dir" ]; then
                continue
            fi
            
            local sample=$(basename "$sample_dir")
            local bins_dir="${sample_dir}/metawrap_50_10_bins"
            
            # Check if refined bins exist
            if [ -d "$bins_dir" ] && [ "$(ls -A "$bins_dir"/*.fa 2>/dev/null)" ]; then
                local bin_count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
                
                # Check if combined reads exist
                local read1="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_1.fq.gz"
                local read2="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_2.fq.gz"
                
                if [ -f "$read1" ] && [ -f "$read2" ]; then
                    echo "${treatment}|${sample}|${bins_dir}|${read1}|${read2}" >> "$sample_file"
                    log "  Found: ${treatment}/${sample} ($bin_count bins)"
                else
                    log "  WARNING: Bins found but reads missing for ${treatment}/${sample}"
                fi
            fi
        done
    done
    
    local total=$(wc -l < "$sample_file")
    log "Discovered $total samples ready for CoverM analysis"
    echo "$sample_file"
}

# Get sample info by array index
get_sample_by_index() {
    local index="$1"
    local sample_file="$2"
    local line_num=$((index + 1))
    sed -n "${line_num}p" "$sample_file"
}

# Validate bins
validate_bins() {
    local bins_dir="$1"
    local min_size=10000
    
    local valid_bins=0
    for bin_file in "$bins_dir"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
        local contig_count=$(grep -c "^>" "$bin_file")
        
        if [ $bin_size -ge $min_size ] && [ $contig_count -ge 1 ]; then
            ((valid_bins++))
        fi
    done
    
    echo "$valid_bins"
}

# Run CoverM
run_coverm() {
    local treatment="$1"
    local sample="$2"
    local bins_dir="$3"
    local read1="$4"
    local read2="$5"
    local output_dir="$6"
    
    log "Running CoverM for ${treatment}/${sample}..."
    
    # Validate bins
    local valid_bins=$(validate_bins "$bins_dir")
    if [ $valid_bins -eq 0 ]; then
        log "ERROR: No valid bins found in $bins_dir"
        return 1
    fi
    
    log "  Analyzing $valid_bins valid bins"
    log "  Bins: $bins_dir"
    log "  R1: $read1"
    log "  R2: $read2"
    
    # Activate CoverM environment
    activate_env checkm
    
    if ! command -v coverm &> /dev/null; then
        log "ERROR: CoverM not available"
        conda deactivate
        return 1
    fi
    
    # Create output directory
    mkdir -p "$output_dir/mapping"
    
    # Get list of valid bin files
    local bin_files=()
    for bin_file in "$bins_dir"/*.fa; do
        if [ -f "$bin_file" ]; then
            local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
            if [ $bin_size -ge 10000 ]; then
                bin_files+=("$bin_file")
            fi
        fi
    done
    
    # Run CoverM
    log "Running CoverM genome mode..."
    
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
        --threads ${SLURM_CPUS_PER_TASK:-16} \
        --bam-file-cache-directory "${output_dir}/mapping" \
        2>&1 | tee "${LOG_DIR}/${treatment}_${sample}_coverm.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/abundance.tsv" ]; then
        log "CoverM completed successfully"
        return 0
    else
        log "ERROR: CoverM failed (exit code: $exit_code)"
        return 1
    fi
}

# Generate summary statistics
generate_summary() {
    local treatment="$1"
    local sample="$2"
    local output_dir="$3"
    local abundance_file="${output_dir}/abundance.tsv"
    local summary_file="${output_dir}/abundance_summary.txt"
    
    if [ ! -f "$abundance_file" ]; then
        log "ERROR: Abundance file not found"
        return 1
    fi
    
    log "Generating summary for ${treatment}/${sample}..."
    
    cat > "$summary_file" << EOF
CoverM Abundance Summary
========================

Sample: ${sample}
Treatment: ${treatment}
Date: $(date)
Bins Source: ${BIN_REFINEMENT_DIR}/${treatment}/${sample}/metawrap_50_10_bins

EOF
    
    awk -F'\t' '
        BEGIN {
            print "Top 10 Most Abundant Bins:"
            print "-------------------------"
        }
        NR == 1 {
            for (i = 2; i <= NF; i++) {
                if ($i ~ /Relative Abundance/) rel_col = i
                if ($i ~ /Mean/) mean_col = i
                if ($i ~ /TPM/) tpm_col = i
                if ($i ~ /Covered Fraction/) cov_col = i
            }
        }
        NR > 1 && $1 != "" && $1 != "unmapped" {
            bins[NR-1] = $1
            rel_abund[NR-1] = (rel_col ? $rel_col : 0)
            mean_cov[NR-1] = (mean_col ? $mean_col : 0)
            tpm_val[NR-1] = (tpm_col ? $tpm_col : 0)
            cov_frac[NR-1] = (cov_col ? $cov_col : 0)
            if ($rel_col > 0) detected++
            total++
        }
        END {
            # Sort by relative abundance
            for (i = 1; i <= total; i++) {
                for (j = i+1; j <= total; j++) {
                    if (rel_abund[j] > rel_abund[i]) {
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
            for (i = 1; i <= total && count < 10; i++) {
                if (rel_abund[i] > 0) {
                    rel_percent = rel_abund[i] * 100
                    if (rel_percent > 100) rel_percent = rel_abund[i]
                    printf "%2d. %-40s %6.2f%% (Cov: %.2fx, TPM: %.2f)\n", 
                        ++count, bins[i], rel_percent, mean_cov[i], tpm_val[i]
                }
            }
            
            print ""
            print "Summary:"
            print "--------"
            printf "Total bins: %d\n", total
            printf "Bins detected: %d (%.1f%%)\n", detected, (total>0 ? detected*100/total : 0)
        }
    ' "$abundance_file" >> "$summary_file"
    
    log "Summary created: $summary_file"
}

# Main execution
main() {
    log "========================================="
    log "CoverM Standalone Analysis"
    log "========================================="
    
    # Initialize
    init_conda
    
    # Discover samples (only once, by first task)
    local sample_file="${BASE_DIR}/coverm_samples.txt"
    
    if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 0 ] || [ ! -f "$sample_file" ]; then
        sample_file=$(discover_samples)
    fi
    
    # Wait for sample file to be ready
    while [ ! -f "$sample_file" ]; do
        sleep 2
    done
    
    # Get sample for this array task
    local task_id=${SLURM_ARRAY_TASK_ID:-0}
    local sample_info=$(get_sample_by_index "$task_id" "$sample_file")
    
    if [ -z "$sample_info" ]; then
        log "No sample found for task ID $task_id"
        exit 0
    fi
    
    # Parse sample info
    IFS='|' read -r treatment sample bins_dir read1 read2 <<< "$sample_info"
    
    log "Processing: ${treatment}/${sample}"
    
    # Check if already completed
    local output_dir="${OUTPUT_DIR}/${treatment}/${sample}"
    if [ -f "${CHECKPOINT_DIR}/${treatment}_${sample}_complete" ]; then
        log "Already completed, skipping..."
        exit 0
    fi
    
    mkdir -p "$output_dir"
    
    # Run CoverM
    if run_coverm "$treatment" "$sample" "$bins_dir" "$read1" "$read2" "$output_dir"; then
        generate_summary "$treatment" "$sample" "$output_dir"
        
        # Create checkpoint
        mkdir -p "$(dirname "${CHECKPOINT_DIR}/${treatment}_${sample}_complete")"
        touch "${CHECKPOINT_DIR}/${treatment}_${sample}_complete"
        
        log "========================================="
        log "CoverM completed for ${treatment}/${sample}"
        log "========================================="
    else
        log "ERROR: CoverM failed for ${treatment}/${sample}"
        exit 1
    fi
}

# Run main
main
