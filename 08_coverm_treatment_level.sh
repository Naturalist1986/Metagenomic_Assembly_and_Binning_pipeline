#!/bin/bash
#SBATCH --job-name=coverm_treatment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap/logs/coverm_treatment/coverm_%A_%a.out
#SBATCH --error=/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap/logs/coverm_treatment/coverm_%A_%a.err

# 08_coverm_treatment_level.sh - Map ALL bins in treatment to each sample's reads

# Configuration
BASE_DIR="/sci/backup/aerez/aerez/moshea/Efrat_Metagenomes_Novogene/new_metawrap"
BIN_REFINEMENT_DIR="${BASE_DIR}/bin_refinement"
COMBINED_READS_DIR="${BASE_DIR}/combined_reads"
OUTPUT_DIR="${BASE_DIR}/coverm_treatment_level"
LOG_DIR="${BASE_DIR}/logs/coverm_treatment"
CHECKPOINT_DIR="${BASE_DIR}/checkpoints/coverm_treatment"
TREATMENT_LIST="${BASE_DIR}/treatment_list.txt"

# Conda setup
CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"

# Create directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$CHECKPOINT_DIR"

# Get task ID
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
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

# Get treatment info by array index
get_treatment_by_index() {
    local index="$1"
    local line_num=$((index + 1))
    sed -n "${line_num}p" "$TREATMENT_LIST" 2>/dev/null
}

# Collect all bins for a treatment
collect_treatment_bins() {
    local treatment="$1"
    local treatment_bins_dir="$2"
    
    log "Collecting all bins for treatment: $treatment"
    mkdir -p "$treatment_bins_dir"
    
    local total_bins=0
    local treatment_dir="${BIN_REFINEMENT_DIR}/${treatment}"
    
    if [ ! -d "$treatment_dir" ]; then
        log "ERROR: Treatment directory not found: $treatment_dir"
        return 1
    fi
    
    # Collect bins from all samples in this treatment
    for sample_dir in "$treatment_dir"/*; do
        if [ ! -d "$sample_dir" ]; then
            continue
        fi
        
        local sample=$(basename "$sample_dir")
        local bins_dir="${sample_dir}/metawrap_50_10_bins"
        
        if [ ! -d "$bins_dir" ]; then
            continue
        fi
        
        # Copy bins with sample prefix to avoid name conflicts
        for bin_file in "$bins_dir"/*.fa; do
            if [ -f "$bin_file" ]; then
                local bin_name=$(basename "$bin_file" .fa)
                local new_name="${sample}.${bin_name}.fa"
                
                # Only include bins >= 10kb
                local bin_size=$(grep -v "^>" "$bin_file" | tr -d '\n' | wc -c)
                if [ $bin_size -ge 10000 ]; then
                    cp "$bin_file" "${treatment_bins_dir}/${new_name}"
                    ((total_bins++))
                fi
            fi
        done
    done
    
    log "Collected $total_bins bins from treatment $treatment"
    echo "$total_bins"
}

# Get all samples for a treatment
get_treatment_samples() {
    local treatment="$1"
    local samples_file="$2"
    
    > "$samples_file"
    
    local treatment_dir="${BIN_REFINEMENT_DIR}/${treatment}"
    
    for sample_dir in "$treatment_dir"/*; do
        if [ ! -d "$sample_dir" ]; then
            continue
        fi
        
        local sample=$(basename "$sample_dir")
        local read1="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_1.fq.gz"
        local read2="${COMBINED_READS_DIR}/${treatment}/${sample}/${sample}_2.fq.gz"
        
        if [ -f "$read1" ] && [ -f "$read2" ]; then
            echo "${sample}|${read1}|${read2}" >> "$samples_file"
        fi
    done
    
    local sample_count=$(wc -l < "$samples_file")
    echo "$sample_count"
}

# Run CoverM for all bins against one sample
run_coverm_for_sample() {
    local treatment="$1"
    local sample="$2"
    local bins_dir="$3"
    local read1="$4"
    local read2="$5"
    local output_dir="$6"
    
    log "Running CoverM: ALL ${treatment} bins vs ${sample} reads"
    
    # Activate CoverM environment
    activate_env checkm
    
    if ! command -v coverm &> /dev/null; then
        log "ERROR: CoverM not available"
        conda deactivate
        return 1
    fi
    
    # Get all bin files
    local bin_files=("$bins_dir"/*.fa)
    local bin_count=${#bin_files[@]}
    
    log "  Mapping $bin_count bins to sample $sample"
    log "  R1: $read1"
    log "  R2: $read2"
    
    # Create output directory
    mkdir -p "$output_dir/mapping"
    
    # Run CoverM
    coverm genome \
        --genome-fasta-files "${bin_files[@]}" \
        --coupled "$read1" "$read2" \
        --output-file "${output_dir}/${sample}_abundance.tsv" \
        --output-format dense \
        --min-read-aligned-percent 0.75 \
        --min-read-percent-identity 0.95 \
        --min-covered-fraction 0 \
        --contig-end-exclusion 75 \
        --trim-min 0.05 \
        --trim-max 0.95 \
        --proper-pairs-only \
        --methods relative_abundance mean trimmed_mean covered_fraction reads_per_base rpkm tpm \
        --threads ${SLURM_CPUS_PER_TASK:-32} \
        --bam-file-cache-directory "${output_dir}/mapping/${sample}" \
        2>&1 | tee -a "${output_dir}/${sample}_coverm.log"
    
    local exit_code=${PIPESTATUS[0]}
    conda deactivate
    
    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/${sample}_abundance.tsv" ]; then
        log "CoverM completed for sample $sample"
        return 0
    else
        log "ERROR: CoverM failed for sample $sample (exit code: $exit_code)"
        return 1
    fi
}

# Generate per-sample summary
generate_sample_summary() {
    local treatment="$1"
    local sample="$2"
    local output_dir="$3"
    local abundance_file="${output_dir}/${sample}_abundance.tsv"
    local summary_file="${output_dir}/${sample}_summary.txt"
    
    if [ ! -f "$abundance_file" ]; then
        return 1
    fi
    
    cat > "$summary_file" << EOF
CoverM Treatment-Level Analysis
================================

Treatment: ${treatment}
Sample: ${sample}
Date: $(date)

Mapping: ALL ${treatment} bins -> ${sample} reads

Top 10 Most Abundant Bins in ${sample}:
---------------------------------------

EOF
    
    awk -F'\t' '
        NR == 1 {
            for (i = 2; i <= NF; i++) {
                if ($i ~ /Relative Abundance/) rel_col = i
                if ($i ~ /Mean/) mean_col = i
            }
            next
        }
        $1 != "" && $1 != "unmapped" {
            bins[NR-1] = $1
            rel_abund[NR-1] = (rel_col ? $rel_col : 0)
            mean_cov[NR-1] = (mean_col ? $mean_col : 0)
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
                    }
                }
            }
            
            # Print top 10
            count = 0
            for (i = 1; i <= total && count < 10; i++) {
                if (rel_abund[i] > 0) {
                    rel_percent = rel_abund[i] * 100
                    if (rel_percent > 100) rel_percent = rel_abund[i]
                    printf "%2d. %-50s %6.2f%% (Cov: %.2fx)\n", 
                        ++count, bins[i], rel_percent, mean_cov[i]
                }
            }
            
            print ""
            printf "Total bins in treatment: %d\n", total
            printf "Bins detected in this sample: %d (%.1f%%)\n", detected, (total>0 ? detected*100/total : 0)
        }
    ' "$abundance_file" >> "$summary_file"
    
    log "Summary created: $summary_file"
}

# Create treatment-wide summary combining all samples
create_treatment_summary() {
    local treatment="$1"
    local output_dir="$2"
    local summary_file="${output_dir}/TREATMENT_SUMMARY.txt"
    
    log "Creating treatment-wide summary for $treatment..."
    
    cat > "$summary_file" << EOF
========================================
TREATMENT-WIDE COVERM SUMMARY
========================================

Treatment: ${treatment}
Date: $(date)

Analysis: All bins from all samples in ${treatment} 
          mapped to each individual sample's reads

========================================

EOF
    
    # Find all abundance files for this treatment
    local abundance_files=("${output_dir}"/*_abundance.tsv)
    local sample_count=${#abundance_files[@]}
    
    echo "Samples analyzed: $sample_count" >> "$summary_file"
    echo "" >> "$summary_file"
    
    # Count total unique bins
    if [ -f "${abundance_files[0]}" ]; then
        local total_bins=$(tail -n +2 "${abundance_files[0]}" | grep -v "^unmapped" | wc -l)
        echo "Total bins in treatment: $total_bins" >> "$summary_file"
        echo "" >> "$summary_file"
    fi
    
    # Summary table
    echo "Per-Sample Detection Summary:" >> "$summary_file"
    echo "$(printf '=%.0s' {1..70})" >> "$summary_file"
    printf "%-20s %10s %10s %10s\n" "Sample" "Total_Bins" "Detected" "Rate%" >> "$summary_file"
    echo "$(printf '=%.0s' {1..70})" >> "$summary_file"
    
    for abund_file in "${abundance_files[@]}"; do
        if [ ! -f "$abund_file" ]; then
            continue
        fi
        
        local sample=$(basename "$abund_file" "_abundance.tsv")
        local total=$(tail -n +2 "$abund_file" | grep -v "^unmapped" | wc -l)
        local detected=$(tail -n +2 "$abund_file" | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)if($i~/Relative Abundance/)col=i; next} $col>0' | wc -l)
        local rate=$(echo "scale=1; ($detected * 100) / $total" | bc -l 2>/dev/null || echo "0")
        
        printf "%-20s %10d %10d %9.1f%%\n" "$sample" "$total" "$detected" "$rate" >> "$summary_file"
    done
    
    echo "" >> "$summary_file"
    echo "Individual sample summaries available:" >> "$summary_file"
    for abund_file in "${abundance_files[@]}"; do
        local sample=$(basename "$abund_file" "_abundance.tsv")
        echo "  - ${sample}_summary.txt" >> "$summary_file"
    done
    
    log "Treatment summary created: $summary_file"
}

# Main execution
main() {
    log "========================================="
    log "Treatment-Level CoverM - Task ${TASK_ID}"
    log "========================================="
    
    # Check if treatment list exists
    if [ ! -f "$TREATMENT_LIST" ]; then
        log "ERROR: Treatment list not found: $TREATMENT_LIST"
        log "Please run the submission script first"
        exit 1
    fi
    
    # Initialize
    init_conda
    
    # Get treatment for this task
    local treatment=$(get_treatment_by_index "$TASK_ID")
    
    if [ -z "$treatment" ]; then
        log "No treatment found for task ID $TASK_ID"
        exit 0
    fi
    
    log "Processing treatment: $treatment"
    
    # Check if already completed
    if [ -f "${CHECKPOINT_DIR}/${treatment}_complete" ]; then
        log "Treatment $treatment already completed, skipping..."
        exit 0
    fi
    
    # Create output directory for this treatment
    local output_dir="${OUTPUT_DIR}/${treatment}"
    mkdir -p "$output_dir"
    
    # Step 1: Collect all bins from this treatment
    local treatment_bins_dir="${output_dir}/all_bins"
    local total_bins=$(collect_treatment_bins "$treatment" "$treatment_bins_dir")
    
    if [ "$total_bins" -eq 0 ]; then
        log "ERROR: No bins found for treatment $treatment"
        exit 1
    fi
    
    # Step 2: Get all samples for this treatment
    local samples_file="${output_dir}/samples.txt"
    local sample_count=$(get_treatment_samples "$treatment" "$samples_file")
    
    if [ "$sample_count" -eq 0 ]; then
        log "ERROR: No samples found for treatment $treatment"
        exit 1
    fi
    
    log "Will map $total_bins bins to $sample_count samples"
    
    # Step 3: Run CoverM for each sample
    local success_count=0
    local fail_count=0
    
    while IFS='|' read -r sample read1 read2; do
        log "----------------------------------------"
        log "Processing sample: $sample"
        
        if run_coverm_for_sample "$treatment" "$sample" "$treatment_bins_dir" "$read1" "$read2" "$output_dir"; then
            generate_sample_summary "$treatment" "$sample" "$output_dir"
            ((success_count++))
        else
            log "WARNING: Failed to process sample $sample"
            ((fail_count++))
        fi
    done < "$samples_file"
    
    log "----------------------------------------"
    log "Treatment $treatment processing complete"
    log "  Successful: $success_count samples"
    log "  Failed: $fail_count samples"
    
    # Step 4: Create treatment-wide summary
    create_treatment_summary "$treatment" "$output_dir"
    
    # Create checkpoint if all samples succeeded
    if [ $fail_count -eq 0 ]; then
        touch "${CHECKPOINT_DIR}/${treatment}_complete"
        log "========================================="
        log "TREATMENT $treatment COMPLETED"
        log "Results: ${output_dir}/"
        log "Summary: ${output_dir}/TREATMENT_SUMMARY.txt"
        log "========================================="
    else
        log "WARNING: Some samples failed for treatment $treatment"
        exit 1
    fi
}

# Run main
main
