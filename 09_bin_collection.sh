#!/bin/bash
#SBATCH --job-name=bin_collection
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00

# 09_bin_collection.sh - High-quality genome bin collection and categorization

# Source configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"

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

log "====== Starting Bin Collection for $SAMPLE_NAME ($TREATMENT) ======"

# Check if stage already completed
if check_sample_checkpoint "$SAMPLE_NAME" "bin_collection"; then
    log "Bin collection already completed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Define quality categories
declare -A QUALITY_THRESHOLDS=(
    ["high_min_comp"]=90
    ["high_max_cont"]=5
    ["medium_min_comp"]=50
    ["medium_max_cont"]=10
)

# Function to find best available bins
find_best_bins_source() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Finding best available bins for $sample_name..."
    
    # Priority order: MAGpurify > Reassembly > Refinement
    local sources=(
        "${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}/purified_bins:MAGpurify purified bins"
        "${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}/reassembled_bins:Reassembled bins"
        "${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}/metawrap_50_10_bins:Refined bins"
    )
    
    for source_info in "${sources[@]}"; do
        IFS=':' read -r bins_dir source_name <<< "$source_info"
        
        if [ -d "$bins_dir" ] && [ "$(ls -A "$bins_dir"/*.fa 2>/dev/null)" ]; then
            local bin_count=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
            log "  Found $bin_count bins from: $source_name"
            log "  Using bins from: $bins_dir"
            echo "$bins_dir|$source_name"
            return 0
        fi
    done
    
    log "ERROR: No bins found for $sample_name"
    return 1
}

# Function to load CheckM2 quality data
load_checkm2_data() {
    local sample_name="$1"
    local treatment="$2"
    local checkm2_file="${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}/quality_report.tsv"
    
    log "Loading CheckM2 quality data for $sample_name..."
    
    if [ ! -f "$checkm2_file" ]; then
        log "WARNING: CheckM2 quality report not found: $checkm2_file"
        return 1
    fi
    
    # Create temporary quality lookup file
    local quality_lookup="${TEMP_DIR}/quality_lookup.tsv"
    
    # Process CheckM2 file and create lookup
    awk -F'\t' '
        NR > 1 && $1 != "" && $2 != "" && $3 != "" {
            bin_name = $1
            gsub(/\.fa$/, "", bin_name)  # Remove .fa extension
            completeness = $2
            contamination = $3
            strain_het = $4
            gc = $7
            genome_size = $8
            n50 = $9
            taxonomy = $12
            
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                   bin_name, completeness, contamination, strain_het, gc, genome_size, n50, taxonomy
        }
    ' "$checkm2_file" > "$quality_lookup"
    
    local quality_bins=$(wc -l < "$quality_lookup")
    log "  Loaded quality data for $quality_bins bins"
    
    if [ $quality_bins -eq 0 ]; then
        log "WARNING: No quality data found"
        return 1
    fi
    
    echo "$quality_lookup"
    return 0
}

# Function to load CoverM abundance data
load_coverm_data() {
    local sample_name="$1"
    local treatment="$2"
    local coverm_file="${OUTPUT_DIR}/coverm/${treatment}/${sample_name}/abundance.tsv"
    
    log "Loading CoverM abundance data for $sample_name..."
    
    if [ ! -f "$coverm_file" ]; then
        log "WARNING: CoverM abundance file not found: $coverm_file"
        return 1
    fi
    
    # Create temporary abundance lookup file
    local abundance_lookup="${TEMP_DIR}/abundance_lookup.tsv"
    
    # Process CoverM file and create lookup
    awk -F'\t' '
        NR == 1 {
            # Find column indices
            for (i = 2; i <= NF; i++) {
                if ($i ~ /Relative Abundance/) rel_col = i
                if ($i ~ /Mean/) mean_col = i
                if ($i ~ /TPM/) tpm_col = i
                if ($i ~ /Covered Fraction/) cov_col = i
            }
        }
        NR > 1 && $1 != "" && $1 != "unmapped" {
            bin_name = $1
            gsub(/\.fa$/, "", bin_name)  # Remove .fa extension
            
            rel_abund = (rel_col ? $rel_col : 0)
            mean_cov = (mean_col ? $mean_col : 0)
            tpm = (tpm_col ? $tpm_col : 0)
            cov_frac = (cov_col ? $cov_col : 0)
            
            printf "%s\t%s\t%s\t%s\t%s\n", bin_name, rel_abund, mean_cov, tpm, cov_frac
        }
    ' "$coverm_file" > "$abundance_lookup"
    
    local abundance_bins=$(wc -l < "$abundance_lookup")
    log "  Loaded abundance data for $abundance_bins bins"
    
    if [ $abundance_bins -eq 0 ]; then
        log "WARNING: No abundance data found"
        return 1
    fi
    
    echo "$abundance_lookup"
    return 0
}

# Function to categorize and collect bins
categorize_and_collect_bins() {
    local sample_name="$1"
    local treatment="$2"
    local bins_dir="$3"
    local bins_source="$4"
    local quality_lookup="$5"
    local abundance_lookup="$6"
    local output_dir="$7"
    
    log "Categorizing and collecting bins for $sample_name..."
    
    # Create output directories
    mkdir -p "${output_dir}/high_quality"
    mkdir -p "${output_dir}/medium_quality"
    mkdir -p "${output_dir}/low_quality"
    mkdir -p "${output_dir}/all_bins"
    
    # Create summary files
    local all_bins_summary="${output_dir}/all_bins_summary.tsv"
    local hq_summary="${output_dir}/high_quality_bins_summary.tsv"
    local mq_summary="${output_dir}/medium_quality_bins_summary.tsv"
    local lq_summary="${output_dir}/low_quality_bins_summary.tsv"
    
    # Header for summary files
    local header="Sample\tBin\tBin_Source\tCompleteness\tContamination\tStrain_Het\tGC\tGenome_Size\tN50\tTaxonomy\tQuality_Category\tRelative_Abundance\tMean_Coverage\tTPM\tCovered_Fraction"
    
    echo -e "$header" > "$all_bins_summary"
    echo -e "$header" > "$hq_summary"
    echo -e "$header" > "$mq_summary"
    echo -e "$header" > "$lq_summary"
    
    # Initialize counters
    local total_bins=0
    local hq_bins=0
    local mq_bins=0
    local lq_bins=0
    
    # Process each bin
    for bin_file in "$bins_dir"/*.fa; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi
        
        local bin_filename=$(basename "$bin_file")
        local bin_name=$(basename "$bin_file" .fa)
        ((total_bins++))
        
        log "  Processing bin: $bin_name"
        
        # Get quality data
        local quality_data="NA\tNA\tNA\tNA\tNA\tNA\tNA"
        local quality_category="Unknown"
        
        if [ -f "$quality_lookup" ]; then
            local quality_line=$(grep "^${bin_name}\t" "$quality_lookup" | head -1)
            if [ -n "$quality_line" ]; then
                quality_data="$quality_line"
                
                # Extract completeness and contamination for categorization
                local completeness=$(echo "$quality_line" | cut -f2)
                local contamination=$(echo "$quality_line" | cut -f3)
                
                # Categorize quality
                if [ "$completeness" != "NA" ] && [ "$contamination" != "NA" ]; then
                    if (( $(echo "$completeness >= ${QUALITY_THRESHOLDS[high_min_comp]}" | bc -l) )) && \
                       (( $(echo "$contamination <= ${QUALITY_THRESHOLDS[high_max_cont]}" | bc -l) )); then
                        quality_category="High_Quality"
                        ((hq_bins++))
                    elif (( $(echo "$completeness >= ${QUALITY_THRESHOLDS[medium_min_comp]}" | bc -l) )) && \
                         (( $(echo "$contamination <= ${QUALITY_THRESHOLDS[medium_max_cont]}" | bc -l) )); then
                        quality_category="Medium_Quality"
                        ((mq_bins++))
                    else
                        quality_category="Low_Quality"
                        ((lq_bins++))
                    fi
                fi
            fi
        fi
        
        # Get abundance data
        local abundance_data="NA\tNA\tNA\tNA"
        if [ -f "$abundance_lookup" ]; then
            local abundance_line=$(grep "^${bin_name}\t" "$abundance_lookup" | head -1)
            if [ -n "$abundance_line" ]; then
                abundance_data=$(echo "$abundance_line" | cut -f2-5)
            fi
        fi
        
        # Create new bin name with sample prefix
        local new_bin_name="${treatment}_${sample_name}_${bin_name}.fa"
        
        # Copy bin to appropriate directories
        case "$quality_category" in
            "High_Quality")
                cp "$bin_file" "${output_dir}/high_quality/${new_bin_name}"
                ;;
            "Medium_Quality")
                cp "$bin_file" "${output_dir}/medium_quality/${new_bin_name}"
                ;;
            *)
                cp "$bin_file" "${output_dir}/low_quality/${new_bin_name}"
                ;;
        esac
        
        # Always copy to all_bins
        cp "$bin_file" "${output_dir}/all_bins/${new_bin_name}"
        
        # Create summary line
        local summary_line="${sample_name}\t${bin_name}\t${bins_source}\t${quality_data}\t${quality_category}\t${abundance_data}"
        
        # Add to appropriate summary files
        echo -e "$summary_line" >> "$all_bins_summary"
        
        case "$quality_category" in
            "High_Quality")
                echo -e "$summary_line" >> "$hq_summary"
                ;;
            "Medium_Quality")
                echo -e "$summary_line" >> "$mq_summary"
                ;;
            *)
                echo -e "$summary_line" >> "$lq_summary"
                ;;
        esac
        
        log "    Category: $quality_category"
    done
    
    log "Bin collection complete for $sample_name:"
    log "  Total bins: $total_bins"
    log "  High-quality bins: $hq_bins"
    log "  Medium-quality bins: $mq_bins"
    log "  Low-quality bins: $lq_bins"
    
    # Store results for later use
    echo "$total_bins|$hq_bins|$mq_bins|$lq_bins"
    return 0
}

# Function to create merged FASTA files
create_merged_fastas() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="$3"
    
    log "Creating merged FASTA files for $sample_name..."
    
    for quality in high_quality medium_quality low_quality all_bins; do
        local merged_file="${output_dir}/${quality}_merged.fa"
        > "$merged_file"
        
        local count=0
        for bin_file in "${output_dir}/${quality}"/*.fa; do
            if [ -f "$bin_file" ] && [[ "$(basename "$bin_file")" != *"_merged.fa" ]]; then
                # Add sample prefix to sequence headers
                awk -v sample="$sample_name" -v bin="$(basename "$bin_file" .fa)" '
                    /^>/ { 
                        print ">" sample "_" bin "_" substr($0, 2)
                        next 
                    } 
                    { print }
                ' "$bin_file" >> "$merged_file"
                ((count++))
            fi
        done
        
        log "  Created ${quality}_merged.fa with $count bins"
    done
}

# Function to create HTML report
create_html_report() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="$3"
    local total_bins="$4"
    local hq_bins="$5"
    local mq_bins="$6"
    local lq_bins="$7"
    
    log "Creating HTML report for $sample_name..."
    
    local html_file="${output_dir}/bin_collection_report.html"
    local all_bins_summary="${output_dir}/all_bins_summary.tsv"
    
    cat > "$html_file" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Bin Collection Report - $sample_name ($treatment)</title>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            box-shadow: 0 0 10px rgba(0,0,0,0.1);
        }
        h1 { 
            color: #2c3e50;
            text-align: center;
            padding-bottom: 20px;
            border-bottom: 2px solid #3498db;
        }
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
            padding: 20px;
            background-color: #f8f9fa;
            border-radius: 5px;
        }
        .stat-box {
            text-align: center;
            padding: 15px;
            background-color: white;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .stat-number {
            font-size: 2.5em;
            font-weight: bold;
            color: #3498db;
        }
        .stat-label {
            color: #7f8c8d;
            margin-top: 5px;
        }
        .high { color: #27ae60; font-weight: bold; }
        .medium { color: #f39c12; font-weight: bold; }
        .low { color: #e74c3c; font-weight: bold; }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin-top: 20px; 
        }
        th, td { 
            text-align: left; 
            padding: 8px; 
            border: 1px solid #ddd; 
            font-size: 0.9em;
        }
        th { 
            background-color: #3498db;
            color: white;
            font-weight: bold;
            cursor: pointer;
        }
        tr:nth-child(even) { background-color: #f9f9f9; }
        tr:hover { background-color: #f1f1f1; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Bin Collection Report</h1>
        <h2>$sample_name ($treatment)</h2>
        
        <div class="summary-stats">
            <div class="stat-box">
                <div class="stat-number">$total_bins</div>
                <div class="stat-label">Total Bins</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">$hq_bins</div>
                <div class="stat-label">High Quality</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">$mq_bins</div>
                <div class="stat-label">Medium Quality</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">$lq_bins</div>
                <div class="stat-label">Low Quality</div>
            </div>
        </div>
        
        <h3>Quality Definitions</h3>
        <ul>
            <li><strong>High Quality:</strong> ≥90% complete, ≤5% contamination</li>
            <li><strong>Medium Quality:</strong> ≥50% complete, ≤10% contamination</li>
            <li><strong>Low Quality:</strong> All other bins</li>
        </ul>
        
        <table id="binTable">
            <thead>
                <tr>
                    <th>Bin Name</th>
                    <th>Completeness (%)</th>
                    <th>Contamination (%)</th>
                    <th>GC (%)</th>
                    <th>Size (bp)</th>
                    <th>N50</th>
                    <th>Quality</th>
                    <th>Rel. Abundance (%)</th>
                    <th>Mean Coverage</th>
                </tr>
            </thead>
            <tbody>
EOF
    
    # Add table rows from summary file
    if [ -f "$all_bins_summary" ]; then
        tail -n +2 "$all_bins_summary" | while IFS=$'\t' read -r sample bin bin_source completeness contamination strain_het gc genome_size n50 taxonomy quality_category rel_abund mean_cov tpm cov_frac; do
            # Format numbers
            local comp_fmt=$(printf "%.1f" "$completeness" 2>/dev/null || echo "$completeness")
            local cont_fmt=$(printf "%.1f" "$contamination" 2>/dev/null || echo "$contamination")
            local gc_fmt=$(printf "%.1f" "$gc" 2>/dev/null || echo "$gc")
            local size_fmt=$(printf "%'d" "$genome_size" 2>/dev/null || echo "$genome_size")
            local n50_fmt=$(printf "%'d" "$n50" 2>/dev/null || echo "$n50")
            local abund_fmt=$(printf "%.3f" "$rel_abund" 2>/dev/null || echo "$rel_abund")
            local cov_fmt=$(printf "%.2f" "$mean_cov" 2>/dev/null || echo "$mean_cov")
            
            # Quality class
            local quality_class=""
            case "$quality_category" in
                "High_Quality") quality_class="high" ;;
                "Medium_Quality") quality_class="medium" ;;
                *) quality_class="low" ;;
            esac
            
            cat >> "$html_file" << TABLEROW
                <tr>
                    <td>$bin</td>
                    <td>$comp_fmt</td>
                    <td>$cont_fmt</td>
                    <td>$gc_fmt</td>
                    <td>$size_fmt</td>
                    <td>$n50_fmt</td>
                    <td><span class="$quality_class">$quality_category</span></td>
                    <td>$abund_fmt</td>
                    <td>$cov_fmt</td>
                </tr>
TABLEROW
        done
    fi
    
    cat >> "$html_file" << EOF
            </tbody>
        </table>
    </div>
    
    <script>
        // Add sorting functionality
        document.querySelectorAll('th').forEach(header => {
            header.addEventListener('click', () => {
                const table = document.getElementById('binTable');
                const rows = Array.from(table.rows).slice(1);
                const index = Array.from(header.parentNode.children).indexOf(header);
                const ascending = header.classList.contains('asc');
                
                rows.sort((a, b) => {
                    const aVal = a.cells[index].textContent.replace(/[^0-9.-]/g, '') || a.cells[index].textContent;
                    const bVal = b.cells[index].textContent.replace(/[^0-9.-]/g, '') || b.cells[index].textContent;
                    
                    if (!isNaN(aVal) && !isNaN(bVal)) {
                        return ascending ? aVal - bVal : bVal - aVal;
                    }
                    return ascending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                });
                
                header.classList.toggle('asc');
                const tbody = table.querySelector('tbody');
                tbody.innerHTML = '';
                rows.forEach(row => tbody.appendChild(row));
            });
        });
    </script>
</body>
</html>
EOF
    
    log "HTML report created: $html_file"
}

# Function to create summary report
create_summary_report() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="$3"
    local bins_source="$4"
    local total_bins="$5"
    local hq_bins="$6"
    local mq_bins="$7"
    local lq_bins="$8"
    
    local summary_file="${output_dir}/collection_summary.txt"
    
    cat > "$summary_file" << EOF
Bin Collection Summary for $sample_name
======================================

Date: $(date)
Sample: $sample_name
Treatment: $treatment
Bins Source: $bins_source

Quality Distribution:
- High-quality bins (≥90% complete, ≤5% contamination): $hq_bins
- Medium-quality bins (≥50% complete, ≤10% contamination): $mq_bins
- Low-quality bins: $lq_bins
- Total bins: $total_bins

Quality Percentages:
EOF
    
    if [ $total_bins -gt 0 ]; then
        local hq_pct=$(echo "scale=1; $hq_bins * 100 / $total_bins" | bc -l)
        local mq_pct=$(echo "scale=1; $mq_bins * 100 / $total_bins" | bc -l)
        local lq_pct=$(echo "scale=1; $lq_bins * 100 / $total_bins" | bc -l)
        
        echo "- High-quality: ${hq_pct}%" >> "$summary_file"
        echo "- Medium-quality: ${mq_pct}%" >> "$summary_file"
        echo "- Low-quality: ${lq_pct}%" >> "$summary_file"
    fi
    
    cat >> "$summary_file" << EOF

Output Files:
- All bins: ${output_dir}/all_bins/
- High-quality bins: ${output_dir}/high_quality/
- Medium-quality bins: ${output_dir}/medium_quality/
- Low-quality bins: ${output_dir}/low_quality/
- Summary files: ${output_dir}/*_summary.tsv
- HTML report: ${output_dir}/bin_collection_report.html
- Merged FASTA files: ${output_dir}/*_merged.fa

Data Integration:
- CheckM2 quality assessment: $([ -f "${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}/quality_report.tsv" ] && echo "Available" || echo "Not available")
- CoverM abundance data: $([ -f "${OUTPUT_DIR}/coverm/${treatment}/${sample_name}/abundance.tsv" ] && echo "Available" || echo "Not available")

EOF
    
    log "Summary report created: $summary_file"
}

# Main processing function
stage_bin_collection() {
    local sample_name="$1"
    local treatment="$2"
    
    log "Running bin collection for $sample_name ($treatment)"
    
    local output_dir="${OUTPUT_DIR}/bin_collection/${treatment}/${sample_name}"
    mkdir -p "$output_dir"
    
    # Check if already processed
    if [ -f "${output_dir}/collection_complete.flag" ]; then
        log "Sample $sample_name already processed, skipping..."
        return 0
    fi
    
    # Find best available bins
    local bins_info=$(find_best_bins_source "$sample_name" "$treatment")
    if [ $? -ne 0 ]; then
        log "ERROR: No bins found for $sample_name"
        return 1
    fi
    
    IFS='|' read -r bins_dir bins_source <<< "$bins_info"
    
    # Load quality and abundance data
    local quality_lookup=$(load_checkm2_data "$sample_name" "$treatment")
    local abundance_lookup=$(load_coverm_data "$sample_name" "$treatment")
    
    # Categorize and collect bins
    local results=$(categorize_and_collect_bins "$sample_name" "$treatment" "$bins_dir" "$bins_source" "$quality_lookup" "$abundance_lookup" "$output_dir")
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to categorize and collect bins"
        return 1
    fi
    
    IFS='|' read -r total_bins hq_bins mq_bins lq_bins <<< "$results"
    
    # Create merged FASTA files
    create_merged_fastas "$sample_name" "$treatment" "$output_dir"
    
    # Create reports
    create_html_report "$sample_name" "$treatment" "$output_dir" "$total_bins" "$hq_bins" "$mq_bins" "$lq_bins"
    create_summary_report "$sample_name" "$treatment" "$output_dir" "$bins_source" "$total_bins" "$hq_bins" "$mq_bins" "$lq_bins"
    
    # Create completion flag
    touch "${output_dir}/collection_complete.flag"
    
    log "Bin collection completed for $sample_name"
    return 0
}

# Validation function
validate_bin_collection() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/bin_collection/${treatment}/${sample_name}"
    
    # Check completion flag
    if [ ! -f "${output_dir}/collection_complete.flag" ]; then
        return 1
    fi
    
    # Check that output directories exist
    for dir in all_bins high_quality medium_quality low_quality; do
        if [ ! -d "${output_dir}/${dir}" ]; then
            return 1
        fi
    done
    
    # Check that we have at least one bin
    if [ ! "$(ls -A "${output_dir}/all_bins/"*.fa 2>/dev/null)" ]; then
        return 1
    fi
    
    # Check that summary files exist
    if [ ! -f "${output_dir}/all_bins_summary.tsv" ]; then
        return 1
    fi
    
    return 0
}

# Run the bin collection stage
if stage_bin_collection "$SAMPLE_NAME" "$TREATMENT"; then
    # Validate results
    if validate_bin_collection "$SAMPLE_NAME" "$TREATMENT"; then
        create_sample_checkpoint "$SAMPLE_NAME" "bin_collection"
        log "====== Bin collection completed for $SAMPLE_NAME ======"
    else
        log "ERROR: Bin collection validation failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log "ERROR: Bin collection stage failed for $SAMPLE_NAME"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"