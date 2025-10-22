#!/bin/bash
#SBATCH --job-name=treatment_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# 10_treatment_analysis.sh - Treatment-level analysis combining samples within treatments

# Source configuration and utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"

# Initialize
init_conda
TEMP_DIR=$(setup_temp_dir)

# Get treatment from command line argument or environment
TREATMENT="${1:-${SLURM_ARRAY_TASK_ID}}"
if [ -z "$TREATMENT" ]; then
    echo "Usage: $0 <treatment_name>"
    echo "Available treatments: $(get_treatments)"
    exit 1
fi

# Validate treatment
VALID_TREATMENTS=($(get_treatments))
if [[ ! " ${VALID_TREATMENTS[@]} " =~ " ${TREATMENT} " ]]; then
    echo "ERROR: Invalid treatment '$TREATMENT'"
    echo "Available treatments: ${VALID_TREATMENTS[@]}"
    exit 1
fi

export TREATMENT

log "====== Starting Treatment-Level Analysis for $TREATMENT ======"

# Check if stage already completed
if check_treatment_checkpoint "$TREATMENT" "treatment_analysis"; then
    log "Treatment analysis already completed for $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 0
fi

# Function to wait for all samples to complete bin collection
wait_for_sample_completion() {
    local treatment="$1"
    local max_wait=3600  # 1 hour
    local wait_time=0
    local check_interval=60  # 1 minute
    
    log "Waiting for all samples in $treatment to complete bin collection..."
    
    # Get all samples for this treatment
    local samples=($(get_samples_for_treatment "$treatment"))
    local total_samples=${#samples[@]}
    
    if [ $total_samples -eq 0 ]; then
        log "ERROR: No samples found for treatment $treatment"
        return 1
    fi
    
    log "Found $total_samples samples for treatment $treatment"
    
    while [ $wait_time -lt $max_wait ]; do
        local completed_samples=0
        local pending_samples=()
        
        for sample in "${samples[@]}"; do
            if check_sample_checkpoint "$sample" "bin_collection"; then
                ((completed_samples++))
            else
                pending_samples+=("$sample")
            fi
        done
        
        log "Bin collection status: $completed_samples/$total_samples samples completed"
        
        if [ $completed_samples -eq $total_samples ]; then
            log "All samples completed bin collection!"
            return 0
        fi
        
        if [ ${#pending_samples[@]} -le 5 ]; then
            log "Pending samples: ${pending_samples[*]}"
        fi
        
        sleep $check_interval
        wait_time=$((wait_time + check_interval))
    done
    
    log "WARNING: Timeout waiting for sample completion. Proceeding with available data."
    return 1
}

# Function to get completed samples
get_completed_samples() {
    local treatment="$1"
    local samples=($(get_samples_for_treatment "$treatment"))
    local completed_samples=()
    
    for sample in "${samples[@]}"; do
        if check_sample_checkpoint "$sample" "bin_collection"; then
            completed_samples+=("$sample")
        fi
    done
    
    echo "${completed_samples[@]}"
}

# Function to collect all bins from samples
collect_treatment_bins() {
    local treatment="$1"
    local samples=($2)
    local output_dir="$3"
    
    log "Collecting bins from all samples in $treatment..."
    
    # Create output directories
    mkdir -p "${output_dir}/all_samples"
    mkdir -p "${output_dir}/high_quality"
    mkdir -p "${output_dir}/medium_quality"
    mkdir -p "${output_dir}/low_quality"
    
    # Initialize counters
    local total_bins=0
    local hq_bins=0
    local mq_bins=0
    local lq_bins=0
    
    # Create combined summary file
    local combined_summary="${output_dir}/treatment_bins_summary.tsv"
    local header="Sample\tBin\tBin_Source\tCompleteness\tContamination\tStrain_Het\tGC\tGenome_Size\tN50\tTaxonomy\tQuality_Category\tRelative_Abundance\tMean_Coverage\tTPM\tCovered_Fraction"
    echo -e "$header" > "$combined_summary"
    
    # Process each sample
    for sample in "${samples[@]}"; do
        log "  Processing sample: $sample"
        
        local sample_collection_dir="${OUTPUT_DIR}/bin_collection/${treatment}/${sample}"
        
        if [ ! -d "$sample_collection_dir" ]; then
            log "    WARNING: Collection directory not found for $sample"
            continue
        fi
        
        # Copy bins from each quality category
        for quality in high_quality medium_quality low_quality; do
            local sample_quality_dir="${sample_collection_dir}/${quality}"
            if [ -d "$sample_quality_dir" ]; then
                local bin_count=$(ls -1 "$sample_quality_dir"/*.fa 2>/dev/null | wc -l)
                if [ $bin_count -gt 0 ]; then
                    cp "$sample_quality_dir"/*.fa "${output_dir}/${quality}/" 2>/dev/null
                    log "    Copied $bin_count bins from $quality category"
                    
                    case "$quality" in
                        "high_quality") hq_bins=$((hq_bins + bin_count)) ;;
                        "medium_quality") mq_bins=$((mq_bins + bin_count)) ;;
                        "low_quality") lq_bins=$((lq_bins + bin_count)) ;;
                    esac
                fi
            fi
        done
        
        # Copy all bins
        local sample_all_dir="${sample_collection_dir}/all_bins"
        if [ -d "$sample_all_dir" ]; then
            local all_bin_count=$(ls -1 "$sample_all_dir"/*.fa 2>/dev/null | wc -l)
            if [ $all_bin_count -gt 0 ]; then
                cp "$sample_all_dir"/*.fa "${output_dir}/all_samples/" 2>/dev/null
                total_bins=$((total_bins + all_bin_count))
            fi
        fi
        
        # Append to combined summary
        local sample_summary="${sample_collection_dir}/all_bins_summary.tsv"
        if [ -f "$sample_summary" ]; then
            tail -n +2 "$sample_summary" >> "$combined_summary"
        fi
    done
    
    log "Treatment bin collection complete:"
    log "  Total bins: $total_bins"
    log "  High-quality bins: $hq_bins"
    log "  Medium-quality bins: $mq_bins"
    log "  Low-quality bins: $lq_bins"
    
    echo "$total_bins|$hq_bins|$mq_bins|$lq_bins"
    return 0
}

# Function to perform bin dereplication (optional)
perform_bin_dereplication() {
    local treatment="$1"
    local bins_dir="$2"
    local output_dir="$3"
    local similarity_threshold="${4:-95}"
    
    log "Performing bin dereplication for $treatment (${similarity_threshold}% similarity)..."
    
    # Check if dRep is available
    activate_env metagenome_assembly
    
    if ! command -v dRep &> /dev/null; then
        log "WARNING: dRep not available, skipping dereplication"
        log "  Copying all bins without dereplication"
        cp "$bins_dir"/*.fa "$output_dir/" 2>/dev/null
        conda deactivate
        return 0
    fi
    
    # Create dRep work directory
    local drep_work="${TEMP_DIR}/drep_${treatment}"
    mkdir -p "$drep_work"
    
    # Count input bins
    local input_bins=$(ls -1 "$bins_dir"/*.fa 2>/dev/null | wc -l)
    if [ $input_bins -eq 0 ]; then
        log "No bins to dereplicate"
        conda deactivate
        return 0
    fi
    
    log "  Dereplicating $input_bins bins..."
    
    # Run dRep
    dRep dereplicate "$drep_work" \
        -g "$bins_dir"/*.fa \
        -comp 50 \
        -con 25 \
        -sa "$similarity_threshold" \
        -nc 0.30 \
        -p $SLURM_CPUS_PER_TASK \
        2>&1 | tee "${LOG_DIR}/${treatment}/drep_dereplication.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ] && [ -d "${drep_work}/dereplicated_genomes" ]; then
        # Copy dereplicated bins
        cp "${drep_work}/dereplicated_genomes"/*.fa "$output_dir/" 2>/dev/null
        
        local output_bins=$(ls -1 "$output_dir"/*.fa 2>/dev/null | wc -l)
        local reduction=$((input_bins - output_bins))
        
        log "  Dereplication complete: $input_bins → $output_bins bins ($reduction removed)"
        
        # Copy dRep results
        if [ -f "${drep_work}/data_tables/genomeInformation.csv" ]; then
            cp "${drep_work}/data_tables/genomeInformation.csv" "${output_dir}/../drep_genome_info.csv"
        fi
        
        if [ -f "${drep_work}/data_tables/Cdb.csv" ]; then
            cp "${drep_work}/data_tables/Cdb.csv" "${output_dir}/../drep_clusters.csv"
        fi
    else
        log "ERROR: dRep failed, copying all bins without dereplication"
        cp "$bins_dir"/*.fa "$output_dir/" 2>/dev/null
    fi
    
    conda deactivate
    return 0
}

# Function to create treatment-level statistics
create_treatment_statistics() {
    local treatment="$1"
    local samples=($2)
    local output_dir="$3"
    local total_bins="$4"
    local hq_bins="$5"
    local mq_bins="$6"
    local lq_bins="$7"
    
    log "Creating treatment-level statistics for $treatment..."
    
    local stats_file="${output_dir}/treatment_statistics.tsv"
    local summary_file="${output_dir}/treatment_summary.txt"
    
    # Create per-sample statistics
    cat > "$stats_file" << EOF
Sample\tTotal_Bins\tHigh_Quality\tMedium_Quality\tLow_Quality\tHQ_Percentage\tMQ_Percentage\tLQ_Percentage
EOF
    
    local sample_count=0
    for sample in "${samples[@]}"; do
        local sample_collection_dir="${OUTPUT_DIR}/bin_collection/${treatment}/${sample}"
        
        if [ -d "$sample_collection_dir" ]; then
            local s_total=$(ls -1 "${sample_collection_dir}/all_bins"/*.fa 2>/dev/null | wc -l)
            local s_hq=$(ls -1 "${sample_collection_dir}/high_quality"/*.fa 2>/dev/null | wc -l)
            local s_mq=$(ls -1 "${sample_collection_dir}/medium_quality"/*.fa 2>/dev/null | wc -l)
            local s_lq=$(ls -1 "${sample_collection_dir}/low_quality"/*.fa 2>/dev/null | wc -l)
            
            local s_hq_pct=0
            local s_mq_pct=0
            local s_lq_pct=0
            
            if [ $s_total -gt 0 ]; then
                s_hq_pct=$(echo "scale=1; $s_hq * 100 / $s_total" | bc -l)
                s_mq_pct=$(echo "scale=1; $s_mq * 100 / $s_total" | bc -l)
                s_lq_pct=$(echo "scale=1; $s_lq * 100 / $s_total" | bc -l)
            fi
            
            echo -e "${sample}\t${s_total}\t${s_hq}\t${s_mq}\t${s_lq}\t${s_hq_pct}\t${s_mq_pct}\t${s_lq_pct}" >> "$stats_file"
            ((sample_count++))
        fi
    done
    
    # Create summary report
    cat > "$summary_file" << EOF
Treatment-Level Analysis Summary for $treatment
==============================================

Date: $(date)
Treatment: $treatment
Samples processed: $sample_count

Bin Quality Distribution:
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

Per-Sample Statistics:
- Average bins per sample: $(echo "scale=1; $total_bins / $sample_count" | bc -l)
- Average HQ bins per sample: $(echo "scale=1; $hq_bins / $sample_count" | bc -l)
- Average MQ bins per sample: $(echo "scale=1; $mq_bins / $sample_count" | bc -l)

Output Files:
- All samples bins: ${output_dir}/all_samples/
- High-quality bins: ${output_dir}/high_quality/
- Medium-quality bins: ${output_dir}/medium_quality/
- Low-quality bins: ${output_dir}/low_quality/
- Treatment statistics: ${output_dir}/treatment_statistics.tsv
- Combined summary: ${output_dir}/treatment_bins_summary.tsv
- HTML report: ${output_dir}/treatment_report.html

Data Sources:
- Samples with CheckM2 data: $(find "${OUTPUT_DIR}/checkm2/${treatment}" -name "quality_report.tsv" | wc -l)
- Samples with CoverM data: $(find "${OUTPUT_DIR}/coverm/${treatment}" -name "abundance.tsv" | wc -l)
- Samples with bin collection: $(find "${OUTPUT_DIR}/bin_collection/${treatment}" -name "collection_complete.flag" | wc -l)

EOF
    
    log "Treatment statistics created: $stats_file"
    log "Treatment summary created: $summary_file"
}

# Function to create treatment-level HTML report
create_treatment_html_report() {
    local treatment="$1"
    local samples=($2)
    local output_dir="$3"
    local total_bins="$4"
    local hq_bins="$5"
    local mq_bins="$6"
    local lq_bins="$7"
    
    log "Creating treatment-level HTML report for $treatment..."
    
    local html_file="${output_dir}/treatment_report.html"
    local combined_summary="${output_dir}/treatment_bins_summary.tsv"
    
    cat > "$html_file" << EOF
<!DOCTYPE html>
<html>
<head>
    <title>Treatment Analysis Report - $treatment</title>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1600px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 0 15px rgba(0,0,0,0.1);
        }
        h1 { 
            color: #2c3e50;
            text-align: center;
            padding-bottom: 20px;
            border-bottom: 3px solid #3498db;
        }
        .summary-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
            padding: 20px;
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            border-radius: 10px;
        }
        .stat-box {
            text-align: center;
            padding: 20px;
            background-color: white;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            transition: transform 0.2s;
        }
        .stat-box:hover {
            transform: translateY(-2px);
        }
        .stat-number {
            font-size: 2.5em;
            font-weight: bold;
            color: #3498db;
        }
        .stat-label {
            color: #7f8c8d;
            margin-top: 5px;
            font-weight: 500;
        }
        .quality-legend {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
        }
        .quality-legend h3 {
            margin: 0 0 10px 0;
            color: #2c3e50;
        }
        .quality-item {
            display: inline-block;
            margin: 5px 10px;
            padding: 5px 10px;
            border-radius: 15px;
            font-size: 0.9em;
            font-weight: bold;
        }
        .high { background-color: #d4edda; color: #155724; }
        .medium { background-color: #fff3cd; color: #856404; }
        .low { background-color: #f8d7da; color: #721c24; }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin-top: 20px;
            font-size: 0.9em;
        }
        th, td { 
            text-align: left; 
            padding: 10px; 
            border: 1px solid #ddd; 
        }
        th { 
            background-color: #3498db;
            color: white;
            font-weight: bold;
            cursor: pointer;
            position: sticky;
            top: 0;
        }
        th:hover {
            background-color: #2980b9;
        }
        tr:nth-child(even) { background-color: #f9f9f9; }
        tr:hover { background-color: #f1f1f1; }
        .sample-tabs {
            display: flex;
            border-bottom: 2px solid #3498db;
            margin: 20px 0;
        }
        .tab {
            padding: 10px 20px;
            background-color: #ecf0f1;
            border: none;
            cursor: pointer;
            margin-right: 5px;
            border-radius: 5px 5px 0 0;
        }
        .tab.active {
            background-color: #3498db;
            color: white;
        }
        .tab-content {
            display: none;
            padding: 20px;
            border: 1px solid #ddd;
            border-top: none;
        }
        .tab-content.active {
            display: block;
        }
        .filter-controls {
            margin: 20px 0;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
        }
        .filter-controls label {
            margin-right: 15px;
            font-weight: bold;
        }
        .filter-controls select, .filter-controls input {
            margin-right: 10px;
            padding: 5px;
            border-radius: 3px;
            border: 1px solid #ddd;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>Treatment Analysis Report: $treatment</h1>
        
        <div class="summary-stats">
            <div class="stat-box">
                <div class="stat-number">${#samples[@]}</div>
                <div class="stat-label">Samples</div>
            </div>
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
        
        <div class="quality-legend">
            <h3>Quality Categories</h3>
            <span class="quality-item high">High Quality: ≥90% complete, ≤5% contamination</span>
            <span class="quality-item medium">Medium Quality: ≥50% complete, ≤10% contamination</span>
            <span class="quality-item low">Low Quality: All other bins</span>
        </div>
        
        <div class="filter-controls">
            <label for="sampleFilter">Filter by Sample:</label>
            <select id="sampleFilter">
                <option value="">All Samples</option>
EOF
    
    # Add sample options
    for sample in "${samples[@]}"; do
        echo "                <option value=\"$sample\">$sample</option>" >> "$html_file"
    done
    
    cat >> "$html_file" << EOF
            </select>
            
            <label for="qualityFilter">Filter by Quality:</label>
            <select id="qualityFilter">
                <option value="">All Qualities</option>
                <option value="High_Quality">High Quality</option>
                <option value="Medium_Quality">Medium Quality</option>
                <option value="Low_Quality">Low Quality</option>
            </select>
            
            <label for="completenessFilter">Min Completeness:</label>
            <input type="number" id="completenessFilter" min="0" max="100" placeholder="0">
            
            <label for="contaminationFilter">Max Contamination:</label>
            <input type="number" id="contaminationFilter" min="0" max="100" placeholder="100">
        </div>
        
        <table id="binTable">
            <thead>
                <tr>
                    <th onclick="sortTable(0)">Sample</th>
                    <th onclick="sortTable(1)">Bin</th>
                    <th onclick="sortTable(2)">Source</th>
                    <th onclick="sortTable(3)">Completeness (%)</th>
                    <th onclick="sortTable(4)">Contamination (%)</th>
                    <th onclick="sortTable(5)">GC (%)</th>
                    <th onclick="sortTable(6)">Size (bp)</th>
                    <th onclick="sortTable(7)">N50</th>
                    <th onclick="sortTable(8)">Quality</th>
                    <th onclick="sortTable(9)">Rel. Abundance (%)</th>
                    <th onclick="sortTable(10)">Mean Coverage</th>
                </tr>
            </thead>
            <tbody id="binTableBody">
EOF
    
    # Add table rows from combined summary
    if [ -f "$combined_summary" ]; then
        tail -n +2 "$combined_summary" | while IFS=$'\t' read -r sample bin bin_source completeness contamination strain_het gc genome_size n50 taxonomy quality_category rel_abund mean_cov tpm cov_frac; do
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
                    <td>$sample</td>
                    <td>$bin</td>
                    <td>$bin_source</td>
                    <td>$comp_fmt</td>
                    <td>$cont_fmt</td>
                    <td>$gc_fmt</td>
                    <td>$size_fmt</td>
                    <td>$n50_fmt</td>
                    <td><span class="quality-item $quality_class">$quality_category</span></td>
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
        // Table sorting functionality
        function sortTable(n) {
            var table = document.getElementById("binTable");
            var rows = Array.from(table.rows).slice(1);
            var ascending = !table.dataset.sortAsc || table.dataset.sortAsc === 'false';
            
            rows.sort(function(a, b) {
                var x = a.cells[n].innerText.replace(/[^0-9.-]/g, '') || a.cells[n].innerText;
                var y = b.cells[n].innerText.replace(/[^0-9.-]/g, '') || b.cells[n].innerText;
                
                if (!isNaN(x) && !isNaN(y)) {
                    return ascending ? x - y : y - x;
                } else {
                    return ascending ? x.localeCompare(y) : y.localeCompare(x);
                }
            });
            
            table.dataset.sortAsc = ascending;
            
            var tbody = table.querySelector('tbody');
            tbody.innerHTML = '';
            rows.forEach(row => tbody.appendChild(row));
        }
        
        // Filter functionality
        function filterTable() {
            var sampleFilter = document.getElementById('sampleFilter').value;
            var qualityFilter = document.getElementById('qualityFilter').value;
            var completenessFilter = document.getElementById('completenessFilter').value;
            var contaminationFilter = document.getElementById('contaminationFilter').value;
            
            var table = document.getElementById('binTable');
            var rows = table.getElementsByTagName('tr');
            
            for (var i = 1; i < rows.length; i++) {
                var row = rows[i];
                var sample = row.cells[0].textContent;
                var completeness = parseFloat(row.cells[3].textContent);
                var contamination = parseFloat(row.cells[4].textContent);
                var quality = row.cells[8].textContent;
                
                var showRow = true;
                
                if (sampleFilter && sample !== sampleFilter) showRow = false;
                if (qualityFilter && quality !== qualityFilter) showRow = false;
                if (completenessFilter && completeness < parseFloat(completenessFilter)) showRow = false;
                if (contaminationFilter && contamination > parseFloat(contaminationFilter)) showRow = false;
                
                row.style.display = showRow ? '' : 'none';
            }
        }
        
        // Add event listeners
        document.getElementById('sampleFilter').addEventListener('change', filterTable);
        document.getElementById('qualityFilter').addEventListener('change', filterTable);
        document.getElementById('completenessFilter').addEventListener('input', filterTable);
        document.getElementById('contaminationFilter').addEventListener('input', filterTable);
    </script>
</body>
</html>
EOF
    
    log "Treatment HTML report created: $html_file"
}

# Function to create merged treatment FASTA files
create_treatment_merged_fastas() {
    local treatment="$1"
    local output_dir="$2"
    
    log "Creating merged treatment FASTA files for $treatment..."
    
    for quality in high_quality medium_quality low_quality all_samples; do
        local quality_dir="${output_dir}/${quality}"
        local merged_file="${output_dir}/${treatment}_${quality}_merged.fa"
        
        if [ -d "$quality_dir" ]; then
            > "$merged_file"
            local count=0
            
            for bin_file in "$quality_dir"/*.fa; do
                if [ -f "$bin_file" ]; then
                    # Add treatment prefix to sequence headers
                    awk -v treatment="$treatment" -v bin="$(basename "$bin_file" .fa)" '
                        /^>/ { 
                            print ">" treatment "_" bin "_" substr($0, 2)
                            next 
                        } 
                        { print }
                    ' "$bin_file" >> "$merged_file"
                    ((count++))
                fi
            done
            
            log "  Created ${treatment}_${quality}_merged.fa with $count bins"
        fi
    done
}

# Main processing function
stage_treatment_analysis() {
    local treatment="$1"
    
    log "Running treatment-level analysis for $treatment"
    
    local output_dir="${OUTPUT_DIR}/treatment_analysis/${treatment}"
    mkdir -p "$output_dir"
    
    # Check if already processed
    if [ -f "${output_dir}/treatment_analysis_complete.flag" ]; then
        log "Treatment analysis already completed for $treatment"
        return 0
    fi
    
    # Wait for all samples to complete (with timeout)
    wait_for_sample_completion "$treatment"
    
    # Get completed samples
    local completed_samples=($(get_completed_samples "$treatment"))
    local sample_count=${#completed_samples[@]}
    
    if [ $sample_count -eq 0 ]; then
        log "ERROR: No completed samples found for treatment $treatment"
        return 1
    fi
    
    log "Processing $sample_count completed samples for $treatment"
    
    # Collect all bins from samples
    local results=$(collect_treatment_bins "$treatment" "${completed_samples[*]}" "$output_dir")
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to collect treatment bins"
        return 1
    fi
    
    IFS='|' read -r total_bins hq_bins mq_bins lq_bins <<< "$results"
    
    # Optional: Perform bin dereplication on high-quality bins
    if [ $hq_bins -gt 1 ]; then
        log "Performing dereplication on high-quality bins..."
        local derep_dir="${output_dir}/dereplicated_hq"
        mkdir -p "$derep_dir"
        perform_bin_dereplication "$treatment" "${output_dir}/high_quality" "$derep_dir"
    fi
    
    # Create treatment-level statistics
    create_treatment_statistics "$treatment" "${completed_samples[*]}" "$output_dir" "$total_bins" "$hq_bins" "$mq_bins" "$lq_bins"
    
    # Create HTML report
    create_treatment_html_report "$treatment" "${completed_samples[*]}" "$output_dir" "$total_bins" "$hq_bins" "$mq_bins" "$lq_bins"
    
    # Create merged FASTA files
    create_treatment_merged_fastas "$treatment" "$output_dir"
    
    # Create completion flag
    touch "${output_dir}/treatment_analysis_complete.flag"
    
    log "Treatment-level analysis completed for $treatment"
    return 0
}

# Validation function
validate_treatment_analysis() {
    local treatment="$1"
    local output_dir="${OUTPUT_DIR}/treatment_analysis/${treatment}"
    
    # Check completion flag
    if [ ! -f "${output_dir}/treatment_analysis_complete.flag" ]; then
        return 1
    fi
    
    # Check that output directories exist
    for dir in all_samples high_quality medium_quality low_quality; do
        if [ ! -d "${output_dir}/${dir}" ]; then
            return 1
        fi
    done
    
    # Check that summary files exist
    if [ ! -f "${output_dir}/treatment_bins_summary.tsv" ]; then
        return 1
    fi
    
    if [ ! -f "${output_dir}/treatment_summary.txt" ]; then
        return 1
    fi
    
    return 0
}

# Run the treatment analysis stage
if stage_treatment_analysis "$TREATMENT"; then
    # Validate results
    if validate_treatment_analysis "$TREATMENT"; then
        create_treatment_checkpoint "$TREATMENT" "treatment_analysis"
        log "====== Treatment-level analysis completed for $TREATMENT ======"
    else
        log "ERROR: Treatment analysis validation failed for $TREATMENT"
        exit 1
    fi
else
    log "ERROR: Treatment analysis stage failed for $TREATMENT"
    cleanup_temp_dir "$TEMP_DIR"
    exit 1
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"