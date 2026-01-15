#!/bin/bash
#SBATCH --job-name=binspreader
#SBATCH --array=0-99%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --account=ofinkel

# 04c_binspreader.sh - Graph-aware bin refinement using BinSPreader
# Runs in stage 4 (refinement) after initial binning to refine bins using assembly graph structure

# Source configuration and utilities
if [ -n "$PIPELINE_SCRIPT_DIR" ]; then
    source "${PIPELINE_SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

# BinSPreader executable path
BINSPREADER="/sci/backup/ofinkel/moshea/Tools/SPAdes-4.2.0-Linux/bin/binspreader"

# ===== FUNCTION DEFINITIONS =====

# Convert FASTA bins to TSV format (contig_name \t bin_name)
bins_fasta_to_tsv() {
    local bins_dir="$1"
    local output_tsv="$2"

    log "Converting FASTA bins to TSV format..."

    > "$output_tsv"  # Clear output file

    local bin_count=0
    for bin_file in "${bins_dir}"/*.fa "${bins_dir}"/*.fasta; do
        if [ ! -f "$bin_file" ]; then
            continue
        fi

        local bin_name=$(basename "$bin_file" | sed 's/\.\(fa\|fasta\)$//')

        # Extract contig names from FASTA headers
        grep "^>" "$bin_file" | sed 's/^>//' | awk -v bn="$bin_name" '{print $1"\t"bn}' >> "$output_tsv"
        ((bin_count++))
    done

    if [ $bin_count -eq 0 ]; then
        log "ERROR: No bins found in $bins_dir"
        return 1
    fi

    log "Converted $bin_count bins to TSV format"
    return 0
}

# Convert TSV bins back to FASTA format
bins_tsv_to_fasta() {
    local input_tsv="$1"
    local assembly_fasta="$2"
    local output_dir="$3"

    log "Converting TSV bins back to FASTA format..."

    mkdir -p "$output_dir"

    # Create index of contigs in assembly
    awk '/^>/ {if(seq) print header"\n"seq; header=$0; seq=""; next} {seq=seq $0} END {print header"\n"seq}' "$assembly_fasta" > "${TEMP_DIR}/assembly_formatted.fasta"

    # Group contigs by bin
    awk '{print $2"\t"$1}' "$input_tsv" | sort -k1,1 > "${TEMP_DIR}/bins_sorted.tsv"

    # Extract contigs for each bin
    local current_bin=""
    local bin_count=0

    while IFS=$'\t' read -r bin_name contig_name; do
        if [ "$bin_name" != "$current_bin" ]; then
            current_bin="$bin_name"
            ((bin_count++))
        fi

        # Extract this contig from assembly
        awk -v contig="$contig_name" '
            /^>/ {
                if ($1 == ">"contig || $1 == ">"contig" ") {
                    print_it=1
                } else {
                    print_it=0
                }
            }
            print_it {print}
        ' "${TEMP_DIR}/assembly_formatted.fasta" >> "${output_dir}/${bin_name}.fa"
    done < "${TEMP_DIR}/bins_sorted.tsv"

    log "Created $bin_count FASTA bins"
    return 0
}

# Create dataset.yaml for BinSPreader
create_dataset_yaml() {
    local read1="$1"
    local read2="$2"
    local output_yaml="$3"

    log "Creating dataset.yaml for BinSPreader..."

    cat > "$output_yaml" << EOF
[
  {
    orientation: "fr",
    type: "paired-end",
    right reads: [
      "$read2"
    ],
    left reads: [
      "$read1"
    ]
  }
]
EOF

    log "dataset.yaml created"
    return 0
}

# Run BinSPreader
run_binspreader() {
    local assembly_graph="$1"
    local bins_tsv="$2"
    local output_dir="$3"
    local dataset_yaml="$4"

    log "Running BinSPreader with graph-aware refinement..."
    log "Assembly graph: $assembly_graph"
    log "Bins TSV: $bins_tsv"
    log "Output directory: $output_dir"

    # BinSPreader creates a "tmp" directory in the current working directory
    # Change to parent directory so it has write permissions to create output_dir
    local output_parent=$(dirname "$output_dir")
    mkdir -p "$output_parent"

    local original_dir=$(pwd)
    cd "$output_parent" || return 1

    # Run BinSPreader with multiple assignment mode
    if [ -n "$dataset_yaml" ] && [ -f "$dataset_yaml" ]; then
        log "Using dataset.yaml for read mapping information"
        "$BINSPREADER" \
            "$assembly_graph" \
            "$bins_tsv" \
            "$output_dir" \
            -m \
            -t ${SLURM_CPUS_PER_TASK:-50} \
            --dataset "$dataset_yaml" \
            2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_binspreader.log"
    else
        log "Running without dataset.yaml"
        "$BINSPREADER" \
            "$assembly_graph" \
            "$bins_tsv" \
            "$output_dir" \
            -m \
            -t ${SLURM_CPUS_PER_TASK:-50} \
            2>&1 | tee "${LOG_DIR}/${TREATMENT}/${SAMPLE_NAME}_binspreader.log"
    fi

    local exit_code=${PIPESTATUS[0]}

    # Change back to original directory
    cd "$original_dir"

    if [ $exit_code -eq 0 ] && [ -f "${output_dir}/bins.tsv" ]; then
        log "BinSPreader completed successfully"
        return 0
    else
        log "BinSPreader failed with exit code: $exit_code"
        return 1
    fi
}

# Find assembly graph (GFA file)
find_assembly_graph() {
    local assembly_dir="$1"
    local sample_name="$2"

    log "Searching for assembly graph in $assembly_dir..."

    # Look for GFA files in common locations
    for possible_graph in \
        "${assembly_dir}/assembly_graph_with_scaffolds.gfa" \
        "${assembly_dir}/assembly_graph.gfa" \
        "${assembly_dir}/${sample_name}_assembly_graph.gfa" \
        "${assembly_dir}/scaffolds.gfa" \
        "${assembly_dir}/graph.gfa"; do

        if [ -f "$possible_graph" ]; then
            log "Found assembly graph: $possible_graph"
            echo "$possible_graph"
            return 0
        fi
    done

    log "ERROR: No assembly graph (GFA) file found"
    log "  Checked in: $assembly_dir"
    log "  BinSPreader requires assembly graph from SPAdes"
    return 1
}

# ===== MAIN EXECUTION =====

# Initialize
init_conda
TEMP_DIR=$(setup_temp_dir)

# Determine if we're in treatment-level or sample-level mode
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

if [ "${ASSEMBLY_MODE}" = "coassembly" ]; then
    # ===== TREATMENT-LEVEL MODE (coassembly) =====
    log "Running in TREATMENT-LEVEL mode (coassembly)"

    if [ ! -f "$TREATMENTS_FILE" ]; then
        echo "ERROR: Treatments file not found at: $TREATMENTS_FILE"
        exit 1
    fi

    mapfile -t TREATMENTS_ARRAY < "$TREATMENTS_FILE"

    if [ $TASK_ID -ge ${#TREATMENTS_ARRAY[@]} ]; then
        echo "No treatment found for array index $TASK_ID"
        exit 0
    fi

    TREATMENT="${TREATMENTS_ARRAY[$TASK_ID]}"
    export TREATMENT

    log "====== Starting BinSPreader for Treatment: $TREATMENT ======"

    # Set up directories
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}"
    binspreader_dir="${binning_dir}/binspreader"
    assembly_dir="${OUTPUT_DIR}/coassembly/${TREATMENT}"

    mkdir -p "$binspreader_dir"

    # Check if already processed (check for any refined bin directories)
    already_done=false
    if [ -d "${binspreader_dir}/comebin_refined" ] || [ -d "${binspreader_dir}/semibin_refined" ] || [ -d "${binspreader_dir}/metawrap_refined" ]; then
        log "BinSPreader already completed for treatment $TREATMENT, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly graph
    assembly_graph=$(find_assembly_graph "$assembly_dir" "$TREATMENT")
    if [ -z "$assembly_graph" ]; then
        log "ERROR: Cannot run BinSPreader without assembly graph"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Find assembly FASTA
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/assembly.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly FASTA file found"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Create dataset.yaml with first sample's reads (shared across all refinements)
    dataset_yaml="${TEMP_DIR}/dataset.yaml"
    total_samples=$(get_total_samples)
    found_reads=false

    for i in $(seq 0 $((total_samples - 1))); do
        sample_info=$(get_sample_info_by_index $i 2>/dev/null)
        if [ -n "$sample_info" ]; then
            IFS='|' read -r sample_name sample_treatment _ _ <<< "$sample_info"

            if [ "$sample_treatment" = "$TREATMENT" ]; then
                quality_dir="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/${sample_name}"
                read1="${quality_dir}/filtered_1.fastq.gz"
                read2="${quality_dir}/filtered_2.fastq.gz"

                if [ -f "$read1" ] && [ -f "$read2" ]; then
                    create_dataset_yaml "$read1" "$read2" "$dataset_yaml"
                    found_reads=true
                    break
                fi
            fi
        fi
    done

    # Refine ALL available bin sets (COMEBin, SemiBin, MetaWRAP)
    # Each binner gets its own refined output directory
    refined_any=false

    # Check for COMEBin bins
    if [ "$USE_COMEBIN" = "true" ] && [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        comebin_bins="${binning_dir}/comebin/comebin_res_bins"
        comebin_refined="${binspreader_dir}/comebin_refined"

        log "Refining COMEBin bins..."
        log "Input: $comebin_bins"
        log "Output: $comebin_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/comebin_bins.tsv"
        if bins_fasta_to_tsv "$comebin_bins" "$bins_tsv"; then
            # Run BinSPreader on COMEBin bins
            binspreader_output="${TEMP_DIR}/binspreader_comebin"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$comebin_refined"; then
                    log "COMEBin bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert COMEBin refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on COMEBin bins"
            fi
        else
            log "WARNING: Failed to convert COMEBin bins to TSV"
        fi
    fi

    # Check for SemiBin bins
    if [ "$USE_SEMIBIN" = "true" ] && [ -d "${binning_dir}/semibin/output_bins" ]; then
        semibin_bins="${binning_dir}/semibin/output_bins"
        semibin_refined="${binspreader_dir}/semibin_refined"

        log "Refining SemiBin bins..."
        log "Input: $semibin_bins"
        log "Output: $semibin_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/semibin_bins.tsv"
        if bins_fasta_to_tsv "$semibin_bins" "$bins_tsv"; then
            # Run BinSPreader on SemiBin bins
            binspreader_output="${TEMP_DIR}/binspreader_semibin"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$semibin_refined"; then
                    log "SemiBin bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert SemiBin refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on SemiBin bins"
            fi
        else
            log "WARNING: Failed to convert SemiBin bins to TSV"
        fi
    fi

    # Check for MetaWRAP bins (if neither COMEBin nor SemiBin is used)
    if [ "$USE_COMEBIN" != "true" ] && [ "$USE_SEMIBIN" != "true" ] && [ -d "${binning_dir}/metawrap_bins" ]; then
        metawrap_bins="${binning_dir}/metawrap_bins"
        metawrap_refined="${binspreader_dir}/metawrap_refined"

        log "Refining MetaWRAP bins..."
        log "Input: $metawrap_bins"
        log "Output: $metawrap_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/metawrap_bins.tsv"
        if bins_fasta_to_tsv "$metawrap_bins" "$bins_tsv"; then
            # Run BinSPreader on MetaWRAP bins
            binspreader_output="${TEMP_DIR}/binspreader_metawrap"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$metawrap_refined"; then
                    log "MetaWRAP bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert MetaWRAP refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on MetaWRAP bins"
            fi
        else
            log "WARNING: Failed to convert MetaWRAP bins to TSV"
        fi
    fi

    if [ "$refined_any" = false ]; then
        log "ERROR: BinSPreader failed to refine any bin sets"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "====== BinSPreader completed for treatment $TREATMENT ======"

else
    # ===== SAMPLE-LEVEL MODE (individual assembly) =====
    log "Running in SAMPLE-LEVEL mode (individual assembly)"

    SAMPLE_INFO=$(get_sample_info_by_index "$TASK_ID")
    if [ -z "$SAMPLE_INFO" ]; then
        echo "ERROR: No sample found for array index $TASK_ID"
        exit 1
    fi

    IFS='|' read -r SAMPLE_NAME TREATMENT R1_PATH R2_PATH <<< "$SAMPLE_INFO"
    export SAMPLE_NAME TREATMENT

    create_sample_dirs "$SAMPLE_NAME" "$TREATMENT"

    log "====== Starting BinSPreader for $SAMPLE_NAME ($TREATMENT) ======"

    # Set up directories
    binning_dir="${OUTPUT_DIR}/binning/${TREATMENT}/${SAMPLE_NAME}"
    binspreader_dir="${binning_dir}/binspreader"
    assembly_dir="${OUTPUT_DIR}/assembly/${TREATMENT}/${SAMPLE_NAME}"
    quality_dir="${OUTPUT_DIR}/quality_filtering/${TREATMENT}/${SAMPLE_NAME}"

    mkdir -p "$binspreader_dir"

    # Check if already processed (check for any refined bin directories)
    if [ -d "${binspreader_dir}/comebin_refined" ] || [ -d "${binspreader_dir}/semibin_refined" ] || [ -d "${binspreader_dir}/metawrap_refined" ]; then
        log "BinSPreader already completed for $SAMPLE_NAME, skipping..."
        cleanup_temp_dir "$TEMP_DIR"
        exit 0
    fi

    # Find assembly graph
    assembly_graph=$(find_assembly_graph "$assembly_dir" "$SAMPLE_NAME")
    if [ -z "$assembly_graph" ]; then
        log "ERROR: Cannot run BinSPreader without assembly graph"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Find assembly FASTA
    assembly_fasta=""
    for possible_file in \
        "${assembly_dir}/contigs.fasta" \
        "${assembly_dir}/scaffolds.fasta" \
        "${assembly_dir}/final_contigs.fasta" \
        "${assembly_dir}/assembly.fasta" \
        "${assembly_dir}/${SAMPLE_NAME}_contigs.fasta" \
        "${assembly_dir}/${SAMPLE_NAME}_scaffolds.fasta"; do
        if [ -f "$possible_file" ]; then
            assembly_fasta="$possible_file"
            break
        fi
    done

    if [ -z "$assembly_fasta" ]; then
        log "ERROR: No assembly FASTA file found"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    # Create dataset.yaml with quality-filtered reads (shared across all refinements)
    dataset_yaml="${TEMP_DIR}/dataset.yaml"
    read1="${quality_dir}/filtered_1.fastq.gz"
    read2="${quality_dir}/filtered_2.fastq.gz"

    if [ -f "$read1" ] && [ -f "$read2" ]; then
        create_dataset_yaml "$read1" "$read2" "$dataset_yaml"
    else
        log "WARNING: Quality-filtered reads not found, running without dataset.yaml"
        dataset_yaml=""
    fi

    # Refine ALL available bin sets (COMEBin, SemiBin, MetaWRAP)
    # Each binner gets its own refined output directory
    refined_any=false

    # Check for COMEBin bins
    if [ "$USE_COMEBIN" = "true" ] && [ -d "${binning_dir}/comebin/comebin_res_bins" ]; then
        comebin_bins="${binning_dir}/comebin/comebin_res_bins"
        comebin_refined="${binspreader_dir}/comebin_refined"

        log "Refining COMEBin bins..."
        log "Input: $comebin_bins"
        log "Output: $comebin_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/comebin_bins.tsv"
        if bins_fasta_to_tsv "$comebin_bins" "$bins_tsv"; then
            # Run BinSPreader on COMEBin bins
            binspreader_output="${TEMP_DIR}/binspreader_comebin"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$comebin_refined"; then
                    log "COMEBin bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert COMEBin refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on COMEBin bins"
            fi
        else
            log "WARNING: Failed to convert COMEBin bins to TSV"
        fi
    fi

    # Check for SemiBin bins
    if [ "$USE_SEMIBIN" = "true" ] && [ -d "${binning_dir}/semibin/output_bins" ]; then
        semibin_bins="${binning_dir}/semibin/output_bins"
        semibin_refined="${binspreader_dir}/semibin_refined"

        log "Refining SemiBin bins..."
        log "Input: $semibin_bins"
        log "Output: $semibin_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/semibin_bins.tsv"
        if bins_fasta_to_tsv "$semibin_bins" "$bins_tsv"; then
            # Run BinSPreader on SemiBin bins
            binspreader_output="${TEMP_DIR}/binspreader_semibin"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$semibin_refined"; then
                    log "SemiBin bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert SemiBin refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on SemiBin bins"
            fi
        else
            log "WARNING: Failed to convert SemiBin bins to TSV"
        fi
    fi

    # Check for MetaWRAP bins (if neither COMEBin nor SemiBin is used)
    if [ "$USE_COMEBIN" != "true" ] && [ "$USE_SEMIBIN" != "true" ] && [ -d "${binning_dir}/metawrap_50_10_bins" ]; then
        metawrap_bins="${binning_dir}/metawrap_50_10_bins"
        metawrap_refined="${binspreader_dir}/metawrap_refined"

        log "Refining MetaWRAP bins..."
        log "Input: $metawrap_bins"
        log "Output: $metawrap_refined"

        # Convert bins to TSV
        bins_tsv="${TEMP_DIR}/metawrap_bins.tsv"
        if bins_fasta_to_tsv "$metawrap_bins" "$bins_tsv"; then
            # Run BinSPreader on MetaWRAP bins
            binspreader_output="${TEMP_DIR}/binspreader_metawrap"
            if run_binspreader "$assembly_graph" "$bins_tsv" "$binspreader_output" "$dataset_yaml"; then
                # Convert output TSV back to FASTA
                if bins_tsv_to_fasta "${binspreader_output}/bins.tsv" "$assembly_fasta" "$metawrap_refined"; then
                    log "MetaWRAP bins refined successfully"
                    refined_any=true
                else
                    log "WARNING: Failed to convert MetaWRAP refined bins to FASTA"
                fi
            else
                log "WARNING: BinSPreader failed on MetaWRAP bins"
            fi
        else
            log "WARNING: Failed to convert MetaWRAP bins to TSV"
        fi
    fi

    if [ "$refined_any" = false ]; then
        log "ERROR: BinSPreader failed to refine any bin sets"
        cleanup_temp_dir "$TEMP_DIR"
        exit 1
    fi

    log "====== BinSPreader completed for $SAMPLE_NAME ======"
fi

# Cleanup
cleanup_temp_dir "$TEMP_DIR"
