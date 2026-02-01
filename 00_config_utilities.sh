#!/bin/bash
# 00_config_utilities.sh - Configuration and utility functions for the metagenomic pipeline

###############################################################################
# CONFIGURATION
###############################################################################

# Base directories
export INPUT_DIR="${INPUT_DIR:-/path/to/input/fastq/files}"
export OUTPUT_DIR="${OUTPUT_DIR:-/path/to/output/directory}"
export WORK_DIR="${OUTPUT_DIR}/processing_workdir"
export FINAL_BINS_DIR="${OUTPUT_DIR}/high_quality_bins"
export LOG_DIR="${OUTPUT_DIR}/logs"

# Sample sheet configuration
export SAMPLE_SHEET="${SAMPLE_SHEET:-${INPUT_DIR}/sample_sheet.xlsx}"
export SAMPLE_SHEET_TAB="${SAMPLE_SHEET_TAB:-Sheet1}"

# Database paths
export MAGPURIFYDB="/sci/backup/aerez/aerez/moshea/MAGpurify_database"
export VIRSORTER2_DB="/sci/backup/aerez/aerez/moshea/virsorter2_db"
export TRIMMOMATIC_DB="/sci/home/moshea/miniconda3/envs/metagenome_assembly/share/trimmomatic-0.39-2/adapters"

# EukFinder database paths
export EUKFINDER_CENTRIFUGE_DB="${EUKFINDER_CENTRIFUGE_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/Centrifuge_Sept2020}"
export EUKFINDER_PLAST_DB="${EUKFINDER_PLAST_DB:-/sci/backup/ofinkel/moshea/eukfinder_databases/PlastDB.fasta}"
export EUKFINDER_PLAST_ID_MAP="${EUKFINDER_PLAST_ID_MAP:-/sci/backup/ofinkel/moshea/eukfinder_databases/PlastDB_map.txt}"

# Trimmomatic parameters
export TRIMMOMATIC_ADAPTERS="${TRIMMOMATIC_ADAPTERS:-TruSeq3-PE-2.fa}"
export TRIMMOMATIC_LEADING="${TRIMMOMATIC_LEADING:-3}"
export TRIMMOMATIC_TRAILING="${TRIMMOMATIC_TRAILING:-3}"
export TRIMMOMATIC_SLIDINGWINDOW="${TRIMMOMATIC_SLIDINGWINDOW:-4:15}"
export TRIMMOMATIC_MINLEN="${TRIMMOMATIC_MINLEN:-36}"

# Assembly mode: "individual" or "coassembly"
export ASSEMBLY_MODE="${ASSEMBLY_MODE:-individual}"

# Assembly parameters
export ASSEMBLY_THREADS="${ASSEMBLY_THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
export ASSEMBLY_MEMORY="${ASSEMBLY_MEMORY:-250}"  # Memory in GB for individual assembly
export COASSEMBLY_MEMORY="${COASSEMBLY_MEMORY:-500}"  # Memory in GB for coassembly
export BIN_REASSEMBLY_THREADS="${BIN_REASSEMBLY_THREADS:-8}"  # Threads for bin reassembly
export BIN_REASSEMBLY_MEMORY="${BIN_REASSEMBLY_MEMORY:-32}"  # Memory in GB for bin reassembly

# EukFinder parameters
export EUKFINDER_THREADS="${EUKFINDER_THREADS:-48}"  # Number of threads for EukFinder
export EUKFINDER_CHUNKS="${EUKFINDER_CHUNKS:-6}"  # Number of chunks for PLAST processing
export EUKFINDER_TAXONOMY_UPDATE="${EUKFINDER_TAXONOMY_UPDATE:-False}"  # Update taxonomy database (True/False)
export EUKFINDER_EVALUE="${EUKFINDER_EVALUE:-0.01}"  # E-value threshold for PLAST
export EUKFINDER_PID="${EUKFINDER_PID:-60}"  # Percent identity threshold for PLAST
export EUKFINDER_COV="${EUKFINDER_COV:-30}"  # Percent coverage threshold for PLAST
export EUKFINDER_MHLEN="${EUKFINDER_MHLEN:-100}"  # Minimum hit length for Centrifuge

# Pipeline parameters
export MIN_COMPLETENESS="${MIN_COMPLETENESS:-90}"
export MAX_CONTAMINATION="${MAX_CONTAMINATION:-5}"

# Conda setup
export CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"

# SLURM account (for job submission)
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-}"

# File-based sample storage - ONLY initialize if not already set
# Use ${VAR:-} syntax to handle unbound variables when set -u is enabled
if [ -z "${SAMPLE_INFO_FILE:-}" ]; then
    SAMPLE_INFO_FILE=""
fi

if [ -z "${TREATMENTS_FILE:-}" ]; then
    TREATMENTS_FILE=""
fi

if [ -z "${TOTAL_SAMPLES:-}" ]; then
    TOTAL_SAMPLES=0
fi

###############################################################################
# SAMPLE STORAGE FUNCTIONS
###############################################################################

# Initialize sample information storage
init_sample_storage() {
    SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"
    TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
    mkdir -p "${WORK_DIR}"
    > "$SAMPLE_INFO_FILE"  # Clear the file
    > "$TREATMENTS_FILE"   # Clear the file
    TOTAL_SAMPLES=0
}

# Add sample to storage
add_sample_info() {
    local sample_name="$1"
    local treatment="$2" 
    local r1_path="$3"
    local r2_path="$4"
    
    echo "${sample_name}|${treatment}|${r1_path}|${r2_path}" >> "$SAMPLE_INFO_FILE"
    
    # Add treatment if not already present
    if ! grep -q "^${treatment}$" "$TREATMENTS_FILE" 2>/dev/null; then
        echo "$treatment" >> "$TREATMENTS_FILE"
    fi
    
    ((TOTAL_SAMPLES++))
    log "Added sample: $sample_name ($treatment)"
}

# Get sample info by index - ROBUST VERSION
get_sample_info_by_index() {
    local index="$1"
    
    # Ensure SAMPLE_INFO_FILE is set
    if [ -z "$SAMPLE_INFO_FILE" ]; then
        SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"
    fi
    
    if [ ! -f "$SAMPLE_INFO_FILE" ]; then
        echo "ERROR: Sample info file not found: $SAMPLE_INFO_FILE" >&2
        return 1
    fi
    
    # Get the line at index+1 (since sed is 1-based)
    local line_num=$((index + 1))
    local sample_info=$(sed -n "${line_num}p" "$SAMPLE_INFO_FILE" 2>/dev/null)
    
    if [ -z "$sample_info" ]; then
        echo "No sample found for array index $index" >&2
        return 1
    fi
    
    # Only output the actual sample info to stdout
    echo "$sample_info"
    return 0
}

# Get all treatments - CLEAN VERSION
get_treatments() {
    # Ensure TREATMENTS_FILE is set
    if [ -z "$TREATMENTS_FILE" ]; then
        TREATMENTS_FILE="${WORK_DIR}/treatments.txt"
    fi
    
    if [ ! -f "$TREATMENTS_FILE" ]; then
        return 1
    fi
    
    # Output treatments as newline-separated list (not space-separated)
    cat "$TREATMENTS_FILE" 2>/dev/null
}

# Get samples for a treatment
get_samples_for_treatment() {
    local treatment="$1"
    
    # Ensure SAMPLE_INFO_FILE is set
    if [ -z "$SAMPLE_INFO_FILE" ]; then
        SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"
    fi
    
    if [ ! -f "$SAMPLE_INFO_FILE" ]; then
        return 1
    fi
    
    # Get all sample names for this treatment
    grep "|${treatment}|" "$SAMPLE_INFO_FILE" 2>/dev/null | cut -d'|' -f1 | tr '\n' ' '
}

# Get total number of samples - CLEAN VERSION
get_total_samples() {
    # Ensure SAMPLE_INFO_FILE is set
    if [ -z "$SAMPLE_INFO_FILE" ]; then
        SAMPLE_INFO_FILE="${WORK_DIR}/sample_info.txt"
    fi
    
    if [ -f "$SAMPLE_INFO_FILE" ]; then
        wc -l < "$SAMPLE_INFO_FILE" 2>/dev/null || echo "0"
    else
        echo "0"
    fi
}

###############################################################################
# MULTI-LANE DETECTION AND HANDLING
###############################################################################

# Detect and group lanes for the same sample
# Returns: sample_name|treatment|r1_paths|r2_paths (paths are comma-separated)
# Detect and group lanes for the same sample
detect_sample_lanes() {
    local input_dir="$1"
    
    log "Detecting multi-lane samples in: $input_dir"
    
    # Initialize storage
    init_sample_storage
    
    # Find all R1 files
    declare -A sample_lanes_r1
    declare -A sample_lanes_r2
    declare -A sample_treatments
    
    # Pattern matching for common lane formats - FIXED find command
    while IFS= read -r -d '' r1_file; do
        local basename=$(basename "$r1_file")
        
        # Extract sample name (everything before _DKDN, _L, or last _1)
        local sample_name=$(echo "$basename" | sed -E 's/(_DKDN.*|_L[0-9]+.*|_[12]\.fq\.gz$|_[12]\.fastq\.gz$)//')
        
        # Try to find corresponding R2 file
        local r2_file=""
        local r2_patterns=(
            "${r1_file/_1.fq.gz/_2.fq.gz}"
            "${r1_file/_1.fastq.gz/_2.fastq.gz}"
            "${r1_file/_R1/_R2}"
            "${r1_file/_R1./_R2.}"
        )
        
        for pattern in "${r2_patterns[@]}"; do
            if [ -f "$pattern" ]; then
                r2_file="$pattern"
                break
            fi
        done
        
        if [ -z "$r2_file" ]; then
            log "WARNING: No R2 file found for $r1_file"
            continue
        fi
        
        # Try to extract treatment from path or filename
        local treatment="unknown"
        if [[ "$basename" =~ ^([^_]+)_ ]]; then
            treatment="${BASH_REMATCH[1]}"
        fi
        
        # Append to sample's lane list
        if [ -z "${sample_lanes_r1[$sample_name]}" ]; then
            sample_lanes_r1[$sample_name]="$r1_file"
            sample_lanes_r2[$sample_name]="$r2_file"
            sample_treatments[$sample_name]="$treatment"
        else
            sample_lanes_r1[$sample_name]="${sample_lanes_r1[$sample_name]},$r1_file"
            sample_lanes_r2[$sample_name]="${sample_lanes_r2[$sample_name]},$r2_file"
        fi
        
    done < <(find "$input_dir" -type f \( -name "*_1.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1*.fq.gz" \) -print0 | sort -z)
    
    # Add samples to storage
    for sample_name in "${!sample_lanes_r1[@]}"; do
        local r1_paths="${sample_lanes_r1[$sample_name]}"
        local r2_paths="${sample_lanes_r2[$sample_name]}"
        local treatment="${sample_treatments[$sample_name]}"
        
        # Count number of lanes
        local num_lanes=$(echo "$r1_paths" | tr ',' '\n' | wc -l)
        
        add_sample_info "$sample_name" "$treatment" "$r1_paths" "$r2_paths"
        
        if [ $num_lanes -gt 1 ]; then
            log "Sample $sample_name: $num_lanes lanes detected"
        else
            log "Sample $sample_name: single lane"
        fi
    done
    
    TOTAL_SAMPLES=$(get_total_samples)
    log "Total samples detected: $TOTAL_SAMPLES"
    
    return 0
}

# Check if sample has multiple lanes
sample_has_multiple_lanes() {
    local sample_info="$1"
    IFS='|' read -r _ _ r1_paths _ <<< "$sample_info"
    
    # Check if R1 paths contain comma (multiple files)
    if [[ "$r1_paths" == *","* ]]; then
        return 0  # Has multiple lanes
    else
        return 1  # Single lane
    fi
}

# Count lanes for a sample
count_sample_lanes() {
    local sample_info="$1"
    IFS='|' read -r _ _ r1_paths _ <<< "$sample_info"
    
    echo "$r1_paths" | tr ',' '\n' | wc -l
}

###############################################################################
# READ VALIDATION AND REPAIR FUNCTIONS
###############################################################################

# Count reads in FASTQ file
count_reads() {
    local fastq_file="$1"
    
    if [ ! -f "$fastq_file" ]; then
        echo "0"
        return 1
    fi
    
    if [[ "$fastq_file" == *.gz ]]; then
        zcat "$fastq_file" 2>/dev/null | wc -l | awk '{print $1/4}'
    else
        wc -l < "$fastq_file" 2>/dev/null | awk '{print $1/4}'
    fi
}

# Validate paired-end read counts match
validate_read_counts() {
    local r1_file="$1"
    local r2_file="$2"
    local sample_name="${3:-unknown}"
    
    log "Validating read counts for $sample_name..."
    
    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        log "ERROR: Input files not found"
        return 1
    fi
    
    # Count reads in each file
    local r1_count=$(count_reads "$r1_file")
    local r2_count=$(count_reads "$r2_file")
    
    log "  R1 reads: $r1_count"
    log "  R2 reads: $r2_count"
    
    if [ "$r1_count" != "$r2_count" ]; then
        log "ERROR: Read count mismatch! R1=$r1_count, R2=$r2_count"
        log "  Difference: $((r1_count > r2_count ? r1_count - r2_count : r2_count - r1_count)) reads"
        return 1
    fi
    
    if [ "$r1_count" -eq 0 ]; then
        log "ERROR: Both files are empty!"
        return 1
    fi
    
    log "✅ Read counts match: $r1_count reads"
    return 0
}

# Repair paired-end synchronization using BBMap
repair_paired_reads() {
    local r1_file="$1"
    local r2_file="$2"
    local output_r1="$3"
    local output_r2="$4"
    local singletons="$5"
    local sample_name="${6:-unknown}"
    local threads="${SLURM_CPUS_PER_TASK:-4}"
    
    log "Repairing paired-end synchronization for $sample_name..."
    
    # Activate BBMap environment
    activate_env bbmap
    
    # Run repair.sh
    repair.sh \
        in="$r1_file" \
        in2="$r2_file" \
        out="$output_r1" \
        out2="$output_r2" \
        outs="$singletons" \
        threads="$threads" \
        repair=t \
        2>&1 | tee -a "${LOG_DIR}/${TREATMENT:-unknown}/${sample_name}_repair.log"
    
    local exit_code=$?
    
    if [ $exit_code -ne 0 ]; then
        log "ERROR: repair.sh failed with exit code $exit_code"
        conda deactivate
        return 1
    fi
    
    # Validate repaired files
    if validate_read_counts "$output_r1" "$output_r2" "$sample_name"; then
        local paired_count=$(count_reads "$output_r1")
        local singleton_count=$(count_reads "$singletons")
        log "✅ Repair successful: $paired_count paired reads, $singleton_count singletons"
        conda deactivate
        return 0
    else
        log "ERROR: Files still not synchronized after repair"
        conda deactivate
        return 1
    fi
}

###############################################################################
# SAMPLE SHEET PARSING
###############################################################################

# Detect sample sheet format
detect_sample_sheet_format() {
    local sample_sheet="$1"
    local header_line="$2"
    
    # Convert header to lowercase for case-insensitive matching
    local header_lower=$(echo "$header_line" | tr '[:upper:]' '[:lower:]')
    
    # Check for new format (has 'run' and 'group' columns)
    if [[ "$header_lower" =~ (^|[[:space:],])run([[:space:],]|$) ]] && \
       [[ "$header_lower" =~ (^|[[:space:],])group([[:space:],]|$) ]]; then
        echo "multi_run"
        return 0
    fi
    
    # Check for old format (has 'treatment' column)
    if [[ "$header_lower" =~ (^|[[:space:],])treatment([[:space:],]|$) ]]; then
        echo "legacy"
        return 0
    fi
    
    # Default to legacy if can't determine
    echo "legacy"
    return 0
}

# Parse sample sheet - UPDATED for multiple formats
parse_sample_sheet() {
    local sample_sheet="$1"
    
    if [ ! -f "$sample_sheet" ]; then
        log "ERROR: Sample sheet not found: $sample_sheet"
        return 1
    fi
    
    log "Parsing sample sheet: $sample_sheet"
    
    # Initialize storage
    init_sample_storage
    
    # Determine file format and convert if needed
    local working_sheet="$sample_sheet"
    
    # Check if it's Excel or CSV and convert to TSV
    if [[ "$sample_sheet" == *.xlsx ]] || [[ "$sample_sheet" == *.xls ]]; then
        log "Converting Excel file to TSV..."
        python3 -c "
import pandas as pd
import sys
try:
    df = pd.read_excel('$sample_sheet', sheet_name='${SAMPLE_SHEET_TAB}')
    df.to_csv(sys.stdout, sep='\t', index=False)
except Exception as e:
    print(f'Error reading Excel file: {e}', file=sys.stderr)
    sys.exit(1)
" > "${WORK_DIR}/sample_sheet.tsv" 2>/dev/null || {
            log "ERROR: Failed to parse Excel file. Ensure pandas is installed."
            return 1
        }
        working_sheet="${WORK_DIR}/sample_sheet.tsv"
    elif [[ "$sample_sheet" == *.csv ]]; then
        log "Converting CSV file to TSV..."
        python3 -c "
import pandas as pd
import sys
try:
    df = pd.read_csv('$sample_sheet')
    df.to_csv(sys.stdout, sep='\t', index=False)
except Exception as e:
    print(f'Error reading CSV file: {e}', file=sys.stderr)
    sys.exit(1)
" > "${WORK_DIR}/sample_sheet.tsv" 2>/dev/null || {
            log "ERROR: Failed to parse CSV file. Ensure pandas is installed."
            return 1
        }
        working_sheet="${WORK_DIR}/sample_sheet.tsv"
    fi
    
    # Read header line
    local header_line=$(head -n 1 "$working_sheet")
    log "Sample sheet headers: $header_line"
    
    # Detect format
    local format=$(detect_sample_sheet_format "$working_sheet" "$header_line")
    log "Detected sample sheet format: $format"
    
    if [ "$format" = "multi_run" ]; then
        parse_multi_run_format "$working_sheet"
    else
        parse_legacy_format "$working_sheet"
    fi
    
    return $?
}

# Parse legacy format (Sample_Name, Treatment, R1_File, R2_File)
parse_legacy_format() {
    local sample_sheet="$1"
    
    log "Parsing legacy format sample sheet..."
    
    local line_num=0
    
    while IFS=$'\t' read -r line; do
        ((line_num++))
        
        # Skip empty lines
        if [ -z "$line" ]; then
            continue
        fi
        
        # Skip header
        if [ $line_num -eq 1 ]; then
            continue
        fi
        
        # Parse sample line
        local sample_name treatment r1_file r2_file
        IFS=$'\t' read -r sample_name treatment r1_file r2_file <<< "$line"
        
        # Skip if any required field is empty
        if [ -z "$sample_name" ] || [ -z "$treatment" ] || [ -z "$r1_file" ] || [ -z "$r2_file" ]; then
            log "WARNING: Skipping incomplete line $line_num: $line"
            continue
        fi
        
        # Construct full paths
        local r1_path="${INPUT_DIR}/${r1_file}"
        local r2_path="${INPUT_DIR}/${r2_file}"
        
        # Check if files exist
        if [ ! -f "$r1_path" ] || [ ! -f "$r2_path" ]; then
            log "WARNING: Files not found for sample $sample_name:"
            log "  R1: $r1_path"
            log "  R2: $r2_path"
            continue
        fi
        
        # Add sample to storage
        add_sample_info "$sample_name" "$treatment" "$r1_path" "$r2_path"
        
    done < "$sample_sheet"
    
    TOTAL_SAMPLES=$(get_total_samples)
    log "Successfully parsed $TOTAL_SAMPLES samples from legacy format sample sheet"
    log "Treatments found: $(get_treatments)"
    
    return 0
}

# Parse multi-run format (sample, run, group, short_reads_1, short_reads_2, ...)
parse_multi_run_format() {
    local sample_sheet="$1"
    
    log "Parsing multi-run format sample sheet..."
    
    # Use Python to parse and group by sample
    python3 << PYTHON_SCRIPT > "${WORK_DIR}/grouped_samples.tsv"
import pandas as pd
import sys
import os

try:
    # Read the TSV file
    sample_sheet = sys.argv[1] if len(sys.argv) > 1 else "${WORK_DIR}/sample_sheet.tsv"
    df = pd.read_csv(sample_sheet, sep='\t')
    
    # Standardize column names (case-insensitive)
    df.columns = df.columns.str.lower()
    
    # Required columns mapping
    required = {
        'sample': 'sample',
        'group': 'treatment',
        'short_reads_1': 'r1',
        'short_reads_2': 'r2'
    }
    
    # Check for required columns
    for col in required.keys():
        if col not in df.columns:
            print(f"ERROR: Required column '{col}' not found in sample sheet", file=sys.stderr)
            sys.exit(1)
    
    # Group by sample name and group/treatment
    grouped = df.groupby(['sample', 'group']).agg({
        'short_reads_1': lambda x: ','.join(x.astype(str)),
        'short_reads_2': lambda x: ','.join(x.astype(str))
    }).reset_index()
    
    # Rename columns to match expected format
    grouped.columns = ['sample_name', 'treatment', 'r1_files', 'r2_files']
    
    # Output as TSV
    grouped.to_csv(sys.stdout, sep='\t', index=False)
    
except Exception as e:
    print(f"ERROR parsing multi-run sample sheet: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON_SCRIPT
    
    if [ $? -ne 0 ]; then
        log "ERROR: Failed to parse multi-run sample sheet"
        return 1
    fi
    
    # Now read the grouped samples
    local line_num=0
    
    while IFS=$'\t' read -r sample_name treatment r1_files r2_files; do
        ((line_num++))
        
        # Skip header
        if [ $line_num -eq 1 ]; then
            continue
        fi
        
        # Skip empty lines
        if [ -z "$sample_name" ]; then
            continue
        fi
        
        # Process comma-separated file lists
        IFS=',' read -ra R1_ARRAY <<< "$r1_files"
        IFS=',' read -ra R2_ARRAY <<< "$r2_files"
        
        # Validate all files exist and build full paths
        local r1_paths=""
        local r2_paths=""
        local all_found=true
        
        for i in "${!R1_ARRAY[@]}"; do
            local r1_file="${R1_ARRAY[$i]}"
            local r2_file="${R2_ARRAY[$i]}"
            
            # Check if path is absolute (starts with /)
if [[ "$r1_file" == /* ]]; then
    r1_path="$r1_file"  # Use as-is for absolute paths
else
    r1_path="${INPUT_DIR}/${r1_file}"  # Prepend INPUT_DIR for relative paths
fi
if [[ "$r2_file" == /* ]]; then
    r2_path="$r2_file"  # Use as-is for absolute paths
else
    r2_path="${INPUT_DIR}/${r2_file}"  # Prepend INPUT_DIR for relative paths
fi
            
            # Check if files exist
            if [ ! -f "$r1_path" ] || [ ! -f "$r2_path" ]; then
                log "WARNING: Files not found for run $((i+1)) of sample $sample_name:"
                log "  R1: $r1_path"
                log "  R2: $r2_path"
                all_found=false
                break
            fi
            
            # Build comma-separated path lists
            if [ -z "$r1_paths" ]; then
                r1_paths="$r1_path"
                r2_paths="$r2_path"
            else
                r1_paths="${r1_paths},${r1_path}"
                r2_paths="${r2_paths},${r2_path}"
            fi
        done
        
        if [ "$all_found" = false ]; then
            log "WARNING: Skipping sample $sample_name due to missing files"
            continue
        fi
        
        # Count runs for this sample
        local num_runs=${#R1_ARRAY[@]}
        
        # Add sample to storage
        add_sample_info "$sample_name" "$treatment" "$r1_paths" "$r2_paths"
        
        if [ $num_runs -gt 1 ]; then
            log "Sample $sample_name: $num_runs runs detected (will be merged in stage -1)"
        else
            log "Sample $sample_name: single run"
        fi
        
    done < "${WORK_DIR}/grouped_samples.tsv"
    
    TOTAL_SAMPLES=$(get_total_samples)
    log "Successfully parsed $TOTAL_SAMPLES samples from multi-run format sample sheet"
    log "Treatments/groups found: $(get_treatments)"
    
    return 0
}

# Auto-discover samples - UPDATED for file storage
auto_discover_samples() {
    local input_dir="$1"
    local pattern="${2:-*_R1*.fastq.gz}"
    
    log "Auto-discovering samples in: $input_dir"
    log "Pattern: $pattern"
    
    # Initialize storage
    init_sample_storage
    local default_treatment="unknown"
    
    # Find all R1 files
    while IFS= read -r -d '' r1_file; do
        local r1_path="$r1_file"
        
        # Try to find corresponding R2 file
        local r2_path=""
        local base_name=$(basename "$r1_file")
        
        # Common R1/R2 patterns
        local r2_patterns=(
            "${r1_file/_R1/_R2}"
            "${r1_file/_R1./_R2.}"
            "${r1_file/_1./_2.}"
            "${r1_file/_1_/_2_}"
            "${r1_file/.R1./.R2.}"
        )
        
        for pattern in "${r2_patterns[@]}"; do
            if [ -f "$pattern" ]; then
                r2_path="$pattern"
                break
            fi
        done
        
        if [ -z "$r2_path" ]; then
            log "WARNING: No R2 file found for $r1_file"
            continue
        fi
        
        # Extract sample name from filename
        local sample_name=$(basename "$r1_file" | sed -E 's/_R1.*//; s/_1.*//; s/\.R1.*//')
        
        # Try to extract treatment from path or filename
        local treatment="$default_treatment"
        if [[ "$r1_file" =~ /([^/]+)/ ]]; then
            # Use parent directory as treatment if it looks like a treatment name
            local parent_dir=$(dirname "$r1_file" | xargs basename)
            if [[ "$parent_dir" != "." && "$parent_dir" != "$input_dir" ]]; then
                treatment="$parent_dir"
            fi
        fi
        
        # Add sample to storage
        add_sample_info "$sample_name" "$treatment" "$r1_path" "$r2_path"
        
    done < <(find "$input_dir" -name "$pattern" -type f -print0)
    
    TOTAL_SAMPLES=$(get_total_samples)
    log "Auto-discovered $TOTAL_SAMPLES samples"
    log "Treatments found: $(get_treatments)"
    
    return 0
}

# Initialize sample information - UPDATED for multi-lane
init_sample_info() {
    # Create essential directories first
    mkdir -p "${WORK_DIR}"
    mkdir -p "${LOG_DIR}"
    mkdir -p "${LOG_DIR}/slurm"  # Create slurm logs subdirectory
    mkdir -p "${OUTPUT_DIR}/checkpoints"

    if [ -f "$SAMPLE_SHEET" ]; then
        log "Using sample sheet: $SAMPLE_SHEET"
        parse_sample_sheet "$SAMPLE_SHEET"
    else
        log "No sample sheet found, auto-discovering samples with lane detection..."
        detect_sample_lanes "$INPUT_DIR"
    fi
    
    TOTAL_SAMPLES=$(get_total_samples)
    if [ $TOTAL_SAMPLES -eq 0 ]; then
        log "ERROR: No samples found!"
        return 1
    fi
    
    log "Final sample count = $TOTAL_SAMPLES"
    
    # Report which samples have multiple lanes
    log "Checking for multi-lane samples..."
    local multi_lane_count=0
    for i in $(seq 0 $((TOTAL_SAMPLES - 1))); do
        local sample_info=$(get_sample_info_by_index $i 2>/dev/null)
        if [ -n "$sample_info" ]; then
            if sample_has_multiple_lanes "$sample_info"; then
                IFS='|' read -r sample_name treatment _ _ <<< "$sample_info"
                local num_lanes=$(count_sample_lanes "$sample_info")
                log "  $sample_name: $num_lanes lanes"
                ((multi_lane_count++))
            fi
        fi
    done
    
    if [ $multi_lane_count -gt 0 ]; then
        log "Found $multi_lane_count samples with multiple lanes"
        log "Lane merging will be performed in stage -1"
    else
        log "All samples have single lanes"
    fi
    
    return 0
}

###############################################################################
# UTILITY FUNCTIONS
###############################################################################

# Initialize conda
init_conda() {
    if command -v module >/dev/null 2>&1; then
        module purge
    fi
    source ${CONDA_BASE}/etc/profile.d/conda.sh
}

# Logging function - FIXED to avoid stdout pollution
log() {
    local sample_name="${SAMPLE_NAME:-GENERAL}"
    local treatment="${TREATMENT:-GENERAL}"
    local logfile="${LOG_DIR}/${treatment}/${sample_name}_pipeline.log"
    mkdir -p $(dirname "$logfile") 2>/dev/null
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$logfile" >&2
}

# Check if checkpoint exists for a sample
check_sample_checkpoint() {
    local sample_name="$1"
    local stage="$2"
    local treatment="${TREATMENT:-unknown}"
    
    [ -f "${OUTPUT_DIR}/checkpoints/${treatment}/${sample_name}_${stage}_complete" ]
}

# Create checkpoint for a sample
create_sample_checkpoint() {
    local sample_name="$1"
    local stage="$2"
    local treatment="${TREATMENT:-unknown}"
    
    mkdir -p "${OUTPUT_DIR}/checkpoints/${treatment}"
    touch "${OUTPUT_DIR}/checkpoints/${treatment}/${sample_name}_${stage}_complete"
    
    log "Checkpoint created: ${sample_name}_${stage}"
}

# Check if treatment-level checkpoint exists
check_treatment_checkpoint() {
    local treatment="$1"
    local stage="$2"
    
    [ -f "${OUTPUT_DIR}/checkpoints/${treatment}/${stage}_complete" ]
}

# Create treatment-level checkpoint
create_treatment_checkpoint() {
    local treatment="$1"
    local stage="$2"
    
    mkdir -p "${OUTPUT_DIR}/checkpoints/${treatment}"
    touch "${OUTPUT_DIR}/checkpoints/${treatment}/${stage}_complete"
    
    log "Treatment checkpoint created: ${treatment}_${stage}"
}

# Legacy checkpoint functions for backward compatibility
check_checkpoint() {
    local stage="$1"
    local treatment="${TREATMENT:-unknown}"
    check_treatment_checkpoint "$treatment" "$stage"
}

create_checkpoint() {
    local stage="$1"
    local treatment="${TREATMENT:-unknown}"
    create_treatment_checkpoint "$treatment" "$stage"
}

# Activate conda environment
activate_env() {
    local env_name="$1"
    export CONDA_BASE="${CONDA_BASE:-/sci/home/moshea/miniconda3}"
    if [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        log "ERROR: conda.sh not found at $CONDA_BASE"
        exit 1
    fi
    # Temporarily disable nounset (set -u) as conda activation scripts have unbound variables
    local nounset_was_set=false
    if [[ $- == *u* ]]; then
        nounset_was_set=true
        set +u
    fi
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "$env_name" || {
        # Restore nounset before exiting
        if [ "$nounset_was_set" = true ]; then
            set -u
        fi
        log "ERROR: Failed to activate $env_name"
        exit 1
    }
    # Restore nounset if it was previously set
    if [ "$nounset_was_set" = true ]; then
        set -u
    fi
    log "Activated conda env: $env_name"
}

# Deactivate conda environment (handles set -u compatibility)
deactivate_env() {
    # Temporarily disable nounset (set -u) as conda deactivation scripts have unbound variables
    local nounset_was_set=false
    if [[ $- == *u* ]]; then
        nounset_was_set=true
        set +u
    fi
    conda deactivate 2>/dev/null || true
    # Restore nounset if it was previously set
    if [ "$nounset_was_set" = true ]; then
        set -u
    fi
}

# Create directories for a sample
create_sample_dirs() {
    local sample_name="$1"
    local treatment="$2"
    local dirs=(
        "${LOG_DIR}/${treatment}"
        "${WORK_DIR}"
        "${OUTPUT_DIR}/checkpoints/${treatment}"
        "${OUTPUT_DIR}/merged_lanes/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/validated/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/plasmids/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/binning/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/bin_refinement/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/reassembly/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/magpurify/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/checkm2/${treatment}/${sample_name}"
        "${OUTPUT_DIR}/coverm/${treatment}/${sample_name}"
        "${FINAL_BINS_DIR}/${treatment}"
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
    done
}

# Create all necessary directories for a treatment
create_treatment_dirs() {
    local treatment="$1"
    local dirs=(
        "${LOG_DIR}/${treatment}"
        "${WORK_DIR}"
        "${OUTPUT_DIR}/checkpoints/${treatment}"
        "${OUTPUT_DIR}/merged_lanes/${treatment}"
        "${OUTPUT_DIR}/validated/${treatment}"
        "${OUTPUT_DIR}/quality_filtering/${treatment}"
        "${OUTPUT_DIR}/assembly/${treatment}"
        "${OUTPUT_DIR}/coassembly/${treatment}"
        "${OUTPUT_DIR}/plasmids/${treatment}"
        "${OUTPUT_DIR}/binning/${treatment}"
        "${OUTPUT_DIR}/bin_refinement/${treatment}"
        "${OUTPUT_DIR}/reassembly/${treatment}"
        "${OUTPUT_DIR}/magpurify/${treatment}"
        "${OUTPUT_DIR}/checkm2/${treatment}"
        "${OUTPUT_DIR}/coverm/${treatment}"
        "${FINAL_BINS_DIR}/${treatment}"
    )
    
    for dir in "${dirs[@]}"; do
        mkdir -p "$dir"
    done
}

# Get sample directories for a treatment (updated for new structure)
get_sample_dirs() {
    local base_dir="$1"
    find "$base_dir" -mindepth 1 -maxdepth 1 -type d
}

# Setup temporary directory
setup_temp_dir() {
    local sample_name="${SAMPLE_NAME:-UNKNOWN}"
    local treatment="${TREATMENT:-UNKNOWN}"
    local job_id="${SLURM_JOB_ID:-$$}"
    local temp_dir="${WORK_DIR}/temp_${treatment}_${sample_name}_${job_id}"
    mkdir -p "$temp_dir"
    export TMPDIR="$temp_dir"
    echo "$temp_dir"
}

# Cleanup temporary directory
cleanup_temp_dir() {
    local temp_dir="$1"
    if [ -n "$temp_dir" ] && [ -d "$temp_dir" ]; then
        rm -rf "$temp_dir"
    fi
}

# Create sample sheet template
create_sample_sheet_template() {
    local output_file="${1:-${INPUT_DIR}/sample_sheet_template.tsv}"
    
    cat > "$output_file" << EOF
Sample_Name	Treatment	R1_File	R2_File
sample1	treatment1	sample1_R1.fastq.gz	sample1_R2.fastq.gz
sample2	treatment1	sample2_R1.fastq.gz	sample2_R2.fastq.gz
sample3	treatment2	sample3_R1.fastq.gz	sample3_R2.fastq.gz
EOF
    
    log "Sample sheet template created: $output_file"
    log "Columns: Sample_Name, Treatment, R1_File, R2_File"
    log "R1_File and R2_File should be relative to INPUT_DIR"
}

###############################################################################
# VALIDATION FUNCTIONS
###############################################################################

# Validate quality filtering results
validate_quality_filtering() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/quality_filtering/${treatment}/${sample_name}"
    
    if [ -f "${output_dir}/filtered_1.fastq.gz" ] && [ -f "${output_dir}/filtered_2.fastq.gz" ]; then
        # Check files are not empty
        if [ -s "${output_dir}/filtered_1.fastq.gz" ] && [ -s "${output_dir}/filtered_2.fastq.gz" ]; then
            return 0
        fi
    fi
    
    return 1
}

# Validate assembly results
validate_assembly() {
    local sample_name="$1"
    local treatment="$2"
    local output_dir="${OUTPUT_DIR}/assembly/${treatment}/${sample_name}"
    
    if [ -f "${output_dir}/contigs.fasta" ] && [ -s "${output_dir}/contigs.fasta" ]; then
        return 0
    fi
    
    return 1
}

###############################################################################
# ASSEMBLY QUALITY FUNCTIONS
###############################################################################

# Calculate assembly success rate (percentage of reads assembled)
calculate_assembly_success_rate() {
    local contigs_file="$1"
    local r1_file="$2"
    local r2_file="$3"
    local sample_name="${4:-unknown}"
    local output_prefix="${5:-${TEMP_DIR}/assembly_mapping}"
    local threads="${SLURM_CPUS_PER_TASK:-8}"

    log "Calculating assembly success rate for $sample_name..."

    # Validate inputs
    if [ ! -f "$contigs_file" ] || [ ! -s "$contigs_file" ]; then
        log "ERROR: Contigs file not found or empty: $contigs_file"
        echo "0.00"
        return 1
    fi

    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        log "ERROR: Input read files not found"
        echo "0.00"
        return 1
    fi

    # Count total input reads
    local total_reads=$(count_reads "$r1_file")
    log "  Total input reads: $total_reads"

    if [ "$total_reads" -eq 0 ]; then
        log "ERROR: No reads found in input files"
        echo "0.00"
        return 1
    fi

    # Activate BBMap environment
    activate_env bbmap

    # Create output directory
    mkdir -p "$(dirname "$output_prefix")"

    # Map reads back to contigs using BBMap
    local stats_file="${output_prefix}_stats.txt"

    # Create log directory for mapping output
    local mapping_log_dir="$(dirname "$output_prefix")"
    mkdir -p "$mapping_log_dir"
    local mapping_log="${output_prefix}_mapping.log"

    log "  Mapping reads back to contigs..."
    bbmap.sh \
        in="$r1_file" \
        in2="$r2_file" \
        ref="$contigs_file" \
        nodisk \
        threads="$threads" \
        statsfile="$stats_file" \
        scafstats="${output_prefix}_scafstats.txt" \
        covstats="${output_prefix}_covstats.txt" \
        rpkm="${output_prefix}_rpkm.txt" \
        sortscafs=f \
        nzo=f \
        slow=f \
        fast=t \
        maxindel=100 \
        ambiguous=random \
        2>&1 | tee "$mapping_log" | grep -E "^(Percent )?mapped:" >&2

    local exit_code=${PIPESTATUS[0]}

    if [ $exit_code -ne 0 ]; then
        log "ERROR: BBMap failed with exit code $exit_code"
        conda deactivate
        echo "0.00"
        return 1
    fi

    # Parse mapping statistics
    local mapped_reads=0
    local percent_mapped=0.00

    # Extract mapped read count from stats file
    if [ -f "$stats_file" ]; then
        # BBMap reports "Percent mapped:" with the percentage
        # Get the main percentage value (not "properly paired" or other sub-categories)
        percent_mapped=$(grep "^Percent mapped:" "$stats_file" | head -1 | awk '{print $NF}' | sed 's/%//' || echo "0.00")

        if [ -z "$percent_mapped" ] || [ "$percent_mapped" = "0.00" ]; then
            # Try alternative parsing - look for "mapped:" line
            percent_mapped=$(grep "^mapped:" "$stats_file" | head -1 | awk '{print $2}' | sed 's/%//' || echo "0.00")
        fi

        # Clean up any potential newlines or extra characters
        percent_mapped=$(echo "$percent_mapped" | tr -d '\n\r' | awk '{print $1}')

        # Calculate mapped reads from percentage
        mapped_reads=$(awk -v total="$total_reads" -v percent="$percent_mapped" 'BEGIN {printf "%.0f", total * percent / 100}')

        log "  Mapped reads: $mapped_reads ($percent_mapped%)"
        log "  Assembly success rate: ${percent_mapped}%"
    else
        log "WARNING: Stats file not found, cannot calculate mapping rate"
        percent_mapped="0.00"
    fi

    conda deactivate

    # Return the percentage
    echo "$percent_mapped"
    return 0
}

# Calculate assembly success rate for coassembly (with merged reads)
calculate_coassembly_success_rate() {
    local contigs_file="$1"
    local merged_r1="$2"
    local merged_r2="$3"
    local treatment="${4:-unknown}"
    local output_prefix="${5:-${TEMP_DIR}/coassembly_mapping}"

    calculate_assembly_success_rate "$contigs_file" "$merged_r1" "$merged_r2" "$treatment" "$output_prefix"
}

# Export all functions
export -f log check_sample_checkpoint create_sample_checkpoint check_treatment_checkpoint create_treatment_checkpoint
export -f check_checkpoint create_checkpoint activate_env create_sample_dirs create_treatment_dirs
export -f get_sample_dirs setup_temp_dir cleanup_temp_dir init_conda
export -f parse_sample_sheet get_sample_info_by_index get_treatments get_samples_for_treatment
export -f auto_discover_samples get_total_samples init_sample_info create_sample_sheet_template
export -f validate_quality_filtering validate_assembly init_sample_storage add_sample_info
export -f detect_sample_lanes sample_has_multiple_lanes count_sample_lanes
export -f count_reads validate_read_counts repair_paired_reads
export -f calculate_assembly_success_rate calculate_coassembly_success_rate