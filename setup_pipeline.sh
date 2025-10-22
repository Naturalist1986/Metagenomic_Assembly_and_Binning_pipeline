#!/bin/bash

# setup_pipeline.sh - Setup script for the metagenomic pipeline

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Usage function
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Setup script for the metagenomic pipeline.

OPTIONS:
    -c, --conda-base PATH    Path to conda/miniconda installation [default: $HOME/miniconda3]
    -i, --input-dir PATH     Input directory containing FASTQ files
    -o, --output-dir PATH    Output directory for results
    -d, --databases PATH     Directory for databases
    -e, --envs-only         Only create conda environments
    -t, --test              Run pipeline test
    -h, --help              Show this help message

EXAMPLES:
    # Full setup
    $0 -i /path/to/fastq -o /path/to/output -d /path/to/databases

    # Only create conda environments
    $0 --envs-only

    # Test installation
    $0 --test -i /path/to/fastq -o /path/to/output

EOF
}

# Parse arguments
ENVS_ONLY=false
RUN_TEST=false
INPUT_DIR_ARG=""
OUTPUT_DIR_ARG=""
DATABASES_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--conda-base)
            CONDA_BASE="$2"
            shift 2
            ;;
        -i|--input-dir)
            INPUT_DIR_ARG="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR_ARG="$2"
            shift 2
            ;;
        -d|--databases)
            DATABASES_DIR="$2"
            shift 2
            ;;
        -e|--envs-only)
            ENVS_ONLY=true
            shift
            ;;
        -t|--test)
            RUN_TEST=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Check conda installation
check_conda() {
    log_info "Checking conda installation..."
    
    if [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        log_error "Conda not found at: $CONDA_BASE"
        log_error "Please install conda/miniconda or specify correct path with -c"
        exit 1
    fi
    
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    log_info "Conda found at: $CONDA_BASE"
}

# Create conda environment
create_conda_env() {
    local env_name="$1"
    local env_file="$2"
    
    log_info "Creating conda environment: $env_name"
    
    if conda env list | grep -q "^$env_name "; then
        log_warn "Environment $env_name already exists"
        read -p "Do you want to recreate it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            conda env remove -n "$env_name" -y
        else
            log_info "Skipping $env_name"
            return 0
        fi
    fi
    
    if [ -f "$env_file" ]; then
        conda env create -f "$env_file"
    else
        # Create environment with basic packages
        case "$env_name" in
            "trimmomatic")
                conda create -n "$env_name" -c bioconda trimmomatic -y
                ;;
            "bbmap")
                conda create -n "$env_name" -c bioconda bbmap -y
                ;;
            "spades")
                conda create -n "$env_name" -c bioconda spades -y
                ;;
            "scapp")
                conda create -n "$env_name" -c bioconda plasclass mob_suite -y
                ;;
            "mobsuite")
                conda create -n "$env_name" -c bioconda mob_suite -y
                ;;
            "metawrap-env")
                conda create -n "$env_name" -c bioconda metawrap-mg -y
                ;;
            "metagenome_assembly")
                conda create -n "$env_name" -c bioconda magpurify -y
                ;;
            "checkm")
                conda create -n "$env_name" -c bioconda checkm2 coverm -y
                ;;
            *)
                log_error "Unknown environment: $env_name"
                return 1
                ;;
        esac
    fi
    
    log_info "Environment $env_name created successfully"
}

# Setup conda environments
setup_environments() {
    log_info "Setting up conda environments..."
    
    # List of required environments
    local envs=(
        "trimmomatic"
        "bbmap"
        "spades"
        "scapp"
        "mobsuite"
        "metawrap-env"
        "metagenome_assembly"
        "checkm"
    )
    
    for env in "${envs[@]}"; do
        create_conda_env "$env" "${SCRIPT_DIR}/envs/${env}.yml"
    done
    
    log_info "All environments created successfully"
}

# Setup databases
setup_databases() {
    if [ -z "$DATABASES_DIR" ]; then
        log_warn "No databases directory specified, skipping database setup"
        return 0
    fi
    
    log_info "Setting up databases in: $DATABASES_DIR"
    mkdir -p "$DATABASES_DIR"
    
    # Setup Trimmomatic adapters
    log_info "Setting up Trimmomatic adapters..."
    local trimmomatic_dir="${DATABASES_DIR}/trimmomatic"
    mkdir -p "$trimmomatic_dir"
    
    # Copy adapters from conda environment
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate trimmomatic
    
    local conda_prefix=$(conda info --base)/envs/trimmomatic
    if [ -d "${conda_prefix}/share/trimmomatic/adapters" ]; then
        cp -r "${conda_prefix}/share/trimmomatic/adapters"/* "$trimmomatic_dir/"
        log_info "Trimmomatic adapters copied to: $trimmomatic_dir"
    else
        log_warn "Trimmomatic adapters not found in conda environment"
    fi
    
    conda deactivate
    
    # Setup CheckM2 database
    log_info "Setting up CheckM2 database..."
    conda activate checkm
    
    if command -v checkm2 &> /dev/null; then
        checkm2 database --download --path "${DATABASES_DIR}/checkm2"
        log_info "CheckM2 database downloaded to: ${DATABASES_DIR}/checkm2"
    else
        log_warn "CheckM2 not available, skipping database setup"
    fi
    
    conda deactivate
    
    # Setup other databases (MAGpurify, etc.)
    log_info "Other databases need to be set up manually:"
    log_info "  - MAGpurify database: https://github.com/snayfach/MAGpurify"
    log_info "  - VirSorter2 database: https://github.com/jiarong/VirSorter2"
}

# Update configuration file
update_config() {
    local config_file="${SCRIPT_DIR}/00_config_utilities.sh"
    local config_backup="${config_file}.backup.$(date +%Y%m%d_%H%M%S)"
    
    log_info "Updating configuration file..."
    
    # Create backup
    cp "$config_file" "$config_backup"
    log_info "Configuration backup created: $config_backup"
    
    # Update configuration
    if [ -n "$INPUT_DIR_ARG" ]; then
        sed -i "s|^export INPUT_DIR=.*|export INPUT_DIR=\"$INPUT_DIR_ARG\"|" "$config_file"
        log_info "Updated INPUT_DIR: $INPUT_DIR_ARG"
    fi
    
    if [ -n "$OUTPUT_DIR_ARG" ]; then
        sed -i "s|^export OUTPUT_DIR=.*|export OUTPUT_DIR=\"$OUTPUT_DIR_ARG\"|" "$config_file"
        log_info "Updated OUTPUT_DIR: $OUTPUT_DIR_ARG"
    fi
    
    if [ -n "$DATABASES_DIR" ]; then
        sed -i "s|^export MAGPURIFYDB=.*|export MAGPURIFYDB=\"$DATABASES_DIR/magpurify\"|" "$config_file"
        sed -i "s|^export VIRSORTER2_DB=.*|export VIRSORTER2_DB=\"$DATABASES_DIR/virsorter2\"|" "$config_file"
        sed -i "s|^export TRIMMOMATIC_DB=.*|export TRIMMOMATIC_DB=\"$DATABASES_DIR/trimmomatic\"|" "$config_file"
        log_info "Updated database paths"
    fi
    
    # Update conda base
    sed -i "s|^export CONDA_BASE=.*|export CONDA_BASE=\"$CONDA_BASE\"|" "$config_file"
    log_info "Updated CONDA_BASE: $CONDA_BASE"
}

# Test pipeline
test_pipeline() {
    log_info "Testing pipeline installation..."
    
    # Check required scripts
    local required_scripts=(
        "00_quality_filtering.sh"
        "01_assembly.sh"
        "run_pipeline.sh"
    )
    
    for script in "${required_scripts[@]}"; do
        if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
            log_error "Required script not found: $script"
            exit 1
        fi
    done
    
    # Test conda environments
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    
    local test_envs=("trimmomatic" "spades" "checkm")
    for env in "${test_envs[@]}"; do
        if ! conda env list | grep -q "^$env "; then
            log_error "Environment not found: $env"
            exit 1
        fi
    done
    
    # Test Python dependencies
    if ! python3 -c "import pandas" 2>/dev/null; then
        log_warn "pandas not found - Excel support may not work"
    fi
    
    log_info "Pipeline test completed successfully"
}

# Create sample dataset
create_sample_dataset() {
    if [ -z "$INPUT_DIR_ARG" ]; then
        log_warn "No input directory specified, skipping sample dataset creation"
        return 0
    fi
    
    log_info "Creating sample dataset structure..."
    
    mkdir -p "$INPUT_DIR_ARG"
    
    # Create sample sheet template
    cat > "${INPUT_DIR_ARG}/sample_sheet.tsv" << EOF
Sample_Name	Treatment	R1_File	R2_File
sample1	treatment1	sample1_R1.fastq.gz	sample1_R2.fastq.gz
sample2	treatment1	sample2_R1.fastq.gz	sample2_R2.fastq.gz
sample3	treatment2	sample3_R1.fastq.gz	sample3_R2.fastq.gz
EOF
    
    log_info "Sample sheet template created: ${INPUT_DIR_ARG}/sample_sheet.tsv"
    log_info "Add your actual FASTQ files to: $INPUT_DIR_ARG"
}

# Main execution
main() {
    echo "======================================"
    echo "Metagenomic Pipeline Setup"
    echo "======================================"
    echo
    
    # Check conda
    check_conda
    
    # Setup environments
    if [ "$ENVS_ONLY" = false ] || [ "$RUN_TEST" = false ]; then
        setup_environments
    fi
    
    # Setup databases
    if [ "$ENVS_ONLY" = false ]; then
        setup_databases
    fi
    
    # Update configuration
    if [ "$ENVS_ONLY" = false ]; then
        update_config
    fi
    
    # Create sample dataset
    if [ "$ENVS_ONLY" = false ]; then
        create_sample_dataset
    fi
    
    # Test pipeline
    if [ "$RUN_TEST" = true ]; then
        test_pipeline
    fi
    
    echo
    echo "======================================"
    echo "Setup completed successfully!"
    echo "======================================"
    echo
    echo "Next steps:"
    echo "1. Add your FASTQ files to the input directory"
    echo "2. Create/edit the sample sheet"
    echo "3. Run the pipeline with: ./run_pipeline.sh"
    echo
    echo "For help: ./run_pipeline.sh --help"
}

# Run main function
main "$@"