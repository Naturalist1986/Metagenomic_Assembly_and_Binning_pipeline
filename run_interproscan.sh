#!/bin/bash
#SBATCH --job-name=IPR_Local
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --array=0-162          # Set this to 0-(number of files minus 1)

# --- 1. SET THE PATH TO YOUR INPUT FILES ---
# Provide the path where your fasta files are located
INPUT_DIR="/sci/backup/ofinkel/moshea/Efrat_Metagenomes_Novogene/coassembly/binning/mtz/comebin/comebin_res/comebin_res_bins/Proteins"
FILES=($INPUT_DIR/*.faa)

# --- 2. GET CURRENT FILE AND ITS DIRECTORY ---
INPUT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
# Get the absolute path to the directory containing the file
FILE_DIR=$(dirname "$(realpath "$INPUT_FILE")")
FILE_NAME=$(basename "$INPUT_FILE")

# --- 3. DEFINE LOG AND TMP LOCATIONS ---
# This puts the Slurm logs and InterProScan temp files in the same folder
LOG_OUT="$FILE_DIR/${FILE_NAME}.slurm.out"
LOG_ERR="$FILE_DIR/${FILE_NAME}.slurm.err"
IPR_TMP="$FILE_DIR/tmp_ipr_${SLURM_ARRAY_TASK_ID}"

# Since SBATCH directives can't use variables inside the script easily,
# we manually redirect stdout and stderr for the main execution:
exec > "$LOG_OUT" 2> "$LOG_ERR"

echo "Processing file: $INPUT_FILE"
echo "Output and Tmp will be in: $FILE_DIR"

# Ensure temp directory exists
mkdir -p "$IPR_TMP"

# Define a path for the cleaned version of the file
CLEAN_FILE="$FILE_DIR/${FILE_NAME}.cleaned.fasta"

echo "Cleaning asterisks from $FILE_NAME..."

# Use sed to remove '*' only from sequence lines (lines not starting with '>')
# or simply remove all '*' if you are certain they aren't in the headers.
sed 's/\*//g' "$INPUT_FILE" > "$CLEAN_FILE"

# Now run InterProScan using the CLEAN_FILE
/sci/backup/ofinkel/moshea/Tools/interproscan-5.76-107.0/interproscan.sh \
    -i "$CLEAN_FILE" \
    -b "$FILE_DIR/${FILE_NAME%.*}" \
    -T "$IPR_TMP" \
    -cpu $SLURM_CPUS_PER_TASK \
    -dp \
    --goterms --pa

# Cleanup: remove the temporary cleaned file
rm "$CLEAN_FILE"

echo "Job finished at $(date)"