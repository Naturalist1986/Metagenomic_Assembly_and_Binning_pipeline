#!/bin/bash
#SBATCH --job-name=mag_map
#SBATCH --array=0-99
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=ofinkel

# map_to_selected_bins.sh - Competitive mapping: MAGs vs Other Bacterial vs Non-Bacterial
# Builds a combined reference from selected bins (MAGs) + EukFinder-classified contigs,
# maps reads ONCE competitively, and counts reads per category from the SAM stream.
# This ensures categories sum to 100% with no double-counting.

set -euo pipefail

# Source configuration and utilities
if [ -n "${PIPELINE_SCRIPT_DIR:-}" ]; then
    SCRIPT_DIR="$PIPELINE_SCRIPT_DIR"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${SCRIPT_DIR}/00_config_utilities.sh"
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

ARRAY_INDEX=${SLURM_ARRAY_TASK_ID:-0}

log "Starting competitive MAG mapping (Array index: $ARRAY_INDEX)"

# Calculate Java memory allocation
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    MEM_GB=$((SLURM_MEM_PER_NODE / 1024))
    JAVA_MEM=$((MEM_GB * 9 / 10))
else
    JAVA_MEM=96
fi
log "Using ${JAVA_MEM}G memory for Java processes"

# Define directories
SELECTED_BINS_DIR="${OUTPUT_DIR}/selected_bins"
EUKFINDER_DIR="${EUKFINDER_DIR:-${OUTPUT_DIR}/eukfinder_output}"
MAPPING_DIR="${OUTPUT_DIR}/selected_bins_mapping"

# Get list of treatments from selected_bins directory
TREATMENTS=($(ls -1 "$SELECTED_BINS_DIR" 2>/dev/null | sort))

if [ ${#TREATMENTS[@]} -eq 0 ]; then
    log "ERROR: No treatments found in $SELECTED_BINS_DIR"
    exit 1
fi

if [ $ARRAY_INDEX -ge ${#TREATMENTS[@]} ]; then
    log "ERROR: Array index $ARRAY_INDEX exceeds treatment count ${#TREATMENTS[@]}"
    exit 1
fi

TREATMENT="${TREATMENTS[$ARRAY_INDEX]}"
log "Processing treatment: $TREATMENT"

TREATMENT_DIR="${MAPPING_DIR}/${TREATMENT}"
mkdir -p "$TREATMENT_DIR"

# Clean up stale BBMap reference indexes
if [ -d "${SCRIPT_DIR}/ref" ]; then
    rm -rf "${SCRIPT_DIR}/ref" 2>/dev/null || true
fi

# Checkpoint check
SUMMARY_FILE="${TREATMENT_DIR}/mag_mapping_summary.txt"
if [ -f "$SUMMARY_FILE" ]; then
    log "Already completed for $TREATMENT (found summary file)"
    exit 0
fi

# ============================================================
# PHASE 1: Build combined tagged reference
# ============================================================
log "=== PHASE 1: Building combined tagged reference ==="

SELECTED_DIR="${SELECTED_BINS_DIR}/${TREATMENT}"
if [ ! -d "$SELECTED_DIR" ]; then
    log "ERROR: Selected bins directory not found: $SELECTED_DIR"
    exit 1
fi

COMBINED_REF="${TREATMENT_DIR}/combined_reference.fasta"
MAG_IDS="${TREATMENT_DIR}/mag_contig_ids.txt"
EUKFINDER_TAGGED="${TREATMENT_DIR}/eukfinder_tagged_temp.fasta"

> "$COMBINED_REF"
> "$MAG_IDS"

# --- 1a. Collect MAG contigs (tagged with MAG__ prefix) ---
log "Collecting MAG contigs from selected bins..."
bin_count=0

# Handle both coassembly (*.fa directly) and individual (sample_subdir/*.fa)
mag_files=()
for f in "${SELECTED_DIR}"/*.fa; do
    [ -f "$f" ] && mag_files+=("$f")
done

if [ ${#mag_files[@]} -eq 0 ]; then
    for sample_dir in "${SELECTED_DIR}"/*/; do
        [ -d "$sample_dir" ] || continue
        for f in "${sample_dir}"*.fa; do
            [ -f "$f" ] && mag_files+=("$f")
        done
    done
fi

for bin_fa in "${mag_files[@]}"; do
    # Tag contig headers with MAG__ prefix
    awk '/^>/ {print ">MAG__" substr($0,2); next} {print}' "$bin_fa" >> "$COMBINED_REF"
    # Extract original contig IDs for dedup later
    grep "^>" "$bin_fa" | awk '{sub(/^>/,""); print $1}' >> "$MAG_IDS"
    bin_count=$((bin_count + 1))
done

if [ $bin_count -eq 0 ]; then
    log "ERROR: No .fa files found for treatment $TREATMENT"
    exit 1
fi

sort -u -o "$MAG_IDS" "$MAG_IDS"
mag_contig_count=$(wc -l < "$MAG_IDS")
log "Found $bin_count MAG bins with $mag_contig_count unique contigs"

# --- 1b. Collect EukFinder contigs (Bact + non-Bact), tagged and deduped ---
log "Collecting EukFinder-classified contigs (excluding MAG contigs)..."

EUKFINDER_TREATMENT_DIR="${EUKFINDER_DIR}/${TREATMENT}"
> "$EUKFINDER_TAGGED"

REFINEMENT_TOOLS=("binette" "dastool" "binspreader" "selected")

if [ -d "$EUKFINDER_TREATMENT_DIR" ]; then
    for bin_dir in "${EUKFINDER_TREATMENT_DIR}"/*/; do
        [ -d "$bin_dir" ] || continue

        # Skip refinement tools (same filter as map_bacterial_vs_nonbacterial.sh)
        dir_name=$(basename "$bin_dir")
        is_refinement=false
        for rt in "${REFINEMENT_TOOLS[@]}"; do
            if [[ "$dir_name" == *"_${rt}_"* ]]; then
                is_refinement=true
                break
            fi
        done
        if [ "$is_refinement" = true ]; then
            continue
        fi

        results_dir="${bin_dir}Eukfinder_results"
        [ -d "$results_dir" ] || continue

        # Bacterial contigs → BACT__ tag (processed first for priority)
        bact_file=$(find "$results_dir" -maxdepth 1 \( -name "*.Bact.fasta" -o -name "Bact.fasta" \) -type f 2>/dev/null | head -1)
        if [ -n "${bact_file:-}" ] && [ -s "$bact_file" ]; then
            awk '/^>/ {print ">BACT__" substr($0,2); next} {print}' "$bact_file" >> "$EUKFINDER_TAGGED"
        fi

        # Non-bacterial contigs → NONBACT__ tag
        for category in Arch Euk Unk EUnk Misc; do
            cat_file=$(find "$results_dir" -maxdepth 1 \( -name "*.${category}.fasta" -o -name "${category}.fasta" \) -type f 2>/dev/null | head -1)
            if [ -n "${cat_file:-}" ] && [ -s "$cat_file" ]; then
                awk '/^>/ {print ">NONBACT__" substr($0,2); next} {print}' "$cat_file" >> "$EUKFINDER_TAGGED"
            fi
        done
    done

    # Deduplicate: remove contigs already in MAGs + remove duplicates across binners
    # BACT contigs appear first in the file, so they get priority over NONBACT for same contig
    log "Deduplicating EukFinder contigs (excluding MAG contigs)..."
    awk -v mag_file="$MAG_IDS" '
    BEGIN { while ((getline id < mag_file) > 0) seen[id] = 1 }
    /^>/ {
        header = $0
        tagged_first = $1; sub(/^>/, "", tagged_first)
        # Strip tag prefix to get original contig ID
        orig_id = tagged_first
        sub(/^MAG__/, "", orig_id)
        sub(/^BACT__/, "", orig_id)
        sub(/^NONBACT__/, "", orig_id)
        if (orig_id in seen) { skip = 1 }
        else { skip = 0; seen[orig_id] = 1; print header }
        next
    }
    !skip { print }
    ' "$EUKFINDER_TAGGED" >> "$COMBINED_REF"
else
    log "WARNING: No EukFinder results found at $EUKFINDER_TREATMENT_DIR"
    log "Reference will contain only MAG contigs (no Other Bact / Non-Bact)"
fi

# Clean up temp file
rm -f "$EUKFINDER_TAGGED"

# Count reference sequences by category
mag_ref=$(grep -c "^>MAG__" "$COMBINED_REF" || echo 0)
bact_ref=$(grep -c "^>BACT__" "$COMBINED_REF" || echo 0)
nonbact_ref=$(grep -c "^>NONBACT__" "$COMBINED_REF" || echo 0)
total_ref=$((mag_ref + bact_ref + nonbact_ref))

log "Combined reference: $total_ref contigs"
log "  MAG contigs: $mag_ref"
log "  Other Bacterial contigs: $bact_ref"
log "  Non-Bacterial contigs: $nonbact_ref"

# ============================================================
# PHASE 2: Find reads
# ============================================================
log "=== PHASE 2: Locating input reads ==="

READ_LOCATIONS=(
    "${OUTPUT_DIR}/coassembly/${TREATMENT}/merged_reads"
    "${OUTPUT_DIR}/assembly/${TREATMENT}"
    "${OUTPUT_DIR}/quality_filtering/${TREATMENT}"
    "${OUTPUT_DIR}/validate_repair/${TREATMENT}"
    "${INPUT_DIR}/${TREATMENT}"
)

R1_PATH=""
R2_PATH=""

for location in "${READ_LOCATIONS[@]}"; do
    if [ -d "$location" ]; then
        r1=$(find "$location" -maxdepth 1 \( -name "*R1*.fastq.gz" -o -name "*_1.fastq.gz" -o -name "filtered_1.fastq.gz" -o -name "merged_R1.fastq.gz" \) 2>/dev/null | head -1)
        r2=$(find "$location" -maxdepth 1 \( -name "*R2*.fastq.gz" -o -name "*_2.fastq.gz" -o -name "filtered_2.fastq.gz" -o -name "merged_R2.fastq.gz" \) 2>/dev/null | head -1)

        if [ -n "$r1" ] && [ -n "$r2" ] && [ -f "$r1" ] && [ -f "$r2" ]; then
            R1_PATH="$r1"
            R2_PATH="$r2"
            break
        fi
    fi
done

if [ -z "$R1_PATH" ] || [ -z "$R2_PATH" ]; then
    log "ERROR: Could not find R1/R2 reads for treatment $TREATMENT"
    log "Searched locations: ${READ_LOCATIONS[*]}"
    exit 1
fi

log "R1: $R1_PATH"
log "R2: $R2_PATH"

# ============================================================
# PHASE 3: Competitive mapping + per-category counting
# ============================================================
log "=== PHASE 3: Mapping reads to combined reference ==="

activate_env bbmap

COUNTS_FILE="${TREATMENT_DIR}/category_counts.txt"

# Map reads to combined reference, stream SAM through awk to count per category.
# BBMap writes SAM to stdout, statsfile to disk, logs to stderr.
# awk counts each mapped read by its reference contig tag (MAG__ / BACT__ / NONBACT__).
# Unmapped reads have RNAME=* and are not counted by awk (computed from total later).
bbmap.sh \
    in1="$R1_PATH" \
    in2="$R2_PATH" \
    ref="$COMBINED_REF" \
    out=stdout.sam \
    statsfile="${TREATMENT_DIR}/combined_stats.txt" \
    minid=0.95 \
    ambiguous=random \
    nodisk=true \
    threads=${SLURM_CPUS_PER_TASK:-32} \
    -Xmx${JAVA_MEM}g \
    2>"${TREATMENT_DIR}/combined_bbmap.log" \
| awk '
    BEGIN { mag=0; bact=0; nonbact=0 }
    /^@/ { next }
    {
        if ($3 ~ /^MAG__/) mag++
        else if ($3 ~ /^BACT__/) bact++
        else if ($3 ~ /^NONBACT__/) nonbact++
    }
    END {
        print "mag_reads=" mag
        print "bact_reads=" bact
        print "nonbact_reads=" nonbact
    }
' > "$COUNTS_FILE"

deactivate_env

# Verify counts file was produced
if [ ! -s "$COUNTS_FILE" ]; then
    log "ERROR: Category counts not produced. Check: ${TREATMENT_DIR}/combined_bbmap.log"
    exit 1
fi

# ============================================================
# PHASE 4: Compute statistics and write summary
# ============================================================
log "=== PHASE 4: Computing statistics ==="

# Source the per-category counts
source "$COUNTS_FILE"

# Get total reads from BBMap stats
total_reads=$(grep "^Reads Used:" "${TREATMENT_DIR}/combined_stats.txt" | \
    awk '{gsub(/,/,"",$3); sum+=$3} END {print sum+0}')
total_reads=${total_reads:-0}

# Compute unmapped
unmapped_reads=$((total_reads - mag_reads - bact_reads - nonbact_reads))

# Update counts file with all fields (used by summarize_with_mags.sh)
cat > "$COUNTS_FILE" << EOF
total_reads=$total_reads
mag_reads=$mag_reads
bact_reads=$bact_reads
nonbact_reads=$nonbact_reads
unmapped_reads=$unmapped_reads
EOF

# Calculate percentages
if [ "$total_reads" -gt 0 ]; then
    mag_pct=$(awk "BEGIN {printf \"%.2f\", ($mag_reads / $total_reads) * 100}")
    bact_pct=$(awk "BEGIN {printf \"%.2f\", ($bact_reads / $total_reads) * 100}")
    nonbact_pct=$(awk "BEGIN {printf \"%.2f\", ($nonbact_reads / $total_reads) * 100}")
    unmapped_pct=$(awk "BEGIN {printf \"%.2f\", ($unmapped_reads / $total_reads) * 100}")
else
    mag_pct="0.00"
    bact_pct="0.00"
    nonbact_pct="0.00"
    unmapped_pct="0.00"
fi

# Write human-readable summary (also serves as checkpoint)
cat > "$SUMMARY_FILE" << EOF
Competitive MAG Mapping Summary
Treatment: $TREATMENT
Generated: $(date)
===========================================

Input Reads:
  R1: $R1_PATH
  R2: $R2_PATH
  Total Reads: $total_reads

Reference (combined, tagged):
  MAG contigs (selected bins): $mag_ref  ($bin_count bins)
  Other Bacterial contigs (EukFinder Bact, non-MAG): $bact_ref
  Non-Bacterial contigs (EukFinder Arch+Euk+Unk+EUnk+Misc, non-MAG): $nonbact_ref
  Total reference contigs: $total_ref

Results (Single Competitive Mapping):
-----------------------------------

MAGs (Selected Bins):
  Mapped reads: $mag_reads
  Percentage: ${mag_pct}%

Other Bacterial:
  Mapped reads: $bact_reads
  Percentage: ${bact_pct}%

Non-Bacterial:
  Mapped reads: $nonbact_reads
  Percentage: ${nonbact_pct}%

Unmapped:
  Unmapped reads: $unmapped_reads
  Percentage: ${unmapped_pct}%
EOF

# Clean up large reference file
rm -f "$COMBINED_REF" "$MAG_IDS"

log "Competitive mapping complete for treatment $TREATMENT"
log "  Total reads: $total_reads"
log "  MAGs: $mag_reads (${mag_pct}%)"
log "  Other Bacterial: $bact_reads (${bact_pct}%)"
log "  Non-Bacterial: $nonbact_reads (${nonbact_pct}%)"
log "  Unmapped: $unmapped_reads (${unmapped_pct}%)"
log "Summary saved to: $SUMMARY_FILE"
exit 0
