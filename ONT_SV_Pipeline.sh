#!/usr/bin/env bash
# =============================================================================
# Oxford Nanopore Long-Read SV/CNV Pipeline
# -----------------------------------------------------------------------------
# Workflow:
#   FASTQ concatenation → read length filtering (filtlong) →
#   alignment (minimap2) → BAM conversion, sorting & indexing (samtools) →
#   alignment QC (samtools flagstat) → structural variant calling (Sniffles2)
#   → VCF stats summary (bcftools)
#
# Dependencies (conda environment recommended — see Environment Setup below):
#   minimap2, samtools, sniffles, filtlong, bcftools, moreutils (for 'ts')
#
# Usage:
#   chmod +x ont_sv_pipeline.sh
#   ./ont_sv_pipeline.sh
#
# To run in the background (persists after terminal closes):
#   nohup bash ont_sv_pipeline.sh &
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION — edit these variables before running
# =============================================================================

# Directory containing input .fastq or .fastq.gz files
FASTQ_DIR="."

# Reference genome FASTA
REFERENCE="hg38.fa"

# Number of threads for minimap2 and samtools
THREADS=16

# Minimum read length to retain after filtering (bp)
MIN_READ_LENGTH=1000

# Minimum acceptable alignment rate (%) — pipeline exits if below this value
MIN_ALIGNMENT_RATE=70

# Minimum disk space required to start (GB)
MIN_DISK_GB=50

# Output prefix for all generated files
OUT_PREFIX="output"

# =============================================================================
# DERIVED FILE NAMES — do not edit below this line
# =============================================================================

CONCAT_FASTQ="${OUT_PREFIX}_concatenated.fastq"
FILTERED_FASTQ="${OUT_PREFIX}_filtered.fastq"
SAM="${OUT_PREFIX}_mapped.sam"
BAM="${OUT_PREFIX}_genomic.bam"
SORTED_BAM="${OUT_PREFIX}_genomic_sorted.bam"
FLAGSTAT="${OUT_PREFIX}_flagstat.txt"
VCF="${OUT_PREFIX}_structural_variants.vcf"
VCF_STATS="${OUT_PREFIX}_vcf_stats.txt"
LOG="pipeline_$(date +%Y%m%d_%H%M%S).log"

# =============================================================================
# LOGGING SETUP
# Timestamps all stdout and stderr to both terminal and log file
# Requires 'ts' from the moreutils package: conda install moreutils
# =============================================================================

exec &> >(ts '[%Y-%m-%d %H:%M:%S]' | tee "$LOG")

echo "========================================"
echo " ONT Long-Read SV/CNV Pipeline"
echo " Log: $LOG"
echo "========================================"
echo ""

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================

echo "[Pre-flight] Checking dependencies..."

MISSING=()
for tool in minimap2 samtools sniffles filtlong bcftools; do
    if ! command -v "$tool" &>/dev/null; then
        MISSING+=("$tool")
    fi
done

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "ERROR: The following tools were not found on PATH:" >&2
    printf '  - %s\n' "${MISSING[@]}" >&2
    echo "Activate your conda environment first: conda activate ont_sv" >&2
    exit 1
fi

echo "  All dependencies found."

# ── Disk space check ──────────────────────────────────────────────────────────
echo "[Pre-flight] Checking available disk space..."

AVAIL_GB=$(df -BG . | awk 'NR==2 {gsub("G",""); print $4}')
if [[ "$AVAIL_GB" -lt "$MIN_DISK_GB" ]]; then
    echo "ERROR: Only ${AVAIL_GB}GB available. Pipeline requires at least ${MIN_DISK_GB}GB." >&2
    exit 1
fi

echo "  ${AVAIL_GB}GB available (minimum: ${MIN_DISK_GB}GB). OK."

# ── Reference genome check ────────────────────────────────────────────────────
if [[ ! -f "$REFERENCE" ]]; then
    echo "ERROR: Reference genome not found at '$REFERENCE'." >&2
    echo "Download with:" >&2
    echo "  wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz" >&2
    exit 1
fi

echo ""

# =============================================================================
# 1. CONCATENATE FASTQ FILES
# =============================================================================

echo "[Step 1/6] Concatenating FASTQ files..."

if [[ -f "$CONCAT_FASTQ" ]]; then
    echo "  Skipping — $CONCAT_FASTQ already exists."
else
    GZ_COUNT=$(find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq.gz" | wc -l)
    if [[ "$GZ_COUNT" -gt 0 ]]; then
        echo "  Decompressing $GZ_COUNT .fastq.gz files..."
        gunzip "$FASTQ_DIR"/*.fastq.gz
    fi

    FASTQ_COUNT=$(find "$FASTQ_DIR" -maxdepth 1 -name "*.fastq" \
        ! -name "$CONCAT_FASTQ" \
        ! -name "$FILTERED_FASTQ" | wc -l)

    if [[ "$FASTQ_COUNT" -eq 0 ]]; then
        echo "ERROR: No .fastq files found in $FASTQ_DIR" >&2
        exit 1
    fi

    echo "  Concatenating $FASTQ_COUNT FASTQ files → $CONCAT_FASTQ"
    cat "$FASTQ_DIR"/*.fastq > "$CONCAT_FASTQ"
fi

echo ""

# =============================================================================
# 2. READ LENGTH FILTERING (filtlong)
# =============================================================================
# Removes reads shorter than MIN_READ_LENGTH bp. Very short ONT reads
# align poorly and can introduce noise into SV calls.

echo "[Step 2/6] Filtering reads by minimum length (${MIN_READ_LENGTH}bp)..."

if [[ -f "$FILTERED_FASTQ" ]]; then
    echo "  Skipping — $FILTERED_FASTQ already exists."
else
    filtlong --min_length "$MIN_READ_LENGTH" "$CONCAT_FASTQ" > "$FILTERED_FASTQ"

    BEFORE=$(awk 'NR%4==1' "$CONCAT_FASTQ" | wc -l)
    AFTER=$(awk 'NR%4==1' "$FILTERED_FASTQ" | wc -l)
    echo "  Reads before filtering : $BEFORE"
    echo "  Reads after filtering  : $AFTER"
    echo "  Reads discarded        : $(( BEFORE - AFTER ))"
fi

echo ""

# =============================================================================
# 3. ALIGN TO REFERENCE GENOME (minimap2)
# =============================================================================
# --MD           : adds MD tag (required by some downstream annotation tools)
# -a             : output in SAM format
# -t             : number of threads
# --secondary=no : suppress secondary alignments (not used by Sniffles,
#                  reduces BAM size significantly)

echo "[Step 3/6] Aligning reads with minimap2..."

if [[ -f "$SORTED_BAM" ]]; then
    echo "  Skipping — $SORTED_BAM already exists."
else
    echo "  Aligning $FILTERED_FASTQ → $SAM"
    minimap2 \
        --MD \
        --secondary=no \
        -t "$THREADS" \
        -a "$REFERENCE" \
        "$FILTERED_FASTQ" \
        > "$SAM"

    echo ""

    # ── Convert, sort, index ──────────────────────────────────────────────────

    echo "[Step 4/6] Converting SAM → BAM, sorting, and indexing..."

    echo "  Converting SAM → BAM..."
    samtools view -bS -@ "$THREADS" "$SAM" > "$BAM"

    echo "  Sorting BAM → $SORTED_BAM"
    samtools sort -@ "$THREADS" "$BAM" -o "$SORTED_BAM"

    echo "  Indexing sorted BAM..."
    samtools index "$SORTED_BAM"

    echo "  Removing intermediate SAM/BAM files..."
    rm -f "$SAM" "$BAM"
fi

echo ""

# =============================================================================
# 4. ALIGNMENT QC (samtools flagstat)
# =============================================================================
# Saves alignment statistics and gates the pipeline on a minimum alignment
# rate — exits before SV calling if alignment quality is too low.

echo "[Step 5/6] Running alignment QC (samtools flagstat)..."

samtools flagstat "$SORTED_BAM" > "$FLAGSTAT"
cat "$FLAGSTAT"

# Parse overall alignment rate and compare to threshold
ALIGN_RATE=$(grep -m1 "mapped (" "$FLAGSTAT" \
    | grep -oP '\d+\.\d+(?=%)' \
    | head -1)

echo ""
echo "  Overall alignment rate: ${ALIGN_RATE}%"

if (( $(echo "$ALIGN_RATE < $MIN_ALIGNMENT_RATE" | bc -l) )); then
    echo "ERROR: Alignment rate ${ALIGN_RATE}% is below the minimum threshold of ${MIN_ALIGNMENT_RATE}%." >&2
    echo "Check read quality, reference genome, or minimap2 parameters before proceeding." >&2
    exit 1
fi

echo "  Alignment rate OK (threshold: ${MIN_ALIGNMENT_RATE}%)."
echo ""

# =============================================================================
# 5. STRUCTURAL VARIANT CALLING (Sniffles2)
# =============================================================================
# Sniffles2 infers SVs from alignment signatures in the BAM file.
# No reference genome is required at this step.
#
# Sniffles2 syntax (v2.x):
#   --input   : sorted BAM file
#   --vcf     : output VCF path
#   --threads : number of threads
#
# Note: if using Sniffles v1.x, replace with:
#   sniffles -m "$SORTED_BAM" -v "$VCF"

echo "[Step 5/6] Calling structural variants with Sniffles2..."

if [[ -f "$VCF" ]]; then
    echo "  Skipping — $VCF already exists."
else
    sniffles \
        --input "$SORTED_BAM" \
        --vcf "$VCF" \
        --threads "$THREADS"
fi

echo ""

# =============================================================================
# 6. VCF SUMMARY STATISTICS (bcftools)
# =============================================================================
# Generates a per-SV-type breakdown (DEL, INS, DUP, INV, BND) and general
# variant statistics.

echo "[Step 6/6] Generating VCF summary statistics..."

bcftools stats "$VCF" > "$VCF_STATS"

echo "  SV type breakdown:"
grep "^SN" "$VCF_STATS" | awk -F'\t' '{printf "    %-40s %s\n", $3, $4}'

echo ""
echo "  Full stats written to: $VCF_STATS"
echo ""
echo "========================================"
echo " Pipeline complete."
echo " Sorted BAM : $SORTED_BAM"
echo " Flagstat   : $FLAGSTAT"
echo " VCF        : $VCF"
echo " VCF stats  : $VCF_STATS"
echo " Log        : $LOG"
echo "========================================"
