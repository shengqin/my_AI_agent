#!/bin/bash
#SBATCH --job-name=cnvkit_ref
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssu42
#SBATCH --partition=emoding

# ==============================================================================
# Build Pooled Normal Reference for CNVkit (Hybrid Capture)
#
# This script processes ALL normal BAMs to create a single pooled reference .cnn
# file. Run this ONCE before submitting the per-sample mutation_calling_from_bam.sh
# job array.
#
# Usage:
#   sbatch build_cnvkit_reference.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED> [OUTPUT_DIR]
#
# Example:
#   sbatch build_cnvkit_reference.sh \
#       sample2barcode_normal.txt \
#       /oak/.../NovaSeqB23_1/Normal \
#       /oak/.../selectors/BLYMv2.bed \
#       cnvkit_pooled_ref
#
# Then submit the main pipeline with the pooled reference:
#   sbatch --array=1-N mutation_calling_from_bam.sh \
#       sample2barcode_tumor.txt sample2barcode_normal.txt \
#       /path/to/cfdna /path/to/normal results /path/to/beds \
#       cnvkit_pooled_ref/pooled_reference.cnn
# ==============================================================================

set -euo pipefail

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting pooled normal reference build..."

# ==============================================================================
# Input Validation
# ==============================================================================
if [[ "$#" -lt 3 ]]; then
    echo "Usage: sbatch build_cnvkit_reference.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED> [OUTPUT_DIR]"
    echo ""
    echo "Arguments:"
    echo "  NORMAL_SAMPLE_FILE  Tab-delimited file with normal sample names in column 1"
    echo "  NORMAL_DIR          Directory containing normal BAM files"
    echo "  TARGET_BED          Target BED file (e.g., BLYMv2 panel)"
    echo "  OUTPUT_DIR          Output directory for reference (default: cnvkit_pooled_ref)"
    exit 1
fi

NORMAL_SAMPLE_FILE="$1"
NORMAL_DIR="$2"
TARGET_BED="$3"
OUTPUT_DIR="${4:-cnvkit_pooled_ref}"

# Reference genome (same as mutation_calling_from_bam.sh)
REF_GENOME="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"

# BAM naming convention (must match mutation_calling_from_bam.sh)
NORMAL_PREFIX="Sample_"
NORMAL_SUFFIX="_Normal.dualindex-deduped.sorted.bam"

# Output reference file
POOLED_REF="$OUTPUT_DIR/pooled_reference.cnn"

# ==============================================================================
# Verify Inputs
# ==============================================================================
if [[ ! -f "$NORMAL_SAMPLE_FILE" ]]; then
    echo "ERROR: Normal sample file not found: $NORMAL_SAMPLE_FILE"
    exit 1
fi

if [[ ! -d "$NORMAL_DIR" ]]; then
    echo "ERROR: Normal BAM directory not found: $NORMAL_DIR"
    exit 1
fi

if [[ ! -f "$TARGET_BED" ]]; then
    echo "ERROR: Target BED file not found: $TARGET_BED"
    exit 1
fi

if [[ ! -f "$REF_GENOME" ]]; then
    echo "ERROR: Reference genome not found: $REF_GENOME"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# ==============================================================================
# Collect All Normal BAM Paths
# ==============================================================================
NORMAL_BAMS=()

while IFS=$'\t' read -r NORMAL_NAME _REST || [[ -n "$NORMAL_NAME" ]]; do
    # Strip carriage returns and whitespace
    NORMAL_NAME=$(echo "$NORMAL_NAME" | tr -d '\r' | xargs)
    [[ -z "$NORMAL_NAME" ]] && continue

    BAM_PATH="${NORMAL_DIR}/${NORMAL_PREFIX}${NORMAL_NAME}${NORMAL_SUFFIX}"

    if [[ -f "$BAM_PATH" ]]; then
        NORMAL_BAMS+=("$BAM_PATH")
        echo "  Found: $BAM_PATH"
    else
        echo "  WARNING: Missing BAM for $NORMAL_NAME: $BAM_PATH (skipping)"
    fi
done < "$NORMAL_SAMPLE_FILE"

if [[ ${#NORMAL_BAMS[@]} -lt 2 ]]; then
    echo "ERROR: Need at least 2 normal BAMs to build a pooled reference. Found: ${#NORMAL_BAMS[@]}"
    exit 1
fi

echo ""
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Building pooled reference from ${#NORMAL_BAMS[@]} normal samples..."

# ==============================================================================
# Build Pooled Reference with CNVkit
# ==============================================================================
source ~/miniconda3/etc/profile.d/conda.sh
conda activate cnvkit

if ! command -v cnvkit.py &> /dev/null; then
    echo "ERROR: cnvkit.py not found in conda env. Aborting."
    exit 1
fi

# cnvkit.py batch with --normal and no tumor builds only the reference.
# --output-reference specifies the pooled .cnn output file.
# --method hybrid uses both on-target and off-target bins (correct for BLYMv2).
cnvkit.py batch \
    --normal "${NORMAL_BAMS[@]}" \
    --targets "$TARGET_BED" \
    --fasta "$REF_GENOME" \
    --output-reference "$POOLED_REF" \
    --output-dir "$OUTPUT_DIR" \
    --method hybrid

conda deactivate

# ==============================================================================
# Verify Output
# ==============================================================================
if [[ -f "$POOLED_REF" ]]; then
    echo ""
    echo "============================================="
    echo "  Pooled Normal Reference Built Successfully"
    echo "============================================="
    echo "  Normal samples used: ${#NORMAL_BAMS[@]}"
    echo "  Reference file:      $POOLED_REF"
    echo "  File size:           $(du -h "$POOLED_REF" | cut -f1)"
    echo ""
    echo "  Next step: pass this reference to mutation_calling_from_bam.sh:"
    echo "    sbatch --array=1-N mutation_calling_from_bam.sh \\"
    echo "        sample2barcode_tumor.txt sample2barcode_normal.txt \\"
    echo "        /path/to/cfdna /path/to/normal results /path/to/beds \\"
    echo "        $POOLED_REF"
    echo "============================================="
else
    echo "ERROR: Pooled reference was not created. Check logs above."
    exit 1
fi
