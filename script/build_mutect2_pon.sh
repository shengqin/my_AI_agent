#!/bin/bash
#SBATCH --job-name=build_pon
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssu42
#SBATCH --partition=emoding

# ==============================================================================
# Build Mutect2 Panel of Normals (PoN) + PoN Frequency Summary
#
# Creates two outputs:
#   1. pon.vcf.gz              — GATK PoN for Mutect2 --panel-of-normals
#   2. pon_site_summary.tsv    — Per-site PoN frequency table for R filtering
#
# Usage:
#   sbatch build_mutect2_pon.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED> [OUTPUT_DIR]
#
# Then pass outputs to the main pipeline:
#   sbatch --array=1-N mutation_calling_from_bam.sh \
#       ... <POOLED_CNVKIT_REF> <PON_VCF>
# ==============================================================================

set -euo pipefail

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Mutect2 Panel of Normals build..."

# ==============================================================================
# Environment
# ==============================================================================
ml biology
ml samtools/1.16.1
ml gatk/4.6.0.0

# ==============================================================================
# Input Validation
# ==============================================================================
if [[ "$#" -lt 3 ]]; then
    echo "Usage: sbatch build_mutect2_pon.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED> [OUTPUT_DIR]"
    echo ""
    echo "Arguments:"
    echo "  NORMAL_SAMPLE_FILE  Tab-delimited file with normal sample names in column 1"
    echo "  NORMAL_DIR          Directory containing normal BAM files"
    echo "  TARGET_BED          Target BED file (e.g., BLYMv2 panel)"
    echo "  OUTPUT_DIR          Output directory (default: mutect2_pon)"
    exit 1
fi

NORMAL_SAMPLE_FILE="$1"
NORMAL_DIR="$2"
TARGET_BED="$3"
OUTPUT_DIR="${4:-mutect2_pon}"

REF_GENOME="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"
NORMAL_PREFIX="Sample_"
NORMAL_SUFFIX="_Normal.dualindex-deduped.sorted.bam"

NORMAL_VCFS_DIR="$OUTPUT_DIR/normal_vcfs"
PON_VCF="$OUTPUT_DIR/pon.vcf.gz"
PON_SUMMARY="$OUTPUT_DIR/pon_site_summary.tsv"

# Verify inputs
for FILE in "$NORMAL_SAMPLE_FILE" "$TARGET_BED" "$REF_GENOME"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: File not found: $FILE"
        exit 1
    fi
done
if [[ ! -d "$NORMAL_DIR" ]]; then
    echo "ERROR: Directory not found: $NORMAL_DIR"
    exit 1
fi

mkdir -p "$NORMAL_VCFS_DIR"

# ==============================================================================
# Step 1: Run Mutect2 in Tumor-Only Mode on Each Normal
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 1: Running Mutect2 tumor-only on each normal..."

N_NORMALS=0

while IFS=$'\t' read -r NORMAL_NAME _REST || [[ -n "$NORMAL_NAME" ]]; do
    NORMAL_NAME=$(echo "$NORMAL_NAME" | tr -d '\r' | xargs)
    [[ -z "$NORMAL_NAME" ]] && continue

    BAM_PATH="${NORMAL_DIR}/${NORMAL_PREFIX}${NORMAL_NAME}${NORMAL_SUFFIX}"
    OUT_VCF="${NORMAL_VCFS_DIR}/${NORMAL_NAME}.vcf.gz"

    if [[ ! -f "$BAM_PATH" ]]; then
        echo "  WARNING: BAM not found for $NORMAL_NAME: $BAM_PATH (skipping)"
        continue
    fi

    N_NORMALS=$((N_NORMALS + 1))

    if [[ -f "$OUT_VCF" ]]; then
        echo "  Already done: $NORMAL_NAME"
        continue
    fi

    echo "  Processing: $NORMAL_NAME"

    # Extract SM tag; patch if missing (same logic as main pipeline)
    set +o pipefail
    SM_TAG=$(samtools view -H "$BAM_PATH" | grep -m1 '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/')
    set -o pipefail

    ACTUAL_BAM="$BAM_PATH"
    if [[ -z "$SM_TAG" ]]; then
        SM_TAG="$NORMAL_NAME"
        PATCHED_BAM="${NORMAL_VCFS_DIR}/${NORMAL_NAME}_rg_patched.bam"
        if [[ ! -f "${PATCHED_BAM}.bai" ]]; then
            samtools addreplacerg -r "@RG\tID:${NORMAL_NAME}\tSM:${NORMAL_NAME}\tLB:Normal\tPL:ILLUMINA" \
                -o "$PATCHED_BAM" -@ 4 "$BAM_PATH"
            samtools index -@ 4 "$PATCHED_BAM"
        fi
        ACTUAL_BAM="$PATCHED_BAM"
    fi

    gatk Mutect2 \
        -R "$REF_GENOME" \
        -I "$ACTUAL_BAM" \
        --max-mnp-distance 0 \
        -L "$TARGET_BED" \
        -O "$OUT_VCF"

done < "$NORMAL_SAMPLE_FILE"

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 1 complete. Processed $N_NORMALS normals."

if [[ "$N_NORMALS" -lt 2 ]]; then
    echo "ERROR: Need at least 2 normals to build a PoN. Found: $N_NORMALS"
    exit 1
fi

# ==============================================================================
# Step 2: Create GenomicsDB Workspace
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 2: Creating GenomicsDB workspace..."

GENOMICSDB_DIR="$OUTPUT_DIR/pon_genomicsdb"

V_ARGS=""
for vcf in "$NORMAL_VCFS_DIR"/*.vcf.gz; do
    [[ "$vcf" == *"_rg_patched"* ]] && continue
    V_ARGS="$V_ARGS -V $vcf"
done

# Remove existing workspace if present (GenomicsDBImport requires clean dir)
if [[ -d "$GENOMICSDB_DIR" ]]; then
    echo "  Removing existing GenomicsDB workspace..."
    rm -rf "$GENOMICSDB_DIR"
fi

gatk GenomicsDBImport \
    -R "$REF_GENOME" \
    -L "$TARGET_BED" \
    $V_ARGS \
    --genomicsdb-workspace-path "$GENOMICSDB_DIR"

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 2 complete."

# ==============================================================================
# Step 3: Create Panel of Normals VCF
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 3: Building Panel of Normals..."

gatk CreateSomaticPanelOfNormals \
    -R "$REF_GENOME" \
    -V "gendb://$GENOMICSDB_DIR" \
    -O "$PON_VCF"

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 3 complete: $PON_VCF"

# ==============================================================================
# Step 4: Generate PoN Site Frequency Summary for R-Level Filtering
#
# For each variant site across all normal VCFs, records:
#   - How many PoN normals have the variant (at AF >= 0.001 / 0.1%)
#   - The maximum AF observed across all normals
#   - The fraction of PoN samples with the variant
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 4: Generating PoN frequency summary..."

PON_RAW="$OUTPUT_DIR/pon_raw_variants.tmp"
> "$PON_RAW"

for vcf in "$NORMAL_VCFS_DIR"/*.vcf.gz; do
    [[ "$vcf" == *"_rg_patched"* ]] && continue
    # Extract all variant sites and their AF (unfiltered — we want everything)
    # bcftools query extracts FORMAT/AF for the single sample
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' "$vcf" 2>/dev/null >> "$PON_RAW" || true
done

# Sort and aggregate: count samples per site, track max AF
# Only count samples where AF >= 0.001 (0.1%) per the paper's threshold
sort -k1,1V -k2,2n -k3,3 -k4,4 "$PON_RAW" | \
awk -F'\t' -v n="$N_NORMALS" '
BEGIN {
    print "Chromosome\tPosition\tRef\tAlt\tN_PON_Samples\tMax_PON_AF\tFraction_PON_Samples"
}
{
    key = $1 "\t" $2 "\t" $3 "\t" $4
    af = $5 + 0  # coerce to number

    if (key != prev_key && NR > 1) {
        if (count > 0) {
            printf "%s\t%d\t%.6f\t%.4f\n", prev_key, count, max_af, count / n
        }
        count = 0
        max_af = 0
    }

    # Count only if AF >= 0.1% (paper threshold)
    if (af >= 0.001) count++
    if (af > max_af) max_af = af
    prev_key = key
}
END {
    if (count > 0) {
        printf "%s\t%d\t%.6f\t%.4f\n", prev_key, count, max_af, count / n
    }
}' > "$PON_SUMMARY"

rm -f "$PON_RAW"

PON_SITES=$(tail -n +2 "$PON_SUMMARY" | wc -l | awk '{print $1}')
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Step 4 complete: $PON_SITES sites in $PON_SUMMARY"

# ==============================================================================
# Done
# ==============================================================================
echo ""
echo "============================================="
echo "  Mutect2 Panel of Normals Built Successfully"
echo "============================================="
echo "  Normal samples:       $N_NORMALS"
echo "  PoN VCF:              $PON_VCF"
echo "  PoN site summary:     $PON_SUMMARY (${PON_SITES} sites)"
echo ""
echo "  Next: pass to mutation_calling_from_bam.sh as arg \$8:"
echo "    sbatch --array=1-N mutation_calling_from_bam.sh \\"
echo "        sample2barcode_tumor.txt sample2barcode_normal.txt \\"
echo "        /path/to/cfdna /path/to/normal results /path/to/beds \\"
echo "        <cnvkit_pooled_ref.cnn> $PON_VCF"
echo "============================================="
