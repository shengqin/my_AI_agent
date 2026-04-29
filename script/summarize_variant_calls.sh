#!/bin/bash
# ==============================================================================
# Cohort Summary for Somatic Variant Calling Results
# Aggregates Mutect2, Manta, and CNVkit outputs into reports.
#
# Usage:
#   ./summarize_variant_calls.sh <output_dir> <results_dir1> [results_dir2] [results_dir3] ...
#
# Example:
#   ./summarize_variant_calls.sh combined_summary \
#       /oak/.../NovaSeqB23_1/results_somatic_calling \
#       /oak/.../NovaSeqB23_2/results_somatic_calling \
#       /oak/.../NovaSeqB24_1/results_somatic_calling \
#       /oak/.../NovaSeqB24_2/results_somatic_calling
#
# Output files (created in output_dir):
#   1. cohort_variant_summary.csv       — Per-sample counts overview
#   2. mutect2_variants.tsv             — All variants with gene, effect, VAF
#   3. cnvkit_aberrations.tsv           — Significant CNV segments
#   4. manta_somatic_svs.tsv           — All somatic structural variants
#   5. cnvkit_gene_cnvs.tsv            — Gene-level copy number calls
# ==============================================================================

set -eu

# Portable zcat: macOS uses gzcat, Linux uses zcat
if command -v gzcat &> /dev/null; then
    ZCAT="gzcat"
else
    ZCAT="zcat"
fi

# Input validation
if [[ "$#" -lt 2 ]]; then
    echo "Usage: ./summarize_variant_calls.sh <output_dir> <results_dir1> [results_dir2] ..."
    echo "Example: ./summarize_variant_calls.sh combined_summary /path/to/B23_1/results /path/to/B23_2/results"
    exit 1
fi

OUTPUT_DIR="$1"
shift
RESULTS_DIRS=("$@")

# Validate all input directories exist
for DIR in "${RESULTS_DIRS[@]}"; do
    if [[ ! -d "$DIR" ]]; then
        echo "ERROR: Directory $DIR not found!"
        exit 1
    fi
done

mkdir -p "$OUTPUT_DIR"

# Output files
SUMMARY_CSV="$OUTPUT_DIR/cohort_variant_summary.csv"
MUTECT_DETAILS="$OUTPUT_DIR/mutect2_variants.tsv"
CNV_DETAILS="$OUTPUT_DIR/cnvkit_aberrations.tsv"
SV_DETAILS="$OUTPUT_DIR/manta_somatic_svs.tsv"
CNV_GENE_DETAILS="$OUTPUT_DIR/cnvkit_gene_cnvs.tsv"

echo "============================================="
echo "  Somatic Variant Calling Summary Generator"
echo "============================================="
echo "Output directory: $OUTPUT_DIR"
echo "Scanning ${#RESULTS_DIRS[@]} results director(ies):"
for DIR in "${RESULTS_DIRS[@]}"; do
    echo "  - $DIR"
done
echo ""

# ---- Write headers ----
echo "Sample,Batch,Pipeline_Status,Mutect2_Total,Mutect2_PASS,Manta_PASS_SVs,Manta_Total_SVs,CNVkit_Aberrant_Segments,CNVkit_Amplifications,CNVkit_Deletions" > "$SUMMARY_CSV"

echo -e "Sample\tBatch\tChromosome\tPosition\tRef\tAlt\tFilter\tGene\tVariant_Type\tProtein_Change\tTumor_AF\tTumor_AD_Ref\tTumor_AD_Alt\tTumor_DP\tNormal_AF\tNormal_DP" > "$MUTECT_DETAILS"

echo -e "Sample\tBatch\tChromosome\tStart\tEnd\tLog2_Ratio\tCopy_Number\tDepth\tP_Value\tProbes\tCall_Type" > "$CNV_DETAILS"

echo -e "Sample\tBatch\tChromosome\tPosition\tSV_Type\tFilter\tQual\tMate_Chrom\tMate_Pos" > "$SV_DETAILS"

echo -e "Sample\tBatch\tGene\tChromosome\tStart\tEnd\tLog2_Ratio\tDepth\tWeight\tProbes\tSegment_Log2\tSegment_Probes\tCall_Type" > "$CNV_GENE_DETAILS"

# ---- Process each results directory ----
SAMPLE_COUNT=0
COMPLETE_COUNT=0

for RESULTS_DIR in "${RESULTS_DIRS[@]}"; do
    # Extract batch name from path (e.g., "NovaSeqB23_1" from the path)
    BATCH_NAME=$(basename "$(dirname "$RESULTS_DIR")")

    for SAMPLE_DIR in "$RESULTS_DIR"/*/; do
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")

        # Skip non-sample directories (e.g., logs)
        [[ "$SAMPLE_NAME" == "logs" ]] && continue
        [[ ! -d "$SAMPLE_DIR" ]] && continue

    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

    # Initialize counters
    MUTECT_TOTAL="N/A"
    MUTECT_PASS="N/A"
    SV_PASS="N/A"
    SV_TOTAL="N/A"
    CNV_ABERRANT="N/A"
    CNV_AMP="N/A"
    CNV_DEL="N/A"
    STATUS="Incomplete"

    MANTA_OK=false
    MUTECT_OK=false
    CNVKIT_OK=false

    # ================================================================
    # 1. MUTECT2 — Parse annotated or filtered VCF
    # ================================================================
    # Prefer annotated VCF (has gene names); fall back to filtered VCF
    ANNOTATED_VCF=$(find "$SAMPLE_DIR/mutect2_somatic_run" -name "*_mutect2_annotated.vcf.gz" 2>/dev/null | head -n 1)
    FILTERED_VCF=$(find "$SAMPLE_DIR/mutect2_somatic_run" -name "*_mutect2_filtered.vcf.gz" 2>/dev/null | head -n 1)

    if [[ -f "$ANNOTATED_VCF" ]]; then
        MUTECT_OK=true
        USE_VCF="$ANNOTATED_VCF"

        MUTECT_TOTAL=$($ZCAT "$USE_VCF" | grep -v '^#' | wc -l | awk '{print $1}')
        MUTECT_PASS=$($ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')

        # SnpEff ANN field format: ANN=Allele|Effect|Impact|Gene|GeneID|...
        $ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" '{
            gene = "."; effect = "."; protein = ".";
            n = split($8, info_fields, ";");
            for (i = 1; i <= n; i++) {
                if (substr(info_fields[i], 1, 4) == "ANN=") {
                    ann = substr(info_fields[i], 5);
                    split(ann, ann_parts, "|");
                    effect = ann_parts[2];
                    gene = ann_parts[4];
                    if (length(ann_parts) >= 11 && ann_parts[11] != "") protein = ann_parts[11];
                    break;
                }
            }
            split($10, nf, ":"); split($11, tf, ":");
            split(tf[2], tad, ",");
            print sample"\t"batch"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7"\t"gene"\t"effect"\t"protein"\t"tf[3]"\t"tad[1]"\t"tad[2]"\t"tf[4]"\t"nf[3]"\t"nf[4]
        }' >> "$MUTECT_DETAILS"

    elif [[ -f "$FILTERED_VCF" ]]; then
        MUTECT_OK=true
        USE_VCF="$FILTERED_VCF"

        MUTECT_TOTAL=$($ZCAT "$USE_VCF" | grep -v '^#' | wc -l | awk '{print $1}')
        MUTECT_PASS=$($ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')

        $ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" '{
            split($10, nf, ":"); split($11, tf, ":");
            split(tf[2], tad, ",");
            print sample"\t"batch"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7"\t.\t.\t.\t"tf[3]"\t"tad[1]"\t"tad[2]"\t"tf[4]"\t"nf[3]"\t"nf[4]
        }' >> "$MUTECT_DETAILS"
    fi

    # ================================================================
    # 2. MANTA — Parse somatic SV VCF
    # ================================================================
    SV_VCF="$SAMPLE_DIR/manta_somatic_run/results/variants/somaticSV.vcf.gz"

    if [[ -f "$SV_VCF" ]]; then
        MANTA_OK=true

        SV_TOTAL=$($ZCAT "$SV_VCF" | grep -v '^#' | wc -l | awk '{print $1}')
        SV_PASS=$($ZCAT "$SV_VCF" | grep -v '^#' | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')

        if [[ "$SV_TOTAL" -gt 0 ]]; then
            $ZCAT "$SV_VCF" | grep -v '^#' | awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" '{
                svtype = "UNKNOWN"; mate_chr = "."; mate_pos = ".";
                if (match($8, /SVTYPE=[A-Z]+/)) svtype = substr($8, RSTART+7, RLENGTH-7);
                if (match($8, /CHR2=[^;]+/)) mate_chr = substr($8, RSTART+5, RLENGTH-5);
                if (match($8, /END=[0-9]+/)) mate_pos = substr($8, RSTART+4, RLENGTH-4);
                print sample"\t"batch"\t"$1"\t"$2"\t"svtype"\t"$7"\t"$6"\t"mate_chr"\t"mate_pos
            }' >> "$SV_DETAILS"
        fi
    fi

    # ================================================================
    # 3. CNVKIT — Parse call.cns (integer copy number calls)
    # ================================================================
    CALL_CNS=$(find "$SAMPLE_DIR/cnvkit_run" -name "*.call.cns" 2>/dev/null | head -n 1)

    if [[ -f "$CALL_CNS" ]]; then
        CNVKIT_OK=true

        CNV_AMP=$(awk -F'\t' 'NR>1 && $6 > 2 && $5 > 0.3 {c++} END {print c+0}' "$CALL_CNS")
        CNV_DEL=$(awk -F'\t' 'NR>1 && $6 < 2 && $5 < -0.3 {c++} END {print c+0}' "$CALL_CNS")
        CNV_ABERRANT=$((CNV_AMP + CNV_DEL))

        awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" 'NR>1 {
            call_type = "Neutral";
            if ($6 > 2 && $5 > 0.3) call_type = "Amplification";
            else if ($6 < 2 && $5 < -0.3) call_type = "Deletion";
            else next;
            print sample"\t"batch"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"call_type
        }' "$CALL_CNS" >> "$CNV_DETAILS"
    fi

    # ================================================================
    # 4. CNVKIT GENEMETRICS — Parse gene-level copy number calls
    # ================================================================
    GENEMETRICS_FILE=$(find "$SAMPLE_DIR/cnvkit_run" -name "*_genemetrics.tsv" 2>/dev/null | head -n 1)

    if [[ -f "$GENEMETRICS_FILE" ]]; then
        # genemetrics output columns: chromosome, start, end, gene, log2, depth, weight, probes, segment_log2, segment_probes
        # Filter to genes with abs(log2) > 0.2 (gain/loss threshold)
        awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" 'NR>1 && ($5 > 0.2 || $5 < -0.2) {
            call_type = "Neutral";
            if ($5 > 0.2) call_type = "Amplification";
            else if ($5 < -0.2) call_type = "Deletion";
            print sample"\t"batch"\t"$4"\t"$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"call_type
        }' "$GENEMETRICS_FILE" >> "$CNV_GENE_DETAILS"
    fi

    # ================================================================
    # Check pipeline completion
    # ================================================================
    if $MANTA_OK && $MUTECT_OK && $CNVKIT_OK; then
        STATUS="Complete"
        COMPLETE_COUNT=$((COMPLETE_COUNT + 1))
    fi

    # Append summary row
    echo "${SAMPLE_NAME},${BATCH_NAME},${STATUS},${MUTECT_TOTAL},${MUTECT_PASS},${SV_PASS},${SV_TOTAL},${CNV_ABERRANT},${CNV_AMP},${CNV_DEL}" >> "$SUMMARY_CSV"

    echo "  ✓ [$BATCH_NAME] $SAMPLE_NAME: Mutect2=$MUTECT_PASS PASS / $MUTECT_TOTAL total | Manta=$SV_PASS SVs | CNVkit=$CNV_ABERRANT aberrations ($CNV_AMP amp, $CNV_DEL del) [$STATUS]"

    done  # end sample loop
done  # end batch loop

echo ""
echo "============================================="
echo "  Summary Complete!"
echo "============================================="
echo "  Samples processed: $SAMPLE_COUNT"
echo "  Fully complete:    $COMPLETE_COUNT / $SAMPLE_COUNT"
echo ""
echo "  Output files:"
echo "    1. $SUMMARY_CSV"
echo "    2. $MUTECT_DETAILS"
echo "    3. $CNV_DETAILS"
echo "    4. $SV_DETAILS"
echo "    5. $CNV_GENE_DETAILS"
echo "============================================="
