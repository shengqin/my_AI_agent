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

echo -e "Sample\tBatch\tChromosome\tPosition\tRef\tAlt\tFilter\tGene\tVariant_Type\tProtein_Change\tTumor_AF\tTumor_AD_Ref\tTumor_AD_Alt\tTumor_DP\tNormal_AF\tNormal_DP\tTumor_F1R2\tTumor_F2R1" > "$MUTECT_DETAILS"

echo -e "Sample\tBatch\tChromosome\tStart\tEnd\tLog2_Ratio\tCopy_Number\tDepth\tWeight\tProbes\tP_Value\tCall_Type" > "$CNV_DETAILS"

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
    #
    # ROBUST PARSING: Detects tumor vs normal column order from VCF
    # #CHROM header and parses FORMAT fields by name, not position.
    # This prevents column misalignment when GATK changes output order.
    # ================================================================
    ANNOTATED_VCF=$(find "$SAMPLE_DIR/mutect2_somatic_run" -name "*_mutect2_annotated.vcf.gz" 2>/dev/null | head -n 1)
    FILTERED_VCF=$(find "$SAMPLE_DIR/mutect2_somatic_run" -name "*_mutect2_filtered.vcf.gz" 2>/dev/null | head -n 1)

    # Detect which VCF column is tumor vs normal from the #CHROM header.
    # GATK Mutect2 column order depends on -I argument order and version;
    # parsing the header is the only reliable way.
    _detect_sample_columns() {
        local vcf="$1" tumor_sm="$2" normal_sm="$3"
        local header_line
        header_line=$($ZCAT "$vcf" | grep '^#CHROM' | head -1)
        # VCF columns: 1=CHROM,2=POS,...,9=FORMAT,10+=samples
        # Find which column matches tumor and normal SM names
        TUMOR_COL=$( echo "$header_line" | awk -F'\t' -v sm="$tumor_sm" '{
            for(i=10;i<=NF;i++) if($i==sm){print i; exit}
        }')
        NORMAL_COL=$(echo "$header_line" | awk -F'\t' -v sm="$normal_sm" '{
            for(i=10;i<=NF;i++) if($i==sm){print i; exit}
        }')
        # Fallback: if SM not found (e.g. patched names), assume col10=tumor, col11=normal
        [[ -z "$TUMOR_COL" ]] && TUMOR_COL=10
        [[ -z "$NORMAL_COL" ]] && NORMAL_COL=11
    }

    # Read the tumor and normal SM tags from BAM headers (same as mutation_calling_from_bam.sh)
    # These are needed to match VCF column headers.
    _TUMOR_SM="$SAMPLE_NAME"
    _NORMAL_SM=""
    # Try to detect normal SM from the VCF header itself
    _get_normal_sm_from_vcf() {
        local vcf="$1" tumor_sm="$2"
        $ZCAT "$vcf" | grep '^#CHROM' | head -1 | awk -F'\t' -v tsm="$tumor_sm" '{
            for(i=10;i<=NF;i++) if($i!=tsm){print $i; exit}
        }'
    }

    if [[ -f "$ANNOTATED_VCF" ]]; then
        MUTECT_OK=true
        USE_VCF="$ANNOTATED_VCF"

        MUTECT_TOTAL=$($ZCAT "$USE_VCF" | grep -v '^#' | wc -l | awk '{print $1}')
        MUTECT_PASS=$($ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')

        # Detect sample columns
        _NORMAL_SM=$(_get_normal_sm_from_vcf "$USE_VCF" "$_TUMOR_SM")
        _detect_sample_columns "$USE_VCF" "$_TUMOR_SM" "$_NORMAL_SM"

        # Parse with FORMAT-field-name-based extraction (not positional).
        # Extracts AF, AD, DP by finding their index in the FORMAT string.
        $ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' \
            -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" \
            -v tcol="$TUMOR_COL" -v ncol="$NORMAL_COL" '{
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
            # Parse FORMAT field names to find AF, AD, DP, F1R2, F2R1 indices
            nf_count = split($9, fmt, ":");
            af_idx=0; ad_idx=0; dp_idx=0; f1_idx=0; f2_idx=0;
            for (i=1; i<=nf_count; i++) {
                if (fmt[i]=="AF") af_idx=i;
                if (fmt[i]=="AD") ad_idx=i;
                if (fmt[i]=="DP") dp_idx=i;
                if (fmt[i]=="F1R2") f1_idx=i;
                if (fmt[i]=="F2R1") f2_idx=i;
            }
            # Extract tumor and normal fields using detected column positions
            split($tcol, tf, ":");
            split($ncol, nf, ":");
            # Get values by FORMAT index
            t_af  = (af_idx>0) ? tf[af_idx] : ".";
            t_dp  = (dp_idx>0) ? tf[dp_idx] : ".";
            n_af  = (af_idx>0) ? nf[af_idx] : ".";
            n_dp  = (dp_idx>0) ? nf[dp_idx] : ".";
            t_ad  = (ad_idx>0) ? tf[ad_idx] : ".,."; split(t_ad, tad, ",");
            t_f1  = (f1_idx>0) ? tf[f1_idx] : ".,.";
            t_f2  = (f2_idx>0) ? tf[f2_idx] : ".,.";
            print sample"\t"batch"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7"\t"gene"\t"effect"\t"protein"\t"t_af"\t"tad[1]"\t"tad[2]"\t"t_dp"\t"n_af"\t"n_dp"\t"t_f1"\t"t_f2
        }' >> "$MUTECT_DETAILS"

    elif [[ -f "$FILTERED_VCF" ]]; then
        MUTECT_OK=true
        USE_VCF="$FILTERED_VCF"

        MUTECT_TOTAL=$($ZCAT "$USE_VCF" | grep -v '^#' | wc -l | awk '{print $1}')
        MUTECT_PASS=$($ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')

        _NORMAL_SM=$(_get_normal_sm_from_vcf "$USE_VCF" "$_TUMOR_SM")
        _detect_sample_columns "$USE_VCF" "$_TUMOR_SM" "$_NORMAL_SM"

        $ZCAT "$USE_VCF" | grep -v '^#' | awk -F'\t' \
            -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" \
            -v tcol="$TUMOR_COL" -v ncol="$NORMAL_COL" '{
            nf_count = split($9, fmt, ":");
            af_idx=0; ad_idx=0; dp_idx=0; f1_idx=0; f2_idx=0;
            for (i=1; i<=nf_count; i++) {
                if (fmt[i]=="AF") af_idx=i;
                if (fmt[i]=="AD") ad_idx=i;
                if (fmt[i]=="DP") dp_idx=i;
                if (fmt[i]=="F1R2") f1_idx=i;
                if (fmt[i]=="F2R1") f2_idx=i;
            }
            split($tcol, tf, ":");
            split($ncol, nf, ":");
            t_af  = (af_idx>0) ? tf[af_idx] : ".";
            t_dp  = (dp_idx>0) ? tf[dp_idx] : ".";
            n_af  = (af_idx>0) ? nf[af_idx] : ".";
            n_dp  = (dp_idx>0) ? nf[dp_idx] : ".";
            t_ad  = (ad_idx>0) ? tf[ad_idx] : ".,."; split(t_ad, tad, ",");
            t_f1  = (f1_idx>0) ? tf[f1_idx] : ".,.";
            t_f2  = (f2_idx>0) ? tf[f2_idx] : ".,.";
            print sample"\t"batch"\t"$1"\t"$2"\t"$4"\t"$5"\t"$7"\t.\t.\t.\t"t_af"\t"tad[1]"\t"tad[2]"\t"t_dp"\t"n_af"\t"n_dp"\t"t_f1"\t"t_f2
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
    #
    # ROBUST PARSING: Uses header-based column lookup instead of
    # positional indices to prevent column misalignment.
    # CNVkit .call.cns columns: chromosome,start,end,gene,log2,depth,probes,weight,cn
    # Also merges p-values from .bintest.cns if available.
    # ================================================================
    CALL_CNS=$(find "$SAMPLE_DIR/cnvkit_run" -name "*.call.cns" 2>/dev/null | head -n 1)
    BINTEST_CNS=$(find "$SAMPLE_DIR/cnvkit_run" -name "*.bintest.cns" 2>/dev/null | head -n 1)

    if [[ -f "$CALL_CNS" ]]; then
        CNVKIT_OK=true

        # Build a p-value lookup from bintest.cns (keyed by chr:start:end)
        PVAL_TMPFILE=""
        if [[ -f "$BINTEST_CNS" ]]; then
            PVAL_TMPFILE=$(mktemp)
            awk -F'\t' 'NR==1 {
                for(i=1;i<=NF;i++) col[$i]=i; next
            }
            {
                chr=( "chromosome" in col ? $col["chromosome"] : $1 )
                s=  ( "start"      in col ? $col["start"]      : $2 )
                e=  ( "end"        in col ? $col["end"]        : $3 )
                pv= ( "p-value"    in col ? $col["p-value"]    : "." )
                if (pv == "") pv = "."
                print chr"\t"s"\t"e"\t"pv
            }' "$BINTEST_CNS" > "$PVAL_TMPFILE"
        fi

        # Count aberrations using header-based column lookup
        CNV_AMP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
            $h["cn"]>2 && $h["log2"]>0.3{c++} END{print c+0}' "$CALL_CNS")
        CNV_DEL=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
            $h["cn"]<2 && $h["log2"]<-0.3{c++} END{print c+0}' "$CALL_CNS")
        CNV_ABERRANT=$((CNV_AMP + CNV_DEL))

        # Extract aberrant segments with header-based parsing
        awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" \
            -v pval_file="${PVAL_TMPFILE:-}" '
        NR==1 {
            for(i=1;i<=NF;i++) h[$i]=i
            # Load p-value lookup if file exists
            if (pval_file != "") {
                while ((getline line < pval_file) > 0) {
                    split(line, pf, "\t")
                    pval_lookup[pf[1]"\t"pf[2]"\t"pf[3]] = pf[4]
                }
                close(pval_file)
            }
            next
        }
        {
            chr     = $h["chromosome"]
            s       = $h["start"]
            e       = $h["end"]
            log2val = $h["log2"]
            depth   = $h["depth"]
            probes  = $h["probes"]
            weight  = $h["weight"]
            cn      = $h["cn"]

            call_type = "Neutral"
            if (cn > 2 && log2val > 0.3) call_type = "Amplification"
            else if (cn < 2 && log2val < -0.3) call_type = "Deletion"
            else next

            # Lookup p-value from bintest
            key = chr"\t"s"\t"e
            pv = (key in pval_lookup) ? pval_lookup[key] : "."

            print sample"\t"batch"\t"chr"\t"s"\t"e"\t"log2val"\t"cn"\t"depth"\t"weight"\t"probes"\t"pv"\t"call_type
        }' "$CALL_CNS" >> "$CNV_DETAILS"

        # Cleanup temp file
        [[ -n "$PVAL_TMPFILE" && -f "$PVAL_TMPFILE" ]] && rm -f "$PVAL_TMPFILE"
    fi

    # ================================================================
    # 4. CNVKIT GENEMETRICS — Parse gene-level copy number calls
    #
    # ROBUST PARSING: Uses header-based column lookup.
    # genemetrics columns: chromosome,start,end,gene,log2,depth,weight,probes,
    #                      segment_log2,segment_probes
    # ================================================================
    GENEMETRICS_FILE=$(find "$SAMPLE_DIR/cnvkit_run" -name "*_genemetrics.tsv" 2>/dev/null | head -n 1)

    if [[ -f "$GENEMETRICS_FILE" ]]; then
        awk -F'\t' -v sample="$SAMPLE_NAME" -v batch="$BATCH_NAME" '
        NR==1 {
            for(i=1;i<=NF;i++) h[$i]=i
            next
        }
        {
            gene   = $h["gene"]
            chr    = $h["chromosome"]
            s      = $h["start"]
            e      = $h["end"]
            log2v  = $h["log2"]
            depth  = $h["depth"]
            weight = $h["weight"]
            probes = $h["probes"]
            seg_l2 = ("segment_log2"   in h) ? $h["segment_log2"]   : "."
            seg_pr = ("segment_probes" in h) ? $h["segment_probes"] : "."

            if (log2v > 0.2) call_type = "Amplification"
            else if (log2v < -0.2) call_type = "Deletion"
            else next

            print sample"\t"batch"\t"gene"\t"chr"\t"s"\t"e"\t"log2v"\t"depth"\t"weight"\t"probes"\t"seg_l2"\t"seg_pr"\t"call_type
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
