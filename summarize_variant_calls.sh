#!/bin/bash
# ==============================================================================
# Cohort Aggregator for Somatic Variant Calling
# Extracts counts of high-confidence (PASS) mutations across Manta, Mutect2, and CNVkit.
# ==============================================================================

# Default input directory is results_somatic_calling
OUTPUT_DIR="${1:-results_somatic_calling}"
SUMMARY_FILE="${2:-cohort_variant_summary.csv}"

if [[ ! -d "$OUTPUT_DIR" ]]; then
    echo "ERROR: Directory $OUTPUT_DIR not found!"
    echo "Usage: ./summarize_variant_calls.sh [results_directory] [output_csv_name]"
    exit 1
fi

echo "Generating variant summary..."

# Write CSV Header
echo "Sample,Target_Capture_Pipeline_Status,Mutect2_PASS_SNVs,Manta_PASS_Structural_Variants,CNVkit_Significant_Aberrations" > "$SUMMARY_FILE"

# Iterate through each sample directory
for SAMPLE_DIR in "$OUTPUT_DIR"/*/; do
    
    # Strip the trailing slash to get the sample name
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    
    # Initialize basic counters
    SNV_COUNT="N/A"
    SV_COUNT="N/A"
    CNV_COUNT="N/A"
    STATUS="Incomplete"
    
    # Try to locate the relevant output files
    SNV_VCF=$(find "$SAMPLE_DIR/mutect2_somatic_run" -name "*_mutect2_filtered.vcf.gz" 2>/dev/null | head -n 1)
    SV_VCF="$SAMPLE_DIR/manta_somatic_run/results/variants/somaticSV.vcf.gz"
    CNV_CNS=$(find "$SAMPLE_DIR/cnvkit_run" -name "*.cns" 2>/dev/null | head -n 1)
    
    ALL_EXIST=true
    
    # Mutect2 Parsing (Count PASS records)
    if [[ -f "$SNV_VCF" ]]; then
        SNV_COUNT=$(zgrep -v '^#' "$SNV_VCF" | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')
    else
        ALL_EXIST=false
    fi
    
    # Manta Parsing (Count PASS somatic records)
    if [[ -f "$SV_VCF" ]]; then
        SV_COUNT=$(zgrep -v '^#' "$SV_VCF" | awk -F'\t' '$7 == "PASS" {c++} END {print c+0}')
    else
        ALL_EXIST=false
    fi
    
    # CNVkit Parsing (Count segments where LogRatio translates to notable Amp/Del)
    # E.g., Log2Ratio > 0.5 or < -0.5
    if [[ -f "$CNV_CNS" ]]; then
        CNV_COUNT=$(awk -F'\t' 'NR>1 { if ($5 > 0.5 || $5 < -0.5) count++ } END {print count+0}' "$CNV_CNS")
    else
        ALL_EXIST=false
    fi
    
    # Check if the pipeline completely finished for this sample
    if $ALL_EXIST; then 
        STATUS="Complete"
    fi
    
    # Append the row to the CSV
    echo "${SAMPLE_NAME},${STATUS},${SNV_COUNT},${SV_COUNT},${CNV_COUNT}" >> "$SUMMARY_FILE"

done

echo "Successfully parsed cohort files!"
echo "Saved summary to: $SUMMARY_FILE"
