#!/bin/bash
#SBATCH --job-name=somatic
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssu42
#SBATCH --partition=emoding

# ==============================================================================
# Somatic Variant Calling Workflow (cfDNA + matched normal)
# Designed for Stanford Sherlock HPC Cluster
# 
# Usage: 
#   1. Set the #SBATCH --array directive above to match your sample count.
#      (e.g., if you have 45 samples, use --array=1-45)
#   2. Submit: sbatch mutation_calling_from_bam.sh
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
# Treat unset variables as an error when substituting
# Fail on any error in a pipeline
set -euo pipefail

# ==============================================================================
# Environment Setup
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Loading modules..."

# Load required tools via Sherlock's module system
ml biology
ml samtools/1.16.1
ml bcftools/1.16
ml python/2.7.13
ml gatk/4.6.0.0
# Manta is called directly via absolute path below (no module needed)

# Resolve htslib tools to absolute paths NOW, from the loaded modules, before any conda env
# activation can shadow them on PATH. The snpeff conda env may not ship bgzip/tabix, so we
# call these captured binaries explicitly during annotation.
BGZIP_BIN="$(command -v bgzip || true)"
TABIX_BIN="$(command -v tabix || true)"
if [[ -z "$BGZIP_BIN" || -z "$TABIX_BIN" ]]; then
    echo "ERROR: bgzip/tabix not found after loading modules (need htslib via samtools/bcftools)."
    exit 1
fi

# ==============================================================================
# Input Variables
# ==============================================================================
if [[ "$#" -lt 4 ]]; then
    echo "Usage: sbatch --array=1-N mutation_calling_from_bam.sh <CFDNA_SAMPLE_FILE> <NORMAL_SAMPLE_FILE> <CFDNA_DIR> <NORMAL_DIR> [OUTPUT_DIR] [BED_DIR] [POOLED_REF] [MUTECT2_PON]"
    echo "Example: sbatch --array=1-5 mutation_calling_from_bam.sh sample2barcode_tumor.txt sample2barcode_normal.txt /path/to/cfdna /path/to/normal results /path/to/beds pooled_reference.cnn mutect2_pon/pon.vcf.gz"
    exit 1
fi

# 1. Path to the sample_to_barcode files
CFDNA_SAMPLE_FILE="$1"
NORMAL_SAMPLE_FILE="$2"

# 2. Directory containing cfDNA BAM files
CFDNA_DIR="$3"
CFDNA_PREFIX="Sample_"
CFDNA_SUFFIX="_cfDNA.dualindex-deduped.sorted.bam"

# 3. Directory containing matched normal BAM files
NORMAL_DIR="$4"
NORMAL_PREFIX="Sample_"
NORMAL_SUFFIX="_Normal.dualindex-deduped.sorted.bam"

# 4. Output Directory for all variant calling results
OUTPUT_DIR="${5:-results_somatic_calling}"

# 5. Directory containing BED files for CNV calling
BED_DIR="${6:-/oak/stanford/groups/emoding/sequencing/pipeline/selectors}"

# 6. Optional: Pre-built pooled normal reference for CNVkit (from build_cnvkit_reference.sh)
#    If provided, CNVkit uses this pooled reference instead of the single matched normal.
#    Build it first:  sbatch build_cnvkit_reference.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED>
POOLED_REF="${7:-}"

# 7. Optional: Mutect2 Panel of Normals VCF (from build_mutect2_pon.sh)
#    Suppresses systematic artifacts and recurrent germline variants during calling.
#    Build it first:  sbatch build_mutect2_pon.sh <NORMAL_SAMPLE_FILE> <NORMAL_DIR> <TARGET_BED>
MUTECT2_PON="${8:-}"

# Path to the reference genome FASTA used to align the BAMs
REF_GENOME="/oak/stanford/groups/emoding/sequencing/pipeline/indices/hg19.fa"

# ==============================================================================
# Sample Selection via SLURM Array
# ==============================================================================
# Ensure SLURM_ARRAY_TASK_ID is set (this script should be run via sbatch --array)
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Please submit using: sbatch --array=1-N script.sh"
    exit 1
fi

# Extract the sample name for this specific array task.
# Select the N-th DATA row (non-empty, non-comment), where N = SLURM_ARRAY_TASK_ID.
# This MUST match how run_full_pipeline.sh computes the array size
# (awk 'NF && $1 !~ /^#/'), otherwise a blank/comment line would desync the task
# index from the intended sample and silently skip or mis-pair samples.
SAMPLE_NAME=$(awk -v t="$SLURM_ARRAY_TASK_ID" 'NF && $1 !~ /^#/ {n++; if(n==t){sub(/\r$/,""); print $1; exit}}' "$CFDNA_SAMPLE_FILE")
TARGET_BED_BASENAME=$(awk -v t="$SLURM_ARRAY_TASK_ID" 'NF && $1 !~ /^#/ {n++; if(n==t){sub(/\r$/,""); print $5; exit}}' "$CFDNA_SAMPLE_FILE")

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "ERROR: Blank sample name retrieved for task ID $SLURM_ARRAY_TASK_ID. Check your cfDNA sample file."
    exit 1
fi

# Validate the per-sample target BED basename (column 5). A blank column would
# silently produce TARGET_BED="${BED_DIR}/" and fail later in a confusing way.
if [[ -z "$TARGET_BED_BASENAME" ]]; then
    echo "ERROR: No target BED basename (column 5) for sample '$SAMPLE_NAME' (task $SLURM_ARRAY_TASK_ID) in $CFDNA_SAMPLE_FILE."
    exit 1
fi

# Convert Tumor sample ID to Patient ID (e.g. LP09-T1 -> LP09, LP09-T2 -> LP09)
PATIENT_ID="${SAMPLE_NAME%-T*}"
TUMOR_NAME="${SAMPLE_NAME}"

# Guard: if the sample name has no '-T' tumor suffix, the parameter expansion above
# leaves PATIENT_ID == SAMPLE_NAME, which will likely fail the matched-normal lookup.
if [[ "$PATIENT_ID" == "$SAMPLE_NAME" ]]; then
    echo "WARNING: Sample '$SAMPLE_NAME' has no '-T' tumor suffix; using the full name as patient ID for normal lookup."
fi

# Dynamically find the matched normal sample from the Normal Sample File
# Searches for the PATIENT_ID followed by a dash (e.g. "LP09-")
NORMAL_NAME=$(awk -v pat="^${PATIENT_ID}-" '$1 ~ pat {sub(/\r$/,""); print $1; exit}' "$NORMAL_SAMPLE_FILE")

if [[ -z "$NORMAL_NAME" ]]; then
    echo "ERROR: Could not find a matching Normal sample for patient $PATIENT_ID in $NORMAL_SAMPLE_FILE!"
    exit 1
fi

echo "Running somatic calling for patient $PATIENT_ID | Tumor: $TUMOR_NAME | Matched Normal: $NORMAL_NAME"

# Build full paths to the matched BAMs
TUMOR_BAM="${CFDNA_DIR}/${CFDNA_PREFIX}${TUMOR_NAME}${CFDNA_SUFFIX}"
NORMAL_BAM="${NORMAL_DIR}/${NORMAL_PREFIX}${NORMAL_NAME}${NORMAL_SUFFIX}"
TARGET_BED="${BED_DIR}/${TARGET_BED_BASENAME}"

# Output workspace directories tailored to this sample
MANTA_RUNDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/manta_somatic_run"
MUTECT_OUTDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/mutect2_somatic_run"
CNVKIT_OUTDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/cnvkit_run"

# ==============================================================================
# Pre-run Checks & Directory Setup
# ==============================================================================
# Ensure output directories exist, and create log directory for SLURM logs
mkdir -p logs "$MANTA_RUNDIR" "$MUTECT_OUTDIR" "$CNVKIT_OUTDIR"

# Verify all input files exist before starting heavy compute
for FILE in "$TUMOR_BAM" "$NORMAL_BAM" "$REF_GENOME" "$TARGET_BED"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: Input file $FILE does not exist!"
        exit 1
    fi
done

# ==============================================================================
# 1. Structural Variant (SV) Calling with Manta (Somatic Mode)
# ==============================================================================
MANTA_SV_VCF="$MANTA_RUNDIR/results/variants/somaticSV.vcf.gz"

if [[ -f "$MANTA_SV_VCF" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Manta already completed for $SAMPLE_NAME. Skipping."
else
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Manta SV Calling for $SAMPLE_NAME..."

    # Step 1.1: Configure the Manta workflow for somatic calling.
    ${HOME}/manta-1.6.0.centos6_x86_64/bin/configManta.py \
        --normalBam "$NORMAL_BAM" \
        --tumorBam "$TUMOR_BAM" \
        --referenceFasta "$REF_GENOME" \
        --runDir "$MANTA_RUNDIR"

    # Step 1.2: Execute the generated runWorkflow.py script.
    # We run it locally on this job's node, parallelizing across the allocated CPUS.
    "$MANTA_RUNDIR/runWorkflow.py" \
        -m local \
        -j "${SLURM_CPUS_PER_TASK:-8}"

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Manta SV Calling Complete for $SAMPLE_NAME."
fi

# ==============================================================================
# 2. SNV and Indel Calling with Mutect2
# ==============================================================================
UNFILTERED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_unfiltered.vcf.gz"
FILTERED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_filtered.vcf.gz"

# NOTE: Mutect2 requires the -normal argument to exactly match the Read Group Sample Name (SM tag) 
# present in the normal BAM file. We dynamically extract it to prevent crashes.
# We temporarily disable pipefail because grep -m1 closes the pipe early,
# causing samtools to throw a SIGPIPE error (which silently kills the whole script).
# Parse the SM tag with tab-aware awk instead of `sed 's/...[^\t]*.../'`, whose
# `\t` is not portably a tab and can mis-capture/truncate the SM value.
# awk exits after the first @RG line, so samtools may receive SIGPIPE — keep pipefail off here.
extract_sm_tag() {
    samtools view -H "$1" | awk -F'\t' '$1=="@RG"{for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i; exit}}'
}
set +o pipefail
NORMAL_SM=$(extract_sm_tag "$NORMAL_BAM")
TUMOR_SM=$(extract_sm_tag "$TUMOR_BAM")
set -o pipefail

# If the normal BAM has no RG, we create a temporary patched BAM
if [[ -z "$NORMAL_SM" ]]; then
    echo "WARNING: Normal BAM is missing @RG lines. Generating a Read Group BAM on the fly..."
    NORMAL_SM="${NORMAL_NAME}"
    PATCHED_NORMAL="${MUTECT_OUTDIR}/${NORMAL_NAME}_rg_patched.bam"
    
    if [[ ! -f "${PATCHED_NORMAL}.bai" ]]; then
        samtools addreplacerg -r "@RG\tID:${NORMAL_NAME}\tSM:${NORMAL_NAME}\tLB:Normal\tPL:ILLUMINA" -o "$PATCHED_NORMAL" -@ 4 "$NORMAL_BAM"
        samtools index -@ 4 "$PATCHED_NORMAL"
    fi
    NORMAL_BAM="$PATCHED_NORMAL"
fi

# If the tumor BAM has no RG, we create a temporary patched BAM
if [[ -z "$TUMOR_SM" ]]; then
    echo "WARNING: Tumor BAM is missing @RG lines. Generating a Read Group BAM on the fly..."
    TUMOR_SM="${SAMPLE_NAME}"
    PATCHED_TUMOR="${MUTECT_OUTDIR}/${SAMPLE_NAME}_rg_patched.bam"
    
    if [[ ! -f "${PATCHED_TUMOR}.bai" ]]; then
        samtools addreplacerg -r "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tLB:cfDNA\tPL:ILLUMINA" -o "$PATCHED_TUMOR" -@ 4 "$TUMOR_BAM"
        samtools index -@ 4 "$PATCHED_TUMOR"
    fi
    TUMOR_BAM="$PATCHED_TUMOR"
fi

# Persist the FINAL tumor/normal SM tags (post-patching) as a sidecar so the downstream
# summarizer can match VCF sample columns by exact SM name instead of guessing from the
# directory name. summarize_variant_calls.sh reads this to avoid collapsing tumor/normal
# onto the same VCF column (which would corrupt Normal_AF/Normal_DP).
if [[ "$TUMOR_SM" == "$NORMAL_SM" ]]; then
    echo "ERROR: Tumor and Normal SM tags are identical ('$TUMOR_SM'). Mutect2 cannot distinguish them."
    echo "       Re-tag the BAMs with distinct SM values."
    exit 1
fi
printf 'TUMOR_SM\t%s\nNORMAL_SM\t%s\n' "$TUMOR_SM" "$NORMAL_SM" > "$MUTECT_OUTDIR/${SAMPLE_NAME}_sm_tags.tsv"

# Skip guard: require BOTH filtered VCF AND orientation model to exist.
# If either is missing (e.g., old run without orientation model), re-run everything.
OB_MODEL="$MUTECT_OUTDIR/${SAMPLE_NAME}_read-orientation-model.tar.gz"

if [[ -f "$FILTERED_VCF" && -f "$OB_MODEL" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 already completed (with orientation model) for $SAMPLE_NAME. Skipping."
else
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Mutect2 SNV/Indel Calling for $SAMPLE_NAME..."

    # Step 2.1: Run GATK Mutect2.
    # cfDNA setting: --minimum-allele-fraction sets the LOWER bound on allele fractions
    # Mutect2 considers. GATK's default is 0.0 (no floor); we set 0.001 (0.1%) as an
    # explicit noise floor that sits well below the 0.5% genotyping cutoff applied in test.R.
    # --max-mnp-distance 0 keeps adjacent SNVs as separate records (NOT merged into MNPs),
    # matching how the PoN was built (build_mutect2_pon.sh) so PoN/COSMIC coordinate joins line up.
    # -L restricts calling to BLYMv2 panel targets (matches PoN intervals).
    # If a Panel of Normals is available, it suppresses systematic sequencing artifacts.
    PON_ARGS=""
    if [[ -n "$MUTECT2_PON" ]]; then
        if [[ -f "$MUTECT2_PON" ]]; then
            echo "  Using Panel of Normals: $MUTECT2_PON"
            PON_ARGS="--panel-of-normals $MUTECT2_PON"
        else
            # Fail closed: a PoN path was requested but is missing. Silently calling
            # without it would drop a clinical artifact filter.
            echo "ERROR: Panel of Normals path provided but not found: $MUTECT2_PON"
            exit 1
        fi
    fi

    gatk Mutect2 \
        -R "$REF_GENOME" \
        -I "$TUMOR_BAM" \
        -I "$NORMAL_BAM" \
        -normal "$NORMAL_SM" \
        -L "$TARGET_BED" \
        -O "$UNFILTERED_VCF" \
        --minimum-allele-fraction 0.001 \
        --max-mnp-distance 0 \
        --f1r2-tar-gz "$MUTECT_OUTDIR/${SAMPLE_NAME}_f1r2.tar.gz" \
        $PON_ARGS

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Calling Complete. Starting LearnReadOrientationModel..."

    # Step 2.1b: Learn Read Orientation Model for oxoG artifacts
    gatk LearnReadOrientationModel \
        -I "$MUTECT_OUTDIR/${SAMPLE_NAME}_f1r2.tar.gz" \
        -O "$MUTECT_OUTDIR/${SAMPLE_NAME}_read-orientation-model.tar.gz"

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Filtering..."

    # Step 2.2: Apply standard somatic filters with FilterMutectCalls.
    gatk FilterMutectCalls \
        -R "$REF_GENOME" \
        -V "$UNFILTERED_VCF" \
        --ob-priors "$MUTECT_OUTDIR/${SAMPLE_NAME}_read-orientation-model.tar.gz" \
        -O "$FILTERED_VCF"

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Filtering Complete for $SAMPLE_NAME."
fi

# Step 2.3: Annotate variants with SnpEff (gene names, variant types, protein changes).
# This produces the annotated VCF needed for OncoPrint / MAF generation.
ANNOTATED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_annotated.vcf.gz"
NORM_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_filtered.norm.vcf.gz"

# Skip only when the annotated VCF AND its normalized provenance both exist. An annotated
# VCF from an OLD run (built before the bcftools-norm step was added) lacks NORM_VCF, so we
# re-run to guarantee the annotated VCF reflects split multiallelics.
if [[ -f "$ANNOTATED_VCF" && -f "$NORM_VCF" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] SnpEff annotation already completed (normalized) for $SAMPLE_NAME. Skipping."
else
    # Step 2.2b: Normalize BEFORE annotation. `bcftools norm -m -any` splits multiallelic
    # records into one ALT per line and left-aligns/trims indels, decomposing the per-ALT
    # FORMAT fields (AF/AD, Number=A/R) accordingly. This makes the downstream summarizer
    # and R filters — which read the first AF/AD element — tie metrics to the correct ALT.
    # bcftools runs from the Lmod module (before activating the snpeff conda env to avoid
    # PATH shadowing); -Oz writes bgzipped output without needing an external bgzip.
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Normalizing (split multiallelics) before annotation..."
    bcftools norm -m -any -f "$REF_GENOME" -Oz -o "$NORM_VCF" "$FILTERED_VCF"
    "$TABIX_BIN" -p vcf "$NORM_VCF"

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting SnpEff Annotation for $SAMPLE_NAME..."

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate snpeff

    # Verify snpEff is actually present in the activated env before relying on it.
    if ! command -v snpEff &>/dev/null; then
        echo "ERROR: 'snpEff' not found in the activated 'snpeff' conda env."
        conda deactivate
        exit 1
    fi

    # SnpEff annotates every variant with gene, effect, and protein change.
    # -noStats suppresses HTML report; -canon uses canonical transcripts only.
    # Use the module bgzip/tabix captured above (the conda env may not provide them).
    snpEff ann -noStats -canon hg19 "$NORM_VCF" | "$BGZIP_BIN" > "$ANNOTATED_VCF"
    "$TABIX_BIN" -p vcf "$ANNOTATED_VCF"

    conda deactivate

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] SnpEff Annotation Complete for $SAMPLE_NAME."
fi

# ==============================================================================
# 3. CNV Calling with CNVkit
# ==============================================================================
# find_one <dir> <find-args...>: print the single matching path, or nothing if none.
# Exits with an error if MORE THAN ONE file matches, instead of silently taking
# `head -n1` (which is order-dependent and can pair mismatched/stale CNVkit files).
find_one() {
    local _dir="$1"; shift
    [[ -d "$_dir" ]] || return 0
    local _matches _n
    _matches=$(find "$_dir" "$@" 2>/dev/null || true)
    _n=$(printf '%s' "$_matches" | grep -c . || true)
    if [[ "$_n" -gt 1 ]]; then
        echo "ERROR: expected exactly one file in $_dir for [$*] but found $_n:" >&2
        printf '%s\n' "$_matches" >&2
        exit 1
    fi
    [[ "$_n" -eq 1 ]] && printf '%s\n' "$_matches"
    return 0
}

# Completion guard: treat CNVkit as done ONLY when the full expected output set exists
# (base .cns + derived .call.cns + .bintest.cns + .cnr), all sharing one stem. A partial
# prior run (e.g. .cns but no .call.cns) must re-run rather than be skipped.
EXISTING_CNS=$(find_one "$CNVKIT_OUTDIR" -name "*.cns" ! -name "*.call.cns" ! -name "*.bintest.cns")
CNVKIT_COMPLETE=false
if [[ -n "$EXISTING_CNS" ]]; then
    _STEM="${EXISTING_CNS%.cns}"
    if [[ -f "${_STEM}.call.cns" && -f "${_STEM}.bintest.cns" && -f "${_STEM}.cnr" ]]; then
        CNVKIT_COMPLETE=true
    fi
fi

if [[ "$CNVKIT_COMPLETE" == true ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit already completed for $SAMPLE_NAME. Skipping."
else
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting CNVkit Calling for $SAMPLE_NAME..."

    # Load conda and activate the cnvkit environment
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate cnvkit

    # Check if CNVkit is available
    if command -v cnvkit.py &> /dev/null; then
        if [[ -n "$POOLED_REF" && -f "$POOLED_REF" ]]; then
            # Pooled normal reference provided (from build_cnvkit_reference.sh).
            # Uses the cohort-averaged normal baseline for more stable CNV calls.
            echo "  Using pooled normal reference: $POOLED_REF"
            cnvkit.py batch "$TUMOR_BAM" \
                --reference "$POOLED_REF" \
                --output-dir "$CNVKIT_OUTDIR" \
                --drop-low-coverage
        else
            # Fallback: single matched normal (no pooled reference available).
            echo "  Using single matched normal: $NORMAL_BAM"
            cnvkit.py batch "$TUMOR_BAM" \
                --normal "$NORMAL_BAM" \
                --targets "$TARGET_BED" \
                --fasta "$REF_GENOME" \
                --output-dir "$CNVKIT_OUTDIR" \
                --method hybrid \
                --drop-low-coverage
        fi

        # Step 3.0b: Generate integer copy number calls and statistical p-values.
        # cnvkit.py batch produces .cns (segments) and .cnr (bin ratios), but NOT .call.cns.
        # We must explicitly run 'call' to get integer CN, then 'bintest' for p-values.
        # call.cns/bintest.cns are derived from the SAME stem as the base .cns so they
        # always pair correctly; find_one guarantees a single unambiguous base file.
        echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting CNVkit call + bintest..."
        TMP_CNR=$(find_one "$CNVKIT_OUTDIR" -name "*.cnr")
        TMP_CNS=$(find_one "$CNVKIT_OUTDIR" -name "*.cns" ! -name "*.call.cns" ! -name "*.bintest.cns")
        if [[ -n "$TMP_CNR" && -n "$TMP_CNS" ]]; then
            # Step 3.0b-i: Integer copy number calling
            CALL_OUT="${TMP_CNS%.cns}.call.cns"
            cnvkit.py call "$TMP_CNS" -o "$CALL_OUT"

            # Step 3.0b-ii: Binomial test for segment-level p-values
            BINTEST_OUT="${TMP_CNS%.cns}.bintest.cns"
            cnvkit.py bintest "$TMP_CNR" -s "$CALL_OUT" > "$BINTEST_OUT"
        else
            echo "WARNING: Expected .cnr/.cns not found after cnvkit batch for $SAMPLE_NAME;"
            echo "         skipping call + bintest (CNV summary will be empty for this sample)."
        fi
    else
        echo "WARNING: cnvkit.py not found after attempting to load conda env. Skipping CNV calling."
    fi

    # Deactivate gracefully when done
    conda deactivate

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit Calling Complete for $SAMPLE_NAME."
fi

# Step 3.1: Generate gene-level copy number metrics for OncoPrint.
CNVKIT_CNR=$(find_one "$CNVKIT_OUTDIR" -name "*.cnr")
CNVKIT_CALL_CNS=$(find_one "$CNVKIT_OUTDIR" -name "*.call.cns")
GENEMETRICS_TSV="$CNVKIT_OUTDIR/${SAMPLE_NAME}_genemetrics.tsv"

if [[ -f "$GENEMETRICS_TSV" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit genemetrics already completed for $SAMPLE_NAME. Skipping."
elif [[ -n "$CNVKIT_CNR" && -f "$CNVKIT_CNR" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting CNVkit genemetrics for $SAMPLE_NAME..."

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate cnvkit

    if command -v cnvkit.py &> /dev/null; then
        # genemetrics maps segment-level CN ratios to individual genes.
        # -t 0.1 = log2 threshold for calling gain/loss at gene level (lowered for cfDNA).
        # -m 3  = minimum number of probes covering a gene to report it.
        if [[ -n "$CNVKIT_CALL_CNS" && -f "$CNVKIT_CALL_CNS" ]]; then
            cnvkit.py genemetrics "$CNVKIT_CNR" -s "$CNVKIT_CALL_CNS" -t 0.1 -m 3 > "$GENEMETRICS_TSV"
        else
            cnvkit.py genemetrics "$CNVKIT_CNR" -t 0.1 -m 3 > "$GENEMETRICS_TSV"
        fi
    else
        echo "WARNING: cnvkit.py not found. Skipping genemetrics."
    fi

    conda deactivate

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit genemetrics Complete for $SAMPLE_NAME."
fi

# ==============================================================================
# Pipeline Complete
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] All somatic variant calling steps completed successfully for $SAMPLE_NAME!"
