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
ml python/2.7.13
ml gatk/4.6.0.0
# Manta is called directly via absolute path below (no module needed)

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

# Extract the sample name for this specific array task
# We read the N-th line of the cfDNA sample file, where N = SLURM_ARRAY_TASK_ID
SAMPLE_NAME=$(awk -v task_id="$SLURM_ARRAY_TASK_ID" 'NR==task_id {sub(/\r$/,""); print $1}' "$CFDNA_SAMPLE_FILE")
TARGET_BED_BASENAME=$(awk -v task_id="$SLURM_ARRAY_TASK_ID" 'NR==task_id {sub(/\r$/,""); print $5}' "$CFDNA_SAMPLE_FILE")

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "ERROR: Blank sample name retrieved for task ID $SLURM_ARRAY_TASK_ID. Check your cfDNA sample file."
    exit 1
fi

# Convert Tumor sample ID to Patient ID (e.g. LP09-T1 -> LP09, LP09-T2 -> LP09)
PATIENT_ID="${SAMPLE_NAME%-T*}"
TUMOR_NAME="${SAMPLE_NAME}"

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
set +o pipefail
NORMAL_SM=$(samtools view -H "$NORMAL_BAM" | grep -m1 '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/')
TUMOR_SM=$(samtools view -H "$TUMOR_BAM" | grep -m1 '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/')
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

if [[ -f "$FILTERED_VCF" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 already completed for $SAMPLE_NAME. Skipping."
else
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Mutect2 SNV/Indel Calling for $SAMPLE_NAME..."

    # Step 2.1: Run GATK Mutect2.
    # Extremely important for cfDNA: --minimum-allele-fraction is manually lowered to 0.001 (0.1%).
    # If a Panel of Normals is available, it suppresses systematic sequencing artifacts.
    PON_ARGS=""
    if [[ -n "$MUTECT2_PON" && -f "$MUTECT2_PON" ]]; then
        echo "  Using Panel of Normals: $MUTECT2_PON"
        PON_ARGS="--panel-of-normals $MUTECT2_PON"
    fi

    gatk Mutect2 \
        -R "$REF_GENOME" \
        -I "$TUMOR_BAM" \
        -I "$NORMAL_BAM" \
        -normal "$NORMAL_SM" \
        -O "$UNFILTERED_VCF" \
        --minimum-allele-fraction 0.001 \
        $PON_ARGS

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Calling Complete. Starting Filtering..."

    # Step 2.2: Apply standard somatic filters with FilterMutectCalls.
    gatk FilterMutectCalls \
        -R "$REF_GENOME" \
        -V "$UNFILTERED_VCF" \
        -O "$FILTERED_VCF"

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Filtering Complete for $SAMPLE_NAME."
fi

# Step 2.3: Annotate variants with SnpEff (gene names, variant types, protein changes).
# This produces the annotated VCF needed for OncoPrint / MAF generation.
ANNOTATED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_annotated.vcf.gz"

if [[ -f "$ANNOTATED_VCF" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] SnpEff annotation already completed for $SAMPLE_NAME. Skipping."
else
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting SnpEff Annotation for $SAMPLE_NAME..."

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate snpeff

    # SnpEff annotates every variant with gene, effect, and protein change.
    # -noStats suppresses HTML report; -canon uses canonical transcripts only.
    snpEff ann -noStats -canon hg19 "$FILTERED_VCF" | bgzip > "$ANNOTATED_VCF"
    tabix -p vcf "$ANNOTATED_VCF"

    conda deactivate

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] SnpEff Annotation Complete for $SAMPLE_NAME."
fi

# ==============================================================================
# 3. CNV Calling with CNVkit
# ==============================================================================
CNVKIT_CNS=$(find "$CNVKIT_OUTDIR" -name "*.cns" 2>/dev/null | head -n 1)

if [[ -n "$CNVKIT_CNS" ]]; then
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
    else
        echo "WARNING: cnvkit.py not found after attempting to load conda env. Skipping CNV calling."
    fi

    # Deactivate gracefully when done
    conda deactivate

    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit Calling Complete for $SAMPLE_NAME."
fi

# Step 3.1: Generate gene-level copy number metrics for OncoPrint.
CNVKIT_CNR=$(find "$CNVKIT_OUTDIR" -name "*.cnr" 2>/dev/null | head -n 1)
CNVKIT_CNS=$(find "$CNVKIT_OUTDIR" -name "*.call.cns" 2>/dev/null | head -n 1)
GENEMETRICS_TSV="$CNVKIT_OUTDIR/${SAMPLE_NAME}_genemetrics.tsv"

if [[ -f "$GENEMETRICS_TSV" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit genemetrics already completed for $SAMPLE_NAME. Skipping."
elif [[ -f "$CNVKIT_CNR" ]]; then
    echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting CNVkit genemetrics for $SAMPLE_NAME..."

    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate cnvkit

    if command -v cnvkit.py &> /dev/null; then
        # genemetrics maps segment-level CN ratios to individual genes.
        # -t 0.2 = log2 threshold for calling gain/loss at gene level.
        # -m 3  = minimum number of probes covering a gene to report it.
        if [[ -f "$CNVKIT_CNS" ]]; then
            cnvkit.py genemetrics "$CNVKIT_CNR" -s "$CNVKIT_CNS" -t 0.2 -m 3 > "$GENEMETRICS_TSV"
        else
            cnvkit.py genemetrics "$CNVKIT_CNR" -t 0.2 -m 3 > "$GENEMETRICS_TSV"
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
