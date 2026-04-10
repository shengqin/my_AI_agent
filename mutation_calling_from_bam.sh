#!/bin/bash
#SBATCH --job-name=B25
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
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
ml gatk/4.6.0.0
# Manta is called directly via absolute path below (no module needed)

# ==============================================================================
# Input Variables
# ==============================================================================
if [[ "$#" -lt 3 ]]; then
    echo "Usage: sbatch --array=1-N mutation_calling_from_bam.sh <SAMPLE_FILE> <CFDNA_DIR> <NORMAL_DIR> [OUTPUT_DIR] [BED_DIR]"
    echo "Example: sbatch --array=1-5 mutation_calling_from_bam.sh sample2barcodeB23_1.txt /path/to/cfdna /path/to/normal results /path/to/beds"
    exit 1
fi

# 1. Path to the sample_to_barcode file (Sample name expected in 1st, Bed in 5th)
SAMPLE_FILE="$1"

# 2. Directory containing cfDNA BAM files
CFDNA_DIR="$2"
CFDNA_PREFIX="Sample_"
CFDNA_SUFFIX=".dualindex-deduped.sorted.bam"

# 3. Directory containing matched normal BAM files
NORMAL_DIR="$3"
NORMAL_PREFIX="Sample_"
NORMAL_SUFFIX=".dualindex-deduped.sorted.bam"

# 4. Output Directory for all variant calling results
OUTPUT_DIR="${4:-results_somatic_calling}"

# 5. Directory containing BED files for CNV calling
BED_DIR="${5:-/oak/stanford/groups/emoding/sequencing/pipeline/selectors}"

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
# We read the N-th line of the sample file, where N = SLURM_ARRAY_TASK_ID
# awk '{print $1}' gets the first column (handling tabs/spaces)
SAMPLE_NAME=$(awk -v task_id="$SLURM_ARRAY_TASK_ID" 'NR==task_id {print $1}' "$SAMPLE_FILE")
TARGET_BED_BASENAME=$(awk -v task_id="$SLURM_ARRAY_TASK_ID" 'NR==task_id {print $5}' "$SAMPLE_FILE")

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "ERROR: Blank sample name retrieved for task ID $SLURM_ARRAY_TASK_ID. Check your sample file."
    exit 1
fi

# Convert Tumor sample ID to Normal sample ID (e.g. LP09-T1 -> LP09-N1)
PATIENT_ID="${SAMPLE_NAME%-T1}"
NORMAL_NAME="${PATIENT_ID}-N1_Normal"
TUMOR_NAME="${SAMPLE_NAME}_cfDNA"

echo "Running somatic calling for sample: $SAMPLE_NAME (Matched Normal: $NORMAL_NAME)"

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

# ==============================================================================
# 2. SNV and Indel Calling with Mutect2
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting Mutect2 SNV/Indel Calling for $SAMPLE_NAME..."

UNFILTERED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_unfiltered.vcf.gz"
FILTERED_VCF="$MUTECT_OUTDIR/${SAMPLE_NAME}_mutect2_filtered.vcf.gz"

# NOTE: Mutect2 requires the -normal argument to exactly match the Read Group Sample Name (SM tag) 
# present in the normal BAM file. We dynamically extract it to prevent crashes.
NORMAL_SM=$(samtools view -H "$NORMAL_BAM" | grep -m1 '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/')

if [[ -z "$NORMAL_SM" ]]; then
    echo "WARNING: Could not extract @RG SM tag from normal BAM. Falling back to filename assumption."
    NORMAL_SM="${NORMAL_NAME}"
fi

# Step 2.1: Run GATK Mutect2.
# Extremely important for cfDNA: --minimum-allele-fraction is manually lowered to 0.001 (0.1%).
gatk Mutect2 \
    -R "$REF_GENOME" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal "$NORMAL_SM" \
    -O "$UNFILTERED_VCF" \
    --minimum-allele-fraction 0.001

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Calling Complete. Starting Filtering..."

# Step 2.2: Apply standard somatic filters with FilterMutectCalls.
gatk FilterMutectCalls \
    -R "$REF_GENOME" \
    -V "$UNFILTERED_VCF" \
    -O "$FILTERED_VCF"

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Mutect2 Filtering Complete for $SAMPLE_NAME."

# ==============================================================================
# 3. CNV Calling with CNVkit
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Starting CNVkit Calling for $SAMPLE_NAME..."

# Check if CNVkit is available
if command -v cnvkit.py &> /dev/null; then
    cnvkit.py batch "$TUMOR_BAM" \
        --normal "$NORMAL_BAM" \
        --targets "$TARGET_BED" \
        --fasta "$REF_GENOME" \
        --output-dir "$CNVKIT_OUTDIR" \
        --drop-low-coverage
else
    echo "WARNING: cnvkit.py not found in PATH. Skipping CNV calling."
    # If it fails, ensure CNVkit is properly installed/loaded via module system if needed
fi

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] CNVkit execution logic finished for $SAMPLE_NAME."

# ==============================================================================
# Pipeline Complete
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] All somatic variant calling steps completed successfully for $SAMPLE_NAME!"
