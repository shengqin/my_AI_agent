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
# Input Variables (Modify these paths for your specific run)
# ==============================================================================
# 1. Path to the sample_to_barcode file (Sample name expected in the 1st column)
SAMPLE_FILE="/path/to/sample_to_barcode.txt"

# 2. Directory containing cfDNA BAM files
# We assume the BAMs are named something like: {SAMPLE_NAME}_cfDNA.bam
# Update the suffix if yours are named exactly the same as the sample, or differently
CFDNA_DIR="/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/demultiplexed/barcode-deduped/cfDNA"
CFDNA_PREFIX="Sample_"
CFDNA_SUFFIX="-T1_cfDNA.dualindex-deduped.sorted.bam" # e.g., ".bam", "_cfdna.bam", or "_tumor.bam"

# 3. Directory containing matched normal BAM files
# We assume the BAMs are named something like: Sample_{SAMPLE_NAME}-N1_Normal.dualindex-deduped.sorted.bam
NORMAL_DIR="/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab/NovaSeqB25/demultiplexed/barcode-deduped/normal"
NORMAL_PREFIX="Sample_"
NORMAL_SUFFIX="-N1_Normal.dualindex-deduped.sorted.bam" # Select the specific variant you want to use

# 4. Output Directory for all variant calling results
OUTPUT_DIR="results_somatic_calling"

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

if [[ -z "$SAMPLE_NAME" ]]; then
    echo "ERROR: Blank sample name retrieved for task ID $SLURM_ARRAY_TASK_ID. Check your sample file."
    exit 1
fi

echo "Running somatic calling for sample: $SAMPLE_NAME"

# Build full paths to the matched BAMs
TUMOR_BAM="${CFDNA_DIR}/${CFDNA_PREFIX}${SAMPLE_NAME}${CFDNA_SUFFIX}"
NORMAL_BAM="${NORMAL_DIR}/${NORMAL_PREFIX}${SAMPLE_NAME}${NORMAL_SUFFIX}"

# Output workspace directories tailored to this sample
MANTA_RUNDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/manta_somatic_run"
MUTECT_OUTDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/mutect2_somatic_run"

# ==============================================================================
# Pre-run Checks & Directory Setup
# ==============================================================================
# Ensure output directories exist, and create log directory for SLURM logs
mkdir -p logs "$MANTA_RUNDIR" "$MUTECT_OUTDIR"

# Verify all input files exist before starting heavy compute
for FILE in "$TUMOR_BAM" "$NORMAL_BAM" "$REF_GENOME"; do
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
# present in the normal BAM file. If your BAM's SM tag is NOT equal to $SAMPLE_NAME, standard 
# Mutect2 behavior will fail. You can extract it dynamically if needed:
# NORMAL_SM=$(samtools view -H "$NORMAL_BAM" | grep -m1 '^@RG' | sed 's/.*SM:\([^\t]*\).*/\1/')

# Step 2.1: Run GATK Mutect2.
# Extremely important for cfDNA: --minimum-allele-fraction is manually lowered to 0.001 (0.1%).
gatk Mutect2 \
    -R "$REF_GENOME" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal "$SAMPLE_NAME" \
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
# Pipeline Complete
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] All somatic variant calling steps completed successfully for $SAMPLE_NAME!"
