#!/bin/bash
#SBATCH --job-name=CNV
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssu42
#SBATCH --partition=emoding

# ==============================================================================
# CNV Calling Workflow (cfDNA + matched normal)
# Designed for Stanford Sherlock HPC Cluster
# 
# Usage: 
#   1. Set the #SBATCH --array directive above to match your sample count.
#      (e.g., if you have 45 samples, use --array=1-45)
#   2. Submit: sbatch cnv_calling_from_bam.sh
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
# Treat unset variables as an error when substituting
# Fail on any error in a pipeline
set -euo pipefail

# ==============================================================================
# Environment Setup
# ==============================================================================
echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] Loading modules..."

# Load miniconda and activate the cnvkit environment
module load miniconda
eval "$(conda shell.bash hook)"
conda activate cnvkit

# ==============================================================================
# Input Variables
# ==============================================================================
if [[ "$#" -lt 4 ]]; then
    echo "Usage: sbatch --array=1-N cnv_calling_from_bam.sh <CFDNA_SAMPLE_FILE> <NORMAL_SAMPLE_FILE> <CFDNA_DIR> <NORMAL_DIR> [OUTPUT_DIR] [BED_DIR]"
    echo "Example: sbatch --array=1-5 cnv_calling_from_bam.sh sample2barcode_tumor.txt sample2barcode_normal.txt /path/to/cfdna /path/to/normal results /path/to/beds"
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

echo "Running CNV calling for patient $PATIENT_ID | Tumor: $TUMOR_NAME | Matched Normal: $NORMAL_NAME"

# Build full paths to the matched BAMs
TUMOR_BAM="${CFDNA_DIR}/${CFDNA_PREFIX}${TUMOR_NAME}${CFDNA_SUFFIX}"
NORMAL_BAM="${NORMAL_DIR}/${NORMAL_PREFIX}${NORMAL_NAME}${NORMAL_SUFFIX}"
TARGET_BED="${BED_DIR}/${TARGET_BED_BASENAME}"

# Output workspace directory tailored to this sample
CNVKIT_OUTDIR="${OUTPUT_DIR}/${SAMPLE_NAME}/cnvkit_run"

# ==============================================================================
# Pre-run Checks & Directory Setup
# ==============================================================================
# Ensure output directories exist, and create log directory for SLURM logs
mkdir -p logs "$CNVKIT_OUTDIR"

# Verify all input files exist before starting heavy compute
for FILE in "$TUMOR_BAM" "$NORMAL_BAM" "$REF_GENOME" "$TARGET_BED"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: Input file $FILE does not exist!"
        exit 1
    fi
done

# ==============================================================================
# 1. CNV Calling with CNVkit
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
    echo "ERROR: cnvkit.py not found after attempting to load conda env. Skipping CNV calling."
    exit 1
fi

# Deactivate gracefully when done
conda deactivate

echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')] All CNV calling steps completed successfully for $SAMPLE_NAME!"
