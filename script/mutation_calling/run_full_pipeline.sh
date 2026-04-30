#!/bin/bash
# ==============================================================================
# Master Orchestrator: Full CAPP-Seq Somatic Variant Calling Pipeline
#
# Submits all pipeline steps in the correct order with SLURM dependencies:
#   Step 1: Build CNVkit pooled normal reference  (single job)
#   Step 2: Build Mutect2 Panel of Normals        (single job, parallel with Step 1)
#   Step 3: Run mutation calling per batch         (job arrays, after Steps 1+2)
#   Step 4: Summarize all results                  (single job, after Step 3)
#
# Usage:
#   bash run_full_pipeline.sh                  # Full run (build PoN + calling)
#   SKIP_PON_BUILD=true bash run_full_pipeline.sh  # Skip Steps 1+2, reuse existing PoNs
#
# All paths are configured below. Edit the BATCHES array to add/remove batches.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# Runtime Flags — Override at the command line, e.g.:
#   SKIP_PON_BUILD=true bash run_full_pipeline.sh
# ==============================================================================
SKIP_PON_BUILD="${SKIP_PON_BUILD:-false}"

# ==============================================================================
# Configuration — Edit this section for your cohort
# ==============================================================================
PIPELINE_DIR="/oak/stanford/groups/emoding/sequencing/pipeline/runs/cappseq/BinkleyLab"
SCRIPT_DIR="${PIPELINE_DIR}"   # Where mutation_calling_from_bam.sh lives

# Shared normal samples (all batches use the same normals from B25)
NORMAL_SAMPLE_FILE="${PIPELINE_DIR}/NovaSeqB25/sample2barcodeB25.txt"
NORMAL_DIR="${PIPELINE_DIR}/NovaSeqB25/demultiplexed/barcode-deduped/normal"

# Target BED file for CNVkit and PoN
TARGET_BED="/oak/stanford/groups/emoding/sequencing/pipeline/selectors/design_BLYMv2_20210721_hg19_sorted.bed"
BED_DIR="/oak/stanford/groups/emoding/sequencing/pipeline/selectors"

# Output directories for pooled references
CNVKIT_REF_DIR="${PIPELINE_DIR}/cnvkit_pooled_ref"
MUTECT2_PON_DIR="${PIPELINE_DIR}/mutect2_pon"

# Summary output
SUMMARY_DIR="${PIPELINE_DIR}/combined_summary"

# SLURM defaults for per-sample mutation calling jobs
SLURM_PARTITION="emoding"
SLURM_TIME="1-00:00:00"
SLURM_CPUS=6
SLURM_MEM_PER_CPU="6G"
SLURM_MAIL="ssu42"

# ==============================================================================
# Batch Definitions — One entry per sequencing batch
#
# Format: "BATCH_NAME|SAMPLE_FILE|CFDNA_DIR|OUTPUT_DIR|ARRAY_SIZE"
# ==============================================================================
BATCHES=(
    "B23_1|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_1/sample2barcodeB23_1.txt|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_1/demultiplexed/barcode-deduped/cfdna|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_1/results_somatic_calling|5"
    "B23_2|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_2/sample2barcodeB23_2.txt|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_2/demultiplexed/barcode-deduped/cfdna|${PIPELINE_DIR}/NovaSeqB23/NovaSeqB23_2/results_somatic_calling|6"
    "B24_1|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_1/sample2barcodeB24_1.txt|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_1/demultiplexed/barcode-deduped/cfdna|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_1/results_somatic_calling|5"
    "B24_2|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_2/sample2barcodeB24_2.txt|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_2/demultiplexed/barcode-deduped/cfdna|${PIPELINE_DIR}/NovaSeqB24/NovaSeqB24_2/results_somatic_calling|6"
)

echo "============================================="
echo "  CAPP-Seq Pipeline Orchestrator"
echo "============================================="
echo "  Batches: ${#BATCHES[@]}"
echo "  Normal source: NovaSeqB25"
echo "  SKIP_PON_BUILD: ${SKIP_PON_BUILD}"
echo ""

# ==============================================================================
# Pre-flight: Rename existing output directories to avoid stale results
# Appends current date (e.g., results_somatic_calling_2026-04-29)
# ==============================================================================
DATE_STAMP=$(date +%Y-%m-%d)

rename_if_exists() {
    local DIR="$1"
    if [[ -d "$DIR" ]]; then
        local BACKUP="${DIR}_${DATE_STAMP}"
        # If backup already exists today, add time
        if [[ -d "$BACKUP" ]]; then
            BACKUP="${DIR}_$(date +%Y-%m-%d_%H%M%S)"
        fi
        echo "  Renaming: $(basename $DIR) -> $(basename $BACKUP)"
        mv "$DIR" "$BACKUP"
    fi
}

echo "[Pre-flight] Checking for existing output directories..."

# Rename per-batch results directories
for BATCH_ENTRY in "${BATCHES[@]}"; do
    IFS='|' read -r BATCH_NAME SAMPLE_FILE CFDNA_DIR OUTPUT_DIR ARRAY_SIZE <<< "$BATCH_ENTRY"
    rename_if_exists "$OUTPUT_DIR"
done

# Only rename reference dirs if we are rebuilding them
if [[ "${SKIP_PON_BUILD}" != "true" ]]; then
    rename_if_exists "$CNVKIT_REF_DIR"
    rename_if_exists "$MUTECT2_PON_DIR"
fi

rename_if_exists "$SUMMARY_DIR"

echo ""

# ==============================================================================
# Step 1: Build CNVkit Pooled Normal Reference
# ==============================================================================
if [[ "${SKIP_PON_BUILD}" == "true" ]]; then
    echo "[Step 1] SKIPPED — Reusing existing CNVkit reference: ${CNVKIT_REF_DIR}"
    # Verify the expected output file exists before proceeding
    if [[ ! -f "${CNVKIT_REF_DIR}/pooled_reference.cnn" ]]; then
        echo "ERROR: Expected reference file not found: ${CNVKIT_REF_DIR}/pooled_reference.cnn"
        echo "       Re-run without SKIP_PON_BUILD=true to rebuild."
        exit 1
    fi
    CNVKIT_JOB=""
else
    echo "[Step 1] Submitting CNVkit pooled reference build..."
    CNVKIT_JOB=$(sbatch \
        --job-name=cnvkit_ref \
        --time=1-00:00:00 \
        --ntasks=1 \
        --cpus-per-task=8 \
        --mem-per-cpu=8G \
        --mail-type=ALL \
        --mail-user=${SLURM_MAIL} \
        --partition=${SLURM_PARTITION} \
        --parsable \
        ${SCRIPT_DIR}/build_cnvkit_reference.sh \
        "${NORMAL_SAMPLE_FILE}" \
        "${NORMAL_DIR}" \
        "${TARGET_BED}" \
        "${CNVKIT_REF_DIR}")
    echo "  CNVkit reference job: ${CNVKIT_JOB}"
fi

# ==============================================================================
# Step 2: Build Mutect2 Panel of Normals
# ==============================================================================
if [[ "${SKIP_PON_BUILD}" == "true" ]]; then
    echo "[Step 2] SKIPPED — Reusing existing Mutect2 PoN: ${MUTECT2_PON_DIR}"
    if [[ ! -f "${MUTECT2_PON_DIR}/pon.vcf.gz" ]]; then
        echo "ERROR: Expected PoN file not found: ${MUTECT2_PON_DIR}/pon.vcf.gz"
        echo "       Re-run without SKIP_PON_BUILD=true to rebuild."
        exit 1
    fi
    PON_JOB=""
else
    echo "[Step 2] Submitting Mutect2 PoN build..."
    PON_JOB=$(sbatch \
        --job-name=mutect2_pon \
        --time=2-00:00:00 \
        --ntasks=1 \
        --cpus-per-task=8 \
        --mem-per-cpu=8G \
        --mail-type=ALL \
        --mail-user=${SLURM_MAIL} \
        --partition=${SLURM_PARTITION} \
        --parsable \
        ${SCRIPT_DIR}/build_mutect2_pon.sh \
        "${NORMAL_SAMPLE_FILE}" \
        "${NORMAL_DIR}" \
        "${TARGET_BED}" \
        "${MUTECT2_PON_DIR}")
    echo "  Mutect2 PoN job: ${PON_JOB}"
fi

# ==============================================================================
# Step 3: Submit Per-Batch Mutation Calling (waits for Steps 1+2 if submitted)
# ==============================================================================
echo "[Step 3] Submitting per-batch mutation calling..."

CALLING_JOBS=""

for BATCH_ENTRY in "${BATCHES[@]}"; do
    IFS='|' read -r BATCH_NAME SAMPLE_FILE CFDNA_DIR OUTPUT_DIR ARRAY_SIZE <<< "$BATCH_ENTRY"

    # Build the dependency string only if PoN jobs were submitted
    DEPENDENCY_ARG=""
    if [[ -n "${CNVKIT_JOB}" && -n "${PON_JOB}" ]]; then
        DEPENDENCY_ARG="--dependency=afterok:${CNVKIT_JOB}:${PON_JOB}"
    fi

    BATCH_JOB=$(sbatch \
        --job-name=${BATCH_NAME} \
        --time=${SLURM_TIME} \
        --ntasks=1 \
        --cpus-per-task=${SLURM_CPUS} \
        --mem-per-cpu=${SLURM_MEM_PER_CPU} \
        --mail-type=ALL \
        --mail-user=${SLURM_MAIL} \
        --partition=${SLURM_PARTITION} \
        --array=1-${ARRAY_SIZE} \
        ${DEPENDENCY_ARG} \
        --parsable \
        ${SCRIPT_DIR}/mutation_calling_from_bam.sh \
        "${SAMPLE_FILE}" \
        "${NORMAL_SAMPLE_FILE}" \
        "${CFDNA_DIR}" \
        "${NORMAL_DIR}" \
        "${OUTPUT_DIR}" \
        "${BED_DIR}" \
        "${CNVKIT_REF_DIR}/pooled_reference.cnn" \
        "${MUTECT2_PON_DIR}/pon.vcf.gz")

    echo "  ${BATCH_NAME}: job ${BATCH_JOB} (array 1-${ARRAY_SIZE})"

    # Collect all batch job IDs for the summary dependency
    if [[ -n "$CALLING_JOBS" ]]; then
        CALLING_JOBS="${CALLING_JOBS}:${BATCH_JOB}"
    else
        CALLING_JOBS="${BATCH_JOB}"
    fi
done

# ==============================================================================
# Step 4: Summarize All Results (waits for all Step 3 jobs)
# ==============================================================================
echo "[Step 4] Submitting variant summary (after all calling jobs)..."

# Build the list of results directories from batch definitions
RESULTS_DIRS_ARGS=""
for BATCH_ENTRY in "${BATCHES[@]}"; do
    IFS='|' read -r _ _ _ OUTPUT_DIR _ <<< "$BATCH_ENTRY"
    RESULTS_DIRS_ARGS="${RESULTS_DIRS_ARGS} ${OUTPUT_DIR}"
done

SUMMARY_JOB=$(sbatch \
    --job-name=summarize \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=2 \
    --mem-per-cpu=4G \
    --mail-type=ALL \
    --mail-user=${SLURM_MAIL} \
    --partition=${SLURM_PARTITION} \
    --dependency=afterok:${CALLING_JOBS} \
    --parsable \
    ${SCRIPT_DIR}/summarize_variant_calls.sh \
    "${SUMMARY_DIR}" \
    ${RESULTS_DIRS_ARGS})

echo "  Summary job: ${SUMMARY_JOB}"

# ==============================================================================
# Done — Print Dependency Chain
# ==============================================================================
echo ""
echo "============================================="
echo "  All Jobs Submitted!"
echo "============================================="
echo ""
if [[ "${SKIP_PON_BUILD}" == "true" ]]; then
    echo "  Dependency chain (PoN steps skipped):"
    echo "    Step 3: Calling   [${CALLING_JOBS}]"
    echo "                           ▼"
    echo "    Step 4: Summarize [${SUMMARY_JOB}]"
else
    echo "  Dependency chain:"
    echo "    Step 1: CNVkit ref  [${CNVKIT_JOB}]  ──┐"
    echo "    Step 2: Mutect2 PoN [${PON_JOB}]  ──┤"
    echo "                                          ▼"
    echo "    Step 3: Calling     [${CALLING_JOBS}]"
    echo "                                          ▼"
    echo "    Step 4: Summarize   [${SUMMARY_JOB}]"
fi
echo ""
echo "  Monitor with: squeue -u \$USER"
echo "  Cancel all:   scancel ${CALLING_JOBS} ${SUMMARY_JOB}"
echo "============================================="
