#!/bin/bash
#SBATCH --job-name=mb_kraken2
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/02_kraken2_%A_%a.out
#SBATCH --error=logs/02_kraken2_%A_%a.err
# -----------------------------------------------------------------------------
# STEP 02 -- Kraken2 nucleotide classification of candidate-microbial reads
# -----------------------------------------------------------------------------
#     sbatch --array=1-<N>%4 02_kraken2.sbatch
# mem must exceed the DB size held in RAM (standard DB ~ 80-140 GB).
#
# KEY DIFFERENCE vs the scratch run:
#   --confidence 0.1 and --minimum-hit-groups 3 suppress the false positives
#   that produced ~99% "classified" on host-depleted reads. We also request
#   --report-minimizer-data so step 03 can apply a KrakenUniq-style
#   distinct-minimizer coverage filter, then derive a plain 6-column report
#   that Bracken can read.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
mkdir -p "$D_KRAKEN" logs

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
sample=$(echo "$LINE" | cut -f2)
THREADS=${SLURM_CPUS_PER_TASK:-8}

r1="$D_HOST/${sample}_microbe_R1.fastq"
r2="$D_HOST/${sample}_microbe_R2.fastq"
[ -s "$r1" ] || { echo "ERROR: missing host-removed reads for $sample ($r1)" >&2; exit 1; }

report_std="$D_KRAKEN/${sample}.kreport"
report_mmz="$D_KRAKEN/${sample}.kreport.minimizer"
# Resume only if the artifacts THIS run needs already exist. When minimizer data
# is requested, the .kreport.minimizer must ALSO be present -- otherwise a rerun
# after an earlier no-minimizer run would skip the sample, step 03 would find no
# minimizer report, and the step-05 distinct-minimizer filter would silently
# become a no-op for that sample.
if [ -s "$report_std" ] && { [ "$REPORT_MINIMIZER_DATA" != "true" ] || [ -s "$report_mmz" ]; }; then
    echo "[SKIP] $sample kraken2 done."; exit 0
fi

echo ">>> Kraken2 on $sample (confidence=$KRAKEN2_CONFIDENCE, hit-groups=$KRAKEN2_MIN_HIT_GROUPS)"

mmz_flag=""
report_target="$report_std"
if [ "$REPORT_MINIMIZER_DATA" = "true" ]; then
    mmz_flag="--report-minimizer-data"
    report_target="$report_mmz"
fi

# --output /dev/null: we don't keep the per-read assignments (huge on big cohorts).
# Change to a real path if you need to extract reads for a specific taxon later.
"$CONTAINER_RUN" exec "$KRAKEN2_SIF" kraken2 \
    --db "$KRAKEN2_DB" \
    --threads "$THREADS" \
    --paired "$r1" "$r2" \
    --confidence "$KRAKEN2_CONFIDENCE" \
    --minimum-hit-groups "$KRAKEN2_MIN_HIT_GROUPS" \
    $mmz_flag \
    --report "$report_target" \
    --output /dev/null

# If we made an 8-column minimizer report, derive the plain 6-column report
# (drop the two inserted minimizer columns) so Bracken can consume it.
if [ "$REPORT_MINIMIZER_DATA" = "true" ]; then
    awk -F'\t' 'BEGIN{OFS="\t"} { if (NF>=8) print $1,$2,$3,$6,$7,$8; else print $0 }' \
        "$report_mmz" > "$report_std"
fi

echo ">>> DONE kraken2 $sample"
