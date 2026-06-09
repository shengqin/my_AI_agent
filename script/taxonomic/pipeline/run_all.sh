#!/bin/bash
# =============================================================================
# run_all.sh  --  submit the whole microbiome pipeline with SLURM dependencies
# =============================================================================
# Usage:   bash run_all.sh [maxconcurrent]
#   maxconcurrent : max array tasks running at once (default 8; use 2 for Kaiju)
#
# Submits, in order, with afterok dependencies so each step waits for the last:
#   00 sample sheet (run now, on the login node)
#   01 host removal        (array)
#   02 kraken2             (array, after 01)        03 bracken+merge (after 02)
#   04 kaiju               (array, after 01)        04b kaiju2table  (after 04)
#   05 downstream R        (after 03 AND 04b)
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
cd "$SCRIPT_DIR"
mkdir -p logs

MAXC="${1:-8}"
KAIJU_MAXC=2   # Kaiju loads ~187 GB DB per task -> keep concurrency low

# Show the active cohort up front -- this pipeline is reusable across cohorts via
# environment overrides (see config.sh), so confirm you are running the RIGHT one.
echo "=============================================================="
echo " Microbiome pipeline -- cohort settings"
echo "   TRIM_DIR = $TRIM_DIR"
echo "   RUN_DIR  = $RUN_DIR"
echo "   MANIFEST = ${SAMPLE_MANIFEST:-(none -> glob $TRIM_DIR)}"
echo "   REFS_DIR = $REFS_DIR"
echo "=============================================================="

# Pre-flight: verify modules/tools/DBs/inputs before submitting anything.
# Bypass with RUN_ALL_SKIP_PREFLIGHT=1 if you have already checked.
if [ "${RUN_ALL_SKIP_PREFLIGHT:-0}" != "1" ]; then
    echo ">>> [preflight] checking environment..."
    bash "$SCRIPT_DIR/preflight.sh" || {
        echo "ERROR: preflight failed. Fix the FAILs above, or bypass with" >&2
        echo "       RUN_ALL_SKIP_PREFLIGHT=1 bash run_all.sh" >&2
        exit 1; }
fi

echo ">>> [00] building sample sheet"
bash 00_build_sample_sheet.sh
N=$(wc -l < "$SAMPLE_LIST")
[ "$N" -gt 0 ] || { echo "ERROR: no samples" >&2; exit 1; }
echo ">>> $N samples"

# Export PIPELINE_DIR so each job can find config.sh regardless of SLURM's
# spool dir (sbatch copies the script, so ${BASH_SOURCE} won't point here).
submit() { sbatch --export=ALL,PIPELINE_DIR="$SCRIPT_DIR" "$@" | awk '{print $NF}'; }

j01=$(submit --array=1-"$N"%"$MAXC" 01_host_removal.sbatch)
echo "    01 host removal     : $j01"

# Per-sample dependency (aftercorr): task i of 02/04 waits only for task i of 01,
# so ONE failed host-removal sample no longer cancels classification for the
# whole cohort (plain afterok on the array was all-or-nothing).
j02=$(submit --dependency=aftercorr:"$j01" --array=1-"$N"%"$MAXC" 02_kraken2.sbatch)
echo "    02 kraken2          : $j02  (aftercorr $j01)"

# Aggregation steps depend with afterany: they run once the producer FINISHES
# (any state) and then operate on whatever per-sample outputs exist, reporting
# cohort completeness, instead of being cancelled if a single task failed.
j03=$(submit --dependency=afterany:"$j02" 03_bracken_merge.sbatch)
echo "    03 bracken+merge    : $j03  (afterany $j02)"

j04=$(submit --dependency=aftercorr:"$j01" --array=1-"$N"%"$KAIJU_MAXC" 04_kaiju.sbatch)
echo "    04 kaiju            : $j04  (aftercorr $j01)"

j04b=$(submit --dependency=afterany:"$j04" 04b_kaiju2table.sbatch)
echo "    04b kaiju2table     : $j04b (afterany $j04)"

# 05 needs a usable species BIOM (03 now exits non-zero if it is missing), but
# tolerates a missing/failed Kaiju matrix (R falls back to kraken2_only), so it
# requires 03 to SUCCEED and only waits for 04b to finish.
j05=$(submit --dependency=afterok:"$j03",afterany:"$j04b" 05_downstream.sbatch)
echo "    05 downstream R     : $j05  (afterok $j03 + afterany $j04b)"

echo ">>> Submitted. Track with:  squeue -u $USER"
echo ">>> Final tables -> $D_RESULTS ; figures -> $D_FIGS"
