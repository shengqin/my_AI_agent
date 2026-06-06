#!/bin/bash
#SBATCH --job-name=mb_downstream
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/05_downstream_%j.out
#SBATCH --error=logs/05_downstream_%j.err
# -----------------------------------------------------------------------------
# STEP 05 -- aggregate + filter + cross-validate (R), then visualize (R)
# -----------------------------------------------------------------------------
#     sbatch 05_downstream.sbatch     # after steps 03 and 04b have finished
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
mkdir -p "$D_RESULTS" "$D_FIGS" logs

# --- R environment (R 4.4.2; do NOT 'module --force purge') ---
module load R/4.4.2 2>/dev/null || true
command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found after 'module load R/4.4.2'." >&2; exit 1; }
echo ">>> R: $(Rscript -e 'cat(R.version.string)' 2>/dev/null)  ($(command -v Rscript))"
# Fail fast (with a fix-it pointer) if required packages are missing for THIS R version.
Rscript -e 'req<-c("readr","dplyr","tidyr","tibble","stringr","purrr","ggplot2","scales","phyloseq","biomformat"); m<-req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)]; if(length(m)){cat("MISSING R packages:",paste(m,collapse=", "),"\n"); quit(status=3)}' \
  || { echo "ERROR: install R deps once with:  ml R/4.4.2 && Rscript \"$SCRIPT_DIR/install_r_deps.R\"" >&2; exit 1; }

BIOM="$D_BRACKEN/NLPHL_microbiome_species.biom"
KAIJU_MAT="$D_KAIJU/cohort_full_lineage_matrix.tsv"
MMZ="$D_BRACKEN/minimizer_stats.tsv"
TOTALS="$D_BRACKEN/total_reads_lookup.tsv"
CTRL="${CONTROLS_LIST:-}"; [ -z "$CTRL" ] && CTRL="NONE"

# distinct-minimizer floor (paired with MIN_MINIMIZER_COVERAGE from config)
MIN_DISTINCT="${MIN_DISTINCT_MINIMIZERS:-10}"

echo ">>> Aggregating + filtering + cross-validating"
Rscript "$SCRIPT_DIR/05_aggregate_filter.R" \
    "$D_RESULTS" "$BIOM" "$KAIJU_MAT" "$MMZ" "$TOTALS" \
    "$D_HOST/readcounts" "$CONTAMINANT_BLOCKLIST" \
    "$MIN_READS" "$MIN_RPM" "$MIN_PREVALENCE_FRAC" \
    "$MIN_MINIMIZER_COVERAGE" "$MIN_DISTINCT" \
    "$REQUIRE_BOTH_CLASSIFIERS" "$CTRL"

echo ">>> Visualizing"
Rscript "$SCRIPT_DIR/06_visualize.R" "$D_RESULTS" "$D_FIGS"

echo ">>> STEP 05 complete. Tables in $D_RESULTS, figures in $D_FIGS"
