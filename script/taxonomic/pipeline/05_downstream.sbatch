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
# Kaiju read floor for entering cross-validation (lower than MIN_READS by design)
KJ_MIN="${KAIJU_MIN_READS:-5}"

# --- Provenance manifest -----------------------------------------------------
# Record parameters, tool/DB versions, and the git commit so any run is auditable
# later (important given the FP-control / retraction framing). Runs in a subshell
# with errexit/nounset/pipefail OFF and a trailing `|| true`, so a missing tool or
# unreadable DB can never abort the actual analysis.
MANIFEST="$D_RESULTS/run_manifest.txt"
(
  set +e +u +o pipefail
  echo "NLPHL tumor-microbiome pipeline -- run manifest"
  echo "================================================"
  echo "generated    : $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "host         : $(hostname)"
  echo "user         : ${USER:-n/a}"
  echo "slurm job id : ${SLURM_JOB_ID:-n/a}"
  echo "pipeline dir : $SCRIPT_DIR"
  gitref=$(git -C "$SCRIPT_DIR" rev-parse --short HEAD 2>/dev/null)
  git -C "$SCRIPT_DIR" diff --quiet 2>/dev/null || gitref="${gitref:-n/a} (uncommitted changes present)"
  echo "git commit   : ${gitref:-n/a}"
  echo ""
  echo "--- Parameters (config.sh) ---"
  echo "STAR_TWOPASS=$STAR_TWOPASS  STAR_FFPE_RELAXED=$STAR_FFPE_RELAXED  STAR_SEED_LMAX='${STAR_SEED_LMAX}'"
  echo "LOWCOMPLEXITY_ENTROPY=$LOWCOMPLEXITY_ENTROPY  DEDUP=$DEDUP"
  echo "KRAKEN2_CONFIDENCE=$KRAKEN2_CONFIDENCE  KRAKEN2_MIN_HIT_GROUPS=$KRAKEN2_MIN_HIT_GROUPS  REPORT_MINIMIZER_DATA=$REPORT_MINIMIZER_DATA"
  echo "READ_LEN=$READ_LEN  BRACKEN_THRESH=$BRACKEN_THRESH  BRACKEN_LEVELS='$BRACKEN_LEVELS'"
  echo "KAIJU_MODE=$KAIJU_MODE  KAIJU_MISMATCHES=$KAIJU_MISMATCHES"
  echo "MIN_READS=$MIN_READS  MIN_RPM=$MIN_RPM  KAIJU_MIN_READS=$KAIJU_MIN_READS"
  echo "MIN_PREVALENCE_FRAC=$MIN_PREVALENCE_FRAC  MIN_MINIMIZER_COVERAGE=$MIN_MINIMIZER_COVERAGE  MIN_DISTINCT_MINIMIZERS=$MIN_DISTINCT_MINIMIZERS"
  echo "REQUIRE_BOTH_CLASSIFIERS=$REQUIRE_BOTH_CLASSIFIERS  CONTROLS_LIST='${CONTROLS_LIST:-}'"
  echo "CONTAMINANT_BLOCKLIST=$CONTAMINANT_BLOCKLIST"
  echo ""
  echo "--- Tool versions ---"
  echo "R           : $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
  echo "STAR        : $("$STAR_BIN_DIR/STAR" --version 2>/dev/null)"
  echo "kraken2     : $("$CONTAINER_RUN" exec "$KRAKEN2_SIF" kraken2 --version 2>/dev/null | head -1)"
  echo "bracken     : $("$CONTAINER_RUN" exec "$BRACKEN_SIF" bracken --version 2>/dev/null | head -1)"
  echo "kaiju       : $("$KAIJU_BIN/kaiju" -h 2>&1 | grep -im1 kaiju)"
  echo "kraken-biom : $(command -v kraken-biom >/dev/null 2>&1 && kraken-biom --version 2>/dev/null || echo n/a)"
  echo ""
  echo "--- Databases (path | size | mtime) ---"
  echo "Kraken2 DB : $KRAKEN2_DB"
  for f in hash.k2d taxo.k2d opts.k2d; do
    p="$KRAKEN2_DB/$f"
    [ -e "$p" ] && echo "  $f | $(du -h "$p" | cut -f1) | $(date -r "$p" '+%Y-%m-%d')" || echo "  $f | MISSING"
  done
  echo "Kaiju idx  : $KAIJU_INDEX"
  [ -e "$KAIJU_INDEX" ] && echo "  $(du -h "$KAIJU_INDEX" | cut -f1) | $(date -r "$KAIJU_INDEX" '+%Y-%m-%d')" || echo "  MISSING"
  echo "STAR index : $STAR_GENOME_DIR/SAindex"
  [ -e "$STAR_GENOME_DIR/SAindex" ] && echo "  $(du -h "$STAR_GENOME_DIR/SAindex" | cut -f1) | $(date -r "$STAR_GENOME_DIR/SAindex" '+%Y-%m-%d')" || echo "  MISSING"
) > "$MANIFEST" 2>/dev/null || true
echo ">>> Wrote provenance manifest: $MANIFEST"

echo ">>> Aggregating + filtering + cross-validating"
Rscript "$SCRIPT_DIR/05_aggregate_filter.R" \
    "$D_RESULTS" "$BIOM" "$KAIJU_MAT" "$MMZ" "$TOTALS" \
    "$D_HOST/readcounts" "$CONTAMINANT_BLOCKLIST" \
    "$MIN_READS" "$MIN_RPM" "$MIN_PREVALENCE_FRAC" \
    "$MIN_MINIMIZER_COVERAGE" "$MIN_DISTINCT" \
    "$REQUIRE_BOTH_CLASSIFIERS" "$CTRL" "$KJ_MIN"

echo ">>> Visualizing"
Rscript "$SCRIPT_DIR/06_visualize.R" "$D_RESULTS" "$D_FIGS"

echo ">>> STEP 05 complete. Tables in $D_RESULTS, figures in $D_FIGS"
