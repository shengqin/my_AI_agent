#!/bin/bash
# =============================================================================
# preflight.sh -- verify the environment BEFORE submitting the pipeline.
# =============================================================================
# Runs on a Sherlock login (or compute) node. It does NOT submit any job and
# NEVER purges modules. Each step's modules are loaded in an ISOLATED subshell,
# so the test reflects the real per-job environment (no cross-step conflicts).
#
#     bash preflight.sh
#
# Exit 0 = ready to submit; exit 1 = at least one hard FAIL to resolve.
# (WARN items are optional tools/refs whose absence only skips a sub-step.)
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# Make `module` available even in a non-login shell.
if ! command -v module >/dev/null 2>&1; then
  for i in /etc/profile.d/modules.sh "${LMOD_PKG:-}/init/bash" /usr/share/lmod/lmod/init/bash; do
    [ -f "$i" ] && { source "$i"; break; }
  done
fi
command -v module >/dev/null 2>&1 || { echo "ERROR: 'module' unavailable -- run on a Sherlock node."; exit 1; }

RES=$(mktemp); trap 'rm -f "$RES"' EXIT
ok(){ echo OK >>"$RES"; echo "  [ OK ] $*"; }
no(){ echo NO >>"$RES"; echo "  [FAIL] $*"; }
wn(){ echo WN >>"$RES"; echo "  [WARN] $*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

echo "=== STEP 01 env: STAR + de-hosting tools ==="
(
  module load biology         2>/dev/null
  module load samtools        2>/dev/null
  module load biology bbmap   2>/dev/null
  module load biology bowtie2 2>/dev/null
  export PATH="$STAR_BIN_DIR:$PATH"
  [ -n "${BBMAP_BIN:-}" ] && export PATH="$BBMAP_BIN:$PATH"
  if have STAR; then ok "STAR -> $(STAR --version 2>&1)  ($(command -v STAR))"
  else no "STAR not found (STAR_BIN_DIR=$STAR_BIN_DIR)"; fi
  have samtools    && ok "samtools"    || no "samtools (try: module load biology samtools)"
  have bbduk.sh    && ok "bbduk.sh"    || wn "bbduk.sh absent -> rRNA + entropy sub-steps will skip"
  have clumpify.sh && ok "clumpify.sh" || wn "clumpify.sh absent -> dedup sub-step will skip"
  have bowtie2     && ok "bowtie2"     || wn "bowtie2 absent -> T2T sub-step skips (off by default anyway)"
)

echo "=== STEP 02/03 env: Kraken2/Bracken containers + python ==="
(
  module load python/3.9.0 2>/dev/null
  # Test exactly the runtime the steps will invoke ($CONTAINER_RUN from config.sh),
  # not just "either is present" -- otherwise preflight could pass on singularity
  # while steps 02/03 call the configured runtime and die.
  RUN="$CONTAINER_RUN"
  if have "$RUN"; then ok "container runtime: $RUN (steps use \$CONTAINER_RUN)"
  else no "container runtime '$RUN' not on PATH (try: module load system apptainer)"; RUN=""; fi
  [ -e "$KRAKEN2_SIF" ] && ok "kraken2.sif present" || no "kraken2.sif MISSING: $KRAKEN2_SIF"
  [ -e "$BRACKEN_SIF" ] && ok "bracken.sif present" || no "bracken.sif MISSING: $BRACKEN_SIF"
  if [ -n "$RUN" ] && [ -e "$KRAKEN2_SIF" ]; then
    "$RUN" exec "$KRAKEN2_SIF" kraken2 --version >/dev/null 2>&1 \
      && ok "kraken2 runs inside container" || no "kraken2 not runnable in $KRAKEN2_SIF"
  fi
  have python3 && ok "python3 (for kraken-biom)" || wn "python3 not found in this env"
)

echo "=== STEP 04 env: Kaiju ==="
(
  export PATH="$KAIJU_BIN:$PATH"
  have kaiju              && ok "kaiju ($(command -v kaiju))" || no "kaiju not found (KAIJU_BIN=$KAIJU_BIN)"
  have kaiju2table        && ok "kaiju2table"                || no "kaiju2table not found"
  have kaiju-addTaxonNames && ok "kaiju-addTaxonNames"       || no "kaiju-addTaxonNames not found"
)

echo "=== STEP 05 env: R 4.4.2 + packages ==="
(
  module load R/4.4.2 2>/dev/null
  if have Rscript; then ok "Rscript -> $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
  else no "Rscript not found after 'module load R/4.4.2'"; fi
  if have Rscript; then
    miss=$(Rscript -e 'req<-c("readr","dplyr","tidyr","tibble","stringr","purrr","ggplot2","scales","phyloseq","biomformat"); cat(paste(req[!vapply(req,requireNamespace,logical(1),quietly=TRUE)],collapse=" "))' 2>/dev/null)
    [ -z "$miss" ] && ok "all required R packages present" || no "missing R packages: $miss  (run: Rscript install_r_deps.R)"
  fi
)

echo "=== Databases / references ==="
[ -e "$STAR_GENOME_DIR/SAindex" ]                          && ok "STAR index ($STAR_GENOME_DIR)"            || no "STAR index SAindex MISSING in $STAR_GENOME_DIR"
[ -e "$KRAKEN2_DB/taxo.k2d" ]                              && ok "Kraken2 DB taxo.k2d"                      || no "Kraken2 DB incomplete (no taxo.k2d): $KRAKEN2_DB"
[ -e "$KRAKEN2_DB/hash.k2d" ]                              && ok "Kraken2 DB hash.k2d"                      || no "Kraken2 DB hash.k2d MISSING"
[ -e "$KRAKEN2_DB/database${READ_LEN}mers.kmer_distrib" ] && ok "Bracken kmer_distrib (READ_LEN=$READ_LEN)" || no "Bracken needs database${READ_LEN}mers.kmer_distrib in $KRAKEN2_DB"
[ -e "$KAIJU_INDEX" ]                                      && ok "Kaiju .fmi index"                         || no "Kaiju index MISSING: $KAIJU_INDEX"
[ -e "$KAIJU_NODES" ]                                      && ok "Kaiju nodes.dmp"                          || no "nodes.dmp MISSING: $KAIJU_NODES"
[ -e "$KAIJU_NAMES" ]                                      && ok "Kaiju names.dmp"                          || no "names.dmp MISSING: $KAIJU_NAMES"
[ -e "$CONTAMINANT_BLOCKLIST" ]                            && ok "contaminant blocklist"                   || wn "blocklist missing: $CONTAMINANT_BLOCKLIST"
if [ -n "${HUMAN_RRNA_FASTA:-}" ] && [ -s "$HUMAN_RRNA_FASTA" ]; then
  ok "host rRNA reference present -> step-01 rRNA scrub ACTIVE"
else
  wn "no HUMAN_RRNA_FASTA -> rRNA scrub will SKIP (optional; build with build_human_rrna.sh)"
fi
if [ -n "${T2T_BT2_INDEX:-}" ] && ls "${T2T_BT2_INDEX}".*.bt2* >/dev/null 2>&1; then
  ok "T2T-CHM13 bowtie2 index present -> step-01 B4 ACTIVE"
else
  wn "no T2T index -> step-01 B4 will SKIP (optional; build with build_t2t_index.sh)"
fi

echo "=== Inputs ==="
if [ -n "${SAMPLE_MANIFEST:-}" ] && [ -s "${SAMPLE_MANIFEST:-/nonexistent}" ]; then
  nm=$(grep -cvE '^[[:space:]]*(#|$)' "$SAMPLE_MANIFEST")
  ok "sample manifest: $nm samples ($SAMPLE_MANIFEST) -- sheet is built from this"
fi
if [ -d "$TRIM_DIR" ]; then
  ntop=$(find -L "$TRIM_DIR" -maxdepth 1 -type f -name "${R1_TAG}*${R1_SUFFIX}" 2>/dev/null | wc -l)
  nall=$(find -L "$TRIM_DIR" -type f -name "${R1_TAG}*${R1_SUFFIX}" -not -path '*trash*' -not -path '*backup*' 2>/dev/null | wc -l)
  ok "trimmed R1 FASTQs: $ntop top-level in $TRIM_DIR ($nall incl. clean subfolders)"
else
  no "TRIM_DIR missing: $TRIM_DIR"
fi

P=$(grep -c '^OK' "$RES"); W=$(grep -c '^WN' "$RES"); F=$(grep -c '^NO' "$RES")
echo ""
echo "================  SUMMARY:  $P OK  /  $W WARN  /  $F FAIL  ================"
if [ "$F" -eq 0 ]; then echo ">>> READY. Submit with:  bash run_all.sh"; exit 0
else echo ">>> Resolve the FAIL item(s) above before submitting."; exit 1; fi
