#!/bin/bash
# =============================================================================
# build_human_rrna.sh -- assemble a HUMAN-ONLY rRNA reference for step 01.
# =============================================================================
# Step 01 uses this with bbduk to drop residual HOST ribosomal RNA reads that
# survive STAR (the human rDNA repeat is poorly assembled in GRCh38, so its
# reads often fail to map and would otherwise pollute the microbial pool).
#
# *** HUMAN-ONLY, on purpose. *** bbduk REMOVES reads matching this file, so it
# must contain only human rRNA. Do NOT substitute SILVA / SortMeRNA / a
# pan-kingdom rRNA set here -- that would also delete microbial 16S/23S, i.e.
# the signal you are trying to detect. (Human 18S vs bacterial 16S share almost
# no exact 31-mers, so a human-only reference leaves bacterial reads untouched.)
#
#   bash build_human_rrna.sh           # writes to $HUMAN_RRNA_FASTA from config.sh
#
# Run on a login node (needs outbound HTTPS to NCBI). Then re-run preflight.sh.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
OUT="${1:-$HUMAN_RRNA_FASTA}"
mkdir -p "$(dirname "$OUT")"

# Canonical human cytoplasmic rRNA sequences (NCBI nuccore):
#   U13369.1     complete human rDNA repeat unit -> 18S, 5.8S, 28S (+ spacers)
#   NR_046235.3  RNA45SN1  (45S pre-rRNA)
#   NR_023363.1  RNA5S1    (5S rRNA)
ACCS="U13369.1 NR_046235.3 NR_023363.1"
base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

tmp="$(mktemp)"; : > "$tmp"
for acc in $ACCS; do
    echo ">>> fetching $acc"
    wget -qO - "${base}?db=nuccore&id=${acc}&rettype=fasta&retmode=text" >> "$tmp" \
        || echo "    WARNING: failed to fetch $acc (skipped)" >&2
done

grep -v '^[[:space:]]*$' "$tmp" > "$OUT"      # strip blank separator lines
rm -f "$tmp"

n=$(grep -c '^>' "$OUT" || true)
bp=$(grep -v '^>' "$OUT" | tr -d '\n' | wc -c)
echo ">>> Wrote $OUT"
echo ">>> $n sequences, ${bp} bp (expect ~3 seqs, ~50-60 kb)"
[ "$n" -ge 1 ] || { echo "ERROR: nothing fetched -- check Sherlock outbound HTTPS/proxy." >&2; exit 1; }
echo ">>> config.sh HUMAN_RRNA_FASTA already points here; run 'bash preflight.sh' to confirm the rRNA scrub is ACTIVE."
