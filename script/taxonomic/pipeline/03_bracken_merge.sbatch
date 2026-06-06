#!/bin/bash
#SBATCH --job-name=mb_bracken_merge
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/03_bracken_%j.out
#SBATCH --error=logs/03_bracken_%j.err
# -----------------------------------------------------------------------------
# STEP 03 -- Bracken abundance re-estimation + cohort merge  (single job)
# -----------------------------------------------------------------------------
#     sbatch 03_bracken_merge.sbatch        # run after ALL step-02 array tasks
#
# Bracken redistributes Kraken2's reads to give proper species/genus abundances.
# We then (a) merge per-sample reports into BIOM matrices (for phyloseq/decontam)
# and (b) build a distinct-minimizer table from the Kraken2 reports so step 05
# can drop taxa supported by too few unique k-mers (KrakenUniq-style FP control).
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
mkdir -p "$D_BRACKEN" logs

# ---- A. Bracken per sample, per level ----
shopt -s nullglob
reports=("$D_KRAKEN"/*.kreport)
[ "${#reports[@]}" -gt 0 ] || { echo "ERROR: no .kreport files in $D_KRAKEN" >&2; exit 1; }

bfail=0
for level in $BRACKEN_LEVELS; do
    suffix=$([ "$level" = "S" ] && echo "species" || echo "genus")
    echo ">>> Bracken level=$level ($suffix)"
    for rep in "${reports[@]}"; do
        base=$(basename "$rep" .kreport)
        out_b="$D_BRACKEN/${base}.${suffix}.bracken"
        out_k="$D_BRACKEN/${base}.${suffix}_bracken.kreport"
        [ -s "$out_k" ] && { echo "    [skip] $base ($suffix)"; continue; }
        "$CONTAINER_RUN" exec "$BRACKEN_SIF" bracken \
            -d "$KRAKEN2_DB" -i "$rep" \
            -o "$out_b" -w "$out_k" \
            -r "$READ_LEN" -l "$level" -t "$BRACKEN_THRESH" \
            || { echo "    WARNING: bracken failed for $base ($suffix)"; bfail=$((bfail+1)); }
    done
done
[ "$bfail" -gt 0 ] && echo "    NOTE: $bfail bracken sample-levels failed (see warnings above); cohort matrices are built from the rest."

# ---- B. Merge into BIOM matrices ----
echo ">>> Merging Bracken reports into BIOM matrices (kraken-biom)"
module load python/3.9.0 2>/dev/null || true
pip install --user --quiet kraken-biom 2>/dev/null || true
KB="$HOME/.local/bin/kraken-biom"
command -v kraken-biom >/dev/null 2>&1 && KB="kraken-biom"

for level in $BRACKEN_LEVELS; do
    suffix=$([ "$level" = "S" ] && echo "species" || echo "genus")
    kfiles=("$D_BRACKEN"/*.${suffix}_bracken.kreport)
    if [ "${#kfiles[@]}" -gt 0 ] && [ -x "$KB" -o "$KB" = "kraken-biom" ]; then
        "$KB" "${kfiles[@]}" --fmt json -o "$D_BRACKEN/NLPHL_microbiome_${suffix}.biom" \
            && echo "    wrote NLPHL_microbiome_${suffix}.biom" \
            || echo "    WARNING: kraken-biom failed for $suffix"
    else
        echo "    WARNING: no ${suffix}_bracken.kreport files, or kraken-biom unavailable -- skipped $suffix BIOM"
    fi
done

# Step 05 hard-requires the species BIOM. Fail loudly HERE rather than letting
# the dependency chain "succeed" and then crashing inside R on a missing file.
if [[ " $BRACKEN_LEVELS " == *" S "* ]] && [ ! -s "$D_BRACKEN/NLPHL_microbiome_species.biom" ]; then
    echo "ERROR: species BIOM was not produced ($D_BRACKEN/NLPHL_microbiome_species.biom)." >&2
    echo "       Check the kraken-biom install (pip install --user kraken-biom) and the bracken output above." >&2
    exit 1
fi

# ---- C. Distinct-minimizer table (KrakenUniq-style FP control) ----
# From each 8-column minimizer report, keep species rows:
#   sample, taxid, name, taxon_reads, total_minimizers, distinct_minimizers
echo ">>> Building distinct-minimizer table"
mmz_out="$D_BRACKEN/minimizer_stats.tsv"
printf "sample\ttaxid\tname\ttaxon_reads\ttotal_minimizers\tdistinct_minimizers\n" > "$mmz_out"
for mr in "$D_KRAKEN"/*.kreport.minimizer; do
    [ -e "$mr" ] || { echo "    (no minimizer reports; skipping)"; break; }
    s=$(basename "$mr" .kreport.minimizer)
    # columns: 1%  2clade_reads  3taxon_reads  4tot_minz  5distinct_minz  6rank  7taxid  8name
    awk -F'\t' -v S="$s" 'BEGIN{OFS="\t"} $6=="S" {
        name=$8; gsub(/^ +/,"",name);
        print S,$7,name,$3,$4,$5 }' "$mr" >> "$mmz_out"
done

# ---- D. STAR total-read lookup (for RPM), mirroring the scratch logic ----
echo ">>> Building total_reads_lookup.tsv (from readcount.tsv input_pairs column)"
lookup="$D_BRACKEN/total_reads_lookup.tsv"
: > "$lookup"
for rc in "$D_HOST"/readcounts/*.readcount.tsv; do
    [ -e "$rc" ] || break
    awk -F'\t' 'NR>1 {print $1"\t"$2}' "$rc" >> "$lookup"   # sample, input_pairs (skip header)
done

# ---- E. Kraken2 classification QC: classified vs unclassified per sample ----
# After --confidence/--minimum-hit-groups, reads that fail (or have no DB hit)
# are reported as "unclassified" (taxid 0). Root (taxid 1) clade-reads = classified.
echo ">>> Building kraken2_classification_QC.tsv"
kqc="$D_BRACKEN/kraken2_classification_QC.tsv"
printf "sample\tclassified_pairs\tunclassified_pairs\ttotal_pairs\tpct_unclassified\n" > "$kqc"
for rep in "${reports[@]}"; do
    base=$(basename "$rep" .kreport)
    u=$(awk -F'\t' '$5=="0"{print $2; exit}' "$rep"); u=${u:-0}   # unclassified (taxid 0)
    r=$(awk -F'\t' '$5=="1"{print $2; exit}' "$rep"); r=${r:-0}   # classified  (root, taxid 1)
    tot=$((u + r))
    pct=$(awk -v u="$u" -v t="$tot" 'BEGIN{ if(t>0) printf "%.2f",100*u/t; else printf "NA" }')
    printf "%s\t%s\t%s\t%s\t%s\n" "$base" "$r" "$u" "$tot" "$pct" >> "$kqc"
done

echo ">>> STEP 03 complete."
