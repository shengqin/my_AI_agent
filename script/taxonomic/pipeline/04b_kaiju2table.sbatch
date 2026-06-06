#!/bin/bash
#SBATCH --job-name=mb_kaiju2table
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/04b_kaiju2table_%j.out
#SBATCH --error=logs/04b_kaiju2table_%j.err
# -----------------------------------------------------------------------------
# STEP 04b -- Merge per-sample Kaiju outputs into cohort matrices (single job)
# -----------------------------------------------------------------------------
#     sbatch 04b_kaiju2table.sbatch        # after ALL step-04 array tasks
# -e drops reads classified to taxa not present at the requested rank.
# -l writes the full lineage so step 05 can match Kaiju to Kraken2 by lineage.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
export PATH="$KAIJU_BIN:$PATH"

shopt -s nullglob
kfiles=("$D_KAIJU"/*.kaiju.out)
[ "${#kfiles[@]}" -gt 0 ] || { echo "ERROR: no .kaiju.out in $D_KAIJU" >&2; exit 1; }

echo ">>> kaiju2table species (full lineage)"
kaiju2table -e -t "$KAIJU_NODES" -n "$KAIJU_NAMES" -r species \
    -l superkingdom,phylum,class,order,family,genus,species \
    -o "$D_KAIJU/cohort_full_lineage_matrix.tsv" "${kfiles[@]}"

echo ">>> kaiju2table genus"
kaiju2table -e -t "$KAIJU_NODES" -n "$KAIJU_NAMES" -r genus \
    -l superkingdom,phylum,class,order,family,genus \
    -o "$D_KAIJU/cohort_genus_matrix.tsv" "${kfiles[@]}"

# ---- Kaiju classification QC: classified (C) vs unclassified (U) per sample ----
# Mirrors the Kraken2 QC so you can compare how much each classifier could call.
echo ">>> Building kaiju_classification_QC.tsv"
kjqc="$D_KAIJU/kaiju_classification_QC.tsv"
printf "sample\tkaiju_classified\tkaiju_unclassified\tkaiju_total\tpct_kaiju_classified\n" > "$kjqc"
for f in "${kfiles[@]}"; do
    s=$(basename "$f" .kaiju.out)
    awk -v s="$s" 'BEGIN{OFS="\t"}
        {if($1=="C")c++; else if($1=="U")u++}
        END{ t=c+u; pct=(t>0)?sprintf("%.2f",100*c/t):"NA"; print s, c+0, u+0, t, pct }' "$f" >> "$kjqc"
done

echo ">>> STEP 04b complete -> $D_KAIJU/cohort_full_lineage_matrix.tsv"
