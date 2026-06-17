#!/bin/bash
# ==============================================================================
# Optional two-pass CNVkit purity recall
#
# Workflow:
#   1. Run the full pipeline normally (CNVkit call without --purity).
#   2. Run summarize_variant_calls.sh.
#   3. Run test.R, which writes analysis/tumor_fraction.csv from confident SNV VAFs.
#   4. Run this script to re-call CNVkit segments with per-sample --purity.
#   5. Re-run summarize_variant_calls.sh, then re-run test.R.
#
# CNVkit --purity mainly improves integer copy-number interpretation. The adaptive
# low-tumor-fraction thresholds in test.R are the detection-oriented step.
#
# Usage:
#   recall_cnvkit_with_purity.sh <tumor_fraction.csv> <results_dir1> [results_dir2 ...]
# ==============================================================================

set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") <tumor_fraction.csv> <results_dir1> [results_dir2 ...]" >&2
}

if [[ "$#" -lt 2 ]]; then
    usage
    exit 1
fi

TF_CSV="$1"
shift
RESULTS_DIRS=("$@")

if [[ ! -f "$TF_CSV" ]]; then
    echo "ERROR: tumor-fraction CSV not found: $TF_CSV" >&2
    exit 1
fi

for dir in "${RESULTS_DIRS[@]}"; do
    if [[ ! -d "$dir" ]]; then
        echo "ERROR: results directory not found: $dir" >&2
        exit 1
    fi
done

header=$(head -n 1 "$TF_CSV" | tr -d '\r')
if [[ "$header" != "Sample,N_SNV_for_TF,Tumor_Fraction,TF_Method" ]]; then
    echo "ERROR: unexpected tumor-fraction CSV header:" >&2
    echo "  $header" >&2
    echo "Expected: Sample,N_SNV_for_TF,Tumor_Fraction,TF_Method" >&2
    exit 1
fi

is_numeric_ge_005() {
    local value="$1"
    awk -v x="$value" 'BEGIN { exit !(x ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && x + 0 >= 0.05) }'
}

find_single_sorted() {
    local dir="$1"
    shift
    local matches n
    matches=$(find "$dir" "$@" -print 2>/dev/null | sort || true)
    n=$(printf '%s\n' "$matches" | grep -c . || true)
    if [[ "$n" -eq 1 ]]; then
        printf '%s\n' "$matches"
        return 0
    fi
    return 1
}

find_cnvkit_run_dir() {
    local sample="$1"
    local matches=()
    local results_dir candidate

    for results_dir in "${RESULTS_DIRS[@]}"; do
        candidate="$results_dir/$sample/cnvkit_run"
        if [[ -d "$candidate" ]]; then
            matches+=("$candidate")
        fi
    done

    if [[ "${#matches[@]}" -eq 1 ]]; then
        printf '%s\n' "${matches[0]}"
        return 0
    fi

    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "  SKIP $sample: found multiple cnvkit_run directories; refusing ambiguous recall." >&2
        printf '    %s\n' "${matches[@]}" >&2
    else
        echo "  SKIP $sample: cnvkit_run directory not found under supplied results dirs." >&2
    fi
    return 1
}

run_cnvkit_with_purity() {
    local base_cns="$1"
    local cnr="$2"
    local purity="$3"
    local out_call="$4"
    local out_bintest="$5"

    # Keep the conda activation inside this subprocess so the caller's shell state
    # is unchanged and the CNVkit environment is explicit/reproducible.
    bash -lc '
        set -euo pipefail
        source ~/miniconda3/etc/profile.d/conda.sh
        conda activate cnvkit
        cnvkit.py call "$1" --purity "$2" -o "$3"
        cnvkit.py bintest "$4" -s "$3" > "$5"
    ' _ "$base_cns" "$purity" "$out_call" "$cnr" "$out_bintest"
}

processed=0
skipped=0

while IFS=, read -r sample n_snv tumor_fraction tf_method extra; do
    sample=${sample//$'\r'/}
    tumor_fraction=${tumor_fraction//$'\r'/}

    if [[ -z "$sample" ]]; then
        continue
    fi

    if [[ -n "${extra:-}" ]]; then
        echo "  SKIP $sample: malformed CSV row has extra comma-separated fields." >&2
        skipped=$((skipped + 1))
        continue
    fi

    if ! is_numeric_ge_005 "$tumor_fraction"; then
        echo "  SKIP $sample: Tumor_Fraction is NA, non-numeric, or <0.05 ($tumor_fraction)."
        skipped=$((skipped + 1))
        continue
    fi

    if ! cnvkit_dir=$(find_cnvkit_run_dir "$sample"); then
        skipped=$((skipped + 1))
        continue
    fi

    if ! base_cns=$(find_single_sorted "$cnvkit_dir" -type f -name "*.cns" ! -name "*.call.cns" ! -name "*.bintest.cns"); then
        echo "  SKIP $sample: expected exactly one base .cns in $cnvkit_dir (excluding .call/.bintest)." >&2
        skipped=$((skipped + 1))
        continue
    fi

    if ! cnr=$(find_single_sorted "$cnvkit_dir" -type f -name "*.cnr"); then
        echo "  SKIP $sample: expected exactly one .cnr in $cnvkit_dir." >&2
        skipped=$((skipped + 1))
        continue
    fi

    stem=${base_cns%.cns}
    call_cns="${stem}.call.cns"
    bintest_cns="${stem}.bintest.cns"
    backup="${call_cns}.prepurity.bak"

    if [[ -f "$call_cns" && ! -f "$backup" ]]; then
        cp -p "$call_cns" "$backup"
        echo "  $sample: backed up $(basename "$call_cns") -> $(basename "$backup")"
    elif [[ -f "$backup" ]]; then
        echo "  $sample: backup already exists ($(basename "$backup")); leaving it unchanged."
    fi

    echo "  $sample: CNVkit call --purity $tumor_fraction"
    run_cnvkit_with_purity "$base_cns" "$cnr" "$tumor_fraction" "$call_cns" "$bintest_cns"
    processed=$((processed + 1))
done < <(tail -n +2 "$TF_CSV")

echo "CNVkit purity recall complete: $processed processed, $skipped skipped."
