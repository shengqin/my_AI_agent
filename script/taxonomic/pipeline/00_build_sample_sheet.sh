#!/bin/bash
# =============================================================================
# 00_build_sample_sheet.sh -- discover trimmed FASTQs -> sample_list.tsv
# =============================================================================
# Builds a tab-separated table:  batch <TAB> sample <TAB> r1_path <TAB> r2_path
# (same format the lab's run_sf_array.sh / TRUST4.sh consume).
#
# Run this on the login node BEFORE submitting the array jobs:
#     bash 00_build_sample_sheet.sh
# Then set the array size of the step jobs to the number of lines printed.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

mkdir -p "$LOG_DIR"
: > "$SAMPLE_LIST"

[ -d "$TRIM_DIR" ] || { echo "ERROR: TRIM_DIR does not exist: $TRIM_DIR" >&2; exit 1; }

# --- Preferred path: build from the authoritative manifest (e.g. LP_samples.txt) ---
# Matches the expression pipeline and ignores stray/duplicate FASTQs in subfolders.
if [ -n "${SAMPLE_MANIFEST:-}" ] && [ -s "${SAMPLE_MANIFEST:-/nonexistent}" ]; then
    echo ">>> Building sample sheet from manifest: $SAMPLE_MANIFEST"
    r2_tag="${R1_TAG/R1/R2}"
    n=0; missing=0
    while IFS= read -r line || [ -n "$line" ]; do
        samp="$(printf '%s' "$line" | awk '{print $1}')"
        [ -z "$samp" ] && continue
        case "$samp" in \#*) continue ;; esac        # skip comment lines
        r1="$TRIM_DIR/${R1_TAG}${samp}${R1_SUFFIX}"
        r2="$TRIM_DIR/${r2_tag}${samp}${R1_SUFFIX}"
        if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
            echo "    WARNING: FASTQs missing for '$samp' -- skipping ($r1)" >&2
            missing=$((missing+1)); continue
        fi
        printf "all\t%s\t%s\t%s\n" "$samp" "$r1" "$r2" >> "$SAMPLE_LIST"
        n=$((n+1))
    done < "$SAMPLE_MANIFEST"
    [ "$n" -gt 0 ] || { echo "ERROR: manifest produced 0 usable samples (check paths/naming)." >&2; exit 1; }
    dups="$(awk -F'\t' '{print $2}' "$SAMPLE_LIST" | sort | uniq -d || true)"
    [ -n "$dups" ] && { echo "    WARNING: duplicate sample names:" >&2; echo "$dups" | sed 's/^/        /' >&2; }
    echo ">>> Wrote $n samples from manifest ($missing missing) to: $SAMPLE_LIST"
    echo ">>> Submit the array steps with:  --array=1-${n}%8"
    exit 0
fi

echo ">>> No manifest set; discovering trimmed FASTQs under: $TRIM_DIR"

# Find R1 files. Primary pattern: R1_<sample>.trimmed.fastq.gz
# Recurse so it works whether files are flat in NovaSeq_trimmed/ or under
# per-batch subdirs (NovaSeqB*/trimmed/, etc.).
# NB: prune junk dirs (e.g. temp_trash/) so stray FASTQs are never picked up.
mapfile -t R1S < <(find -L "$TRIM_DIR" -type f -name "${R1_TAG}*${R1_SUFFIX}" \
                     -not -path '*trash*' -not -path '*backup*' | sort)

# Fallback to the common Illumina pattern <sample>_R1*.fastq.gz if none found.
if [ "${#R1S[@]}" -eq 0 ]; then
    echo "    No ${R1_TAG}*${R1_SUFFIX} files; trying *_R1*.fastq.gz ..."
    mapfile -t R1S < <(find -L "$TRIM_DIR" -type f -name "*_R1*.fastq.gz" \
                         -not -path '*trash*' -not -path '*backup*' | sort)
fi

[ "${#R1S[@]}" -gt 0 ] || { echo "ERROR: no R1 FASTQs found in $TRIM_DIR" >&2; exit 1; }

n=0
for r1 in "${R1S[@]}"; do
    # Derive the mate R2 by swapping the first R1 token -> R2.
    if [[ "$(basename "$r1")" == ${R1_TAG}* ]]; then
        r2="$(dirname "$r1")/$(basename "$r1" | sed "s/^${R1_TAG}/R2_/")"
        sample="$(basename "$r1" "$R1_SUFFIX")"; sample="${sample#${R1_TAG}}"
    else
        r2="${r1/_R1/_R2}"
        sample="$(basename "$r1" | sed -E 's/_R1.*//')"
    fi

    if [ ! -f "$r2" ]; then
        echo "    WARNING: R2 missing for $sample (expected $r2) -- skipping" >&2
        continue
    fi

    # Batch = name of the parent dir two levels up if it looks like NovaSeq*, else "all".
    parent="$(basename "$(dirname "$r1")")"
    grandparent="$(basename "$(dirname "$(dirname "$r1")")")"
    if [[ "$grandparent" == NovaSeq* ]]; then batch="$grandparent"
    elif [[ "$parent" == NovaSeq* ]]; then batch="$parent"
    else batch="all"; fi

    printf "%s\t%s\t%s\t%s\n" "$batch" "$sample" "$r1" "$r2" >> "$SAMPLE_LIST"
    n=$((n+1))
done

# Warn on duplicate sample names (would collide in output dirs).
dups="$(awk -F'\t' '{print $2}' "$SAMPLE_LIST" | sort | uniq -d || true)"
if [ -n "$dups" ]; then
    echo "    WARNING: duplicate sample names detected:" >&2
    echo "$dups" | sed 's/^/        /' >&2
fi

echo ">>> Wrote $n samples to: $SAMPLE_LIST"
echo ">>> Submit the array steps with:  --array=1-${n}%8"
