#!/bin/bash
#SBATCH --job-name=mb_build_t2t
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=emoding
#SBATCH --output=logs/build_t2t_%j.out
#SBATCH --error=logs/build_t2t_%j.err
# =============================================================================
# build_t2t_index.sh -- obtain a bowtie2 index of T2T-CHM13v2.0 for step-01 B4.
# =============================================================================
# B4 removes residual human reads that GRCh38 (and therefore STAR) misses, by
# aligning the survivors to the complete T2T-CHM13v2.0 assembly with bowtie2.
# This complements STAR (which is splice-aware and does the bulk RNA host removal);
# the residual reads reaching B4 are mostly unspliced repeat/genomic human, which
# bowtie2 handles well.
#
# Prefers the PREBUILT bowtie2 index from genome-idx (fast: download + unzip).
# Falls back to the FASTA + bowtie2-build (slow -> submit with sbatch).
#
#   bash build_t2t_index.sh      # prebuilt download is light; login node is fine
#   sbatch build_t2t_index.sh    # use this if it must build from FASTA
#
# Writes to $(dirname $T2T_BT2_INDEX); index prefix = $T2T_BT2_INDEX. Then preflight.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
module load biology bowtie2 2>/dev/null || true

refdir="$(dirname "$T2T_BT2_INDEX")"
prefix="$(basename "$T2T_BT2_INDEX")"
mkdir -p "$refdir"; cd "$refdir"

if ls "${prefix}".*.bt2* >/dev/null 2>&1; then
    echo ">>> T2T bowtie2 index already present: $T2T_BT2_INDEX"; exit 0
fi

echo ">>> Trying prebuilt bowtie2 index (genome-idx)..."
if wget -c "https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip" -O chm13v2.0.zip; then
    unzip -o chm13v2.0.zip && rm -f chm13v2.0.zip || true
    # the zip unpacks into a chm13v2.0/ subdir -> lift the .bt2 files up to refdir
    if [ -d "$prefix" ] && ls "$prefix"/*.bt2* >/dev/null 2>&1; then
        mv -f "$prefix"/*.bt2* . && rmdir "$prefix" 2>/dev/null || true
    fi
fi

# If we still don't have a usable index, build from FASTA.
if ! ls "${prefix}".*.bt2* >/dev/null 2>&1; then
    echo ">>> Prebuilt index unavailable; building from FASTA (slow -> use sbatch)..."
    command -v bowtie2-build >/dev/null 2>&1 || { echo "ERROR: bowtie2-build not on PATH (module load biology bowtie2)" >&2; exit 1; }
    wget -c "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz" -O chm13v2.0.fa.gz
    gunzip -c chm13v2.0.fa.gz > chm13v2.0.fa   # -c: portable (old gzip lacks -k)
    bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" chm13v2.0.fa "$prefix"
fi

if ls "${prefix}".*.bt2* >/dev/null 2>&1; then
    echo ">>> Done. config.sh T2T_BT2_INDEX already points here: $T2T_BT2_INDEX"
    echo ">>> Run 'bash preflight.sh' to confirm step-01 B4 is ACTIVE."
else
    echo "ERROR: failed to obtain T2T bowtie2 index." >&2; exit 1
fi
