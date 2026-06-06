#!/bin/bash
#SBATCH --job-name=mb_kaiju
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/04_kaiju_%A_%a.out
#SBATCH --error=logs/04_kaiju_%A_%a.err
# -----------------------------------------------------------------------------
# STEP 04 -- Kaiju protein-level classification (orthogonal to Kraken2)
# -----------------------------------------------------------------------------
#     sbatch --array=1-<N>%2 04_kaiju.sbatch
# nr_euk is held in RAM (~187 GB) -> mem >200G and LOW concurrency (%2).
#
# Kaiju translates reads in 6 frames and searches a protein DB (nr_euk), so it
# is robust to FFPE nucleotide damage and good for divergent taxa/viruses --
# but it cannot classify rRNA reads (non-coding). Comparing Kaiju with Kraken2
# is exactly how we separate believable signal from artifact (step 05).
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"
export PATH="$KAIJU_BIN:$PATH"
mkdir -p "$D_KAIJU" logs

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
sample=$(echo "$LINE" | cut -f2)
THREADS=${SLURM_CPUS_PER_TASK:-16}

r1="$D_HOST/${sample}_microbe_R1.fastq"
r2="$D_HOST/${sample}_microbe_R2.fastq"
[ -s "$r1" ] || { echo "ERROR: missing host-removed reads for $sample" >&2; exit 1; }

out="$D_KAIJU/${sample}.kaiju.out"
[ -s "$out" ] && { echo "[SKIP] $sample kaiju done."; exit 0; }

echo ">>> Kaiju ($KAIJU_MODE, -e $KAIJU_MISMATCHES) on $sample"
kaiju -t "$KAIJU_NODES" -f "$KAIJU_INDEX" \
      -i "$r1" -j "$r2" \
      -o "$out" \
      -z "$THREADS" -a "$KAIJU_MODE" -e "$KAIJU_MISMATCHES"

kaiju-addTaxonNames -t "$KAIJU_NODES" -n "$KAIJU_NAMES" \
      -i "$out" -o "$D_KAIJU/${sample}.names.out" \
      -r superkingdom,phylum,class,order,family,genus,species

echo ">>> DONE kaiju $sample"
