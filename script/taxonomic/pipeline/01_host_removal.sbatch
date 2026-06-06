#!/bin/bash
#SBATCH --job-name=mb_host_removal
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=emoding
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ssu42
#SBATCH --output=logs/01_host_%A_%a.out
#SBATCH --error=logs/01_host_%A_%a.err
# -----------------------------------------------------------------------------
# STEP 01 -- Host (human) read removal  ==> candidate-microbial FASTQs
# -----------------------------------------------------------------------------
# Submit AFTER building the sample sheet:
#     sbatch --array=1-<N>%8 01_host_removal.sbatch
#
# What it does (and WHY), for each sample:
#   A. STAR -> GRCh38 (splice-aware). Reads that map are human; the UNMAPPED
#      reads are the candidate microbial pool. STAR's Log.final.out gives the
#      total input-read count we use later for RPM normalization.
#   B. (recommended, auto-skipped if tools/refs absent) clean the unmapped pool:
#         - bbduk rRNA k-mer removal (host rRNA survives step A and otherwise
#           dominates Pico-kit libraries),
#         - bbduk low-complexity (entropy) filter -- the #1 source of bogus
#           taxonomic hits on FFPE data,
#         - clumpify exact-duplicate (PCR) removal (low-input FFPE -> many dups;
#           optical dedup is not enabled -- see DEDUP note in config.sh),
#         - optional bowtie2 vs T2T-CHM13 to trap residual human GRCh38 misses.
#   C. Emit <sample>_microbe_R1/R2.fastq and a per-sample read-accounting line.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"
source "$SCRIPT_DIR/config.sh"

# ---- modules / tools (faithful to lab convention; optional tools guarded) ----
# NOTE: do NOT `module --force purge` here -- it breaks the Sherlock module env.
module load biology 2>/dev/null || true
module load samtools 2>/dev/null || true
module load biology bbmap 2>/dev/null || true
module load biology bowtie2 2>/dev/null || true

# Use the locally-installed STAR 2.7.11b (matches bs_STAR_align.sh and the index
# build), NOT the older `star` module which would shadow it.
export PATH="$STAR_BIN_DIR:$PATH"
# Optional BBMap install (bbduk.sh, clumpify.sh) if the module didn't provide them.
[ -n "${BBMAP_BIN:-}" ] && export PATH="$BBMAP_BIN:$PATH"
have() { command -v "$1" >/dev/null 2>&1; }
have STAR || { echo "ERROR: STAR not found (check STAR_BIN_DIR=$STAR_BIN_DIR)" >&2; exit 1; }
echo ">>> STAR: $(STAR --version 2>&1) from $(command -v STAR)"

mkdir -p "$D_HOST" "$D_HOST/star" "$D_HOST/readcounts" "$LOG_DIR" logs

# ---- pull this array task's sample from the sheet ----
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
batch=$(echo "$LINE" | cut -f1)
sample=$(echo "$LINE" | cut -f2)
r1=$(echo "$LINE" | cut -f3)
r2=$(echo "$LINE" | cut -f4)
THREADS=${SLURM_CPUS_PER_TASK:-12}

echo "=================================================="
echo " Sample : $sample   (batch $batch, array ${SLURM_ARRAY_TASK_ID:-NA})"
echo "=================================================="

out_r1="$D_HOST/${sample}_microbe_R1.fastq"
out_r2="$D_HOST/${sample}_microbe_R2.fastq"
if [ -s "$out_r1" ] && [ -s "$out_r2" ]; then
    echo "[SKIP] $sample already host-removed."; exit 0
fi

work="$D_HOST/star/${sample}_"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/${sample}.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

# ---- STEP A: STAR de-hosting ----
# Mirrors bs_STAR_align.sh (2-pass, FFPE-relaxed) but adds --outReadsUnmapped
# Fastx (the unmapped mates are the microbial pool) and omits TranscriptomeSAM.
# The prebuilt index already contains the annotation, so a GTF is only added if
# STAR_GTF is set. For host removal, relaxed filters are ON by default so that
# more short FFPE host reads are mapped (= fewer human reads leak downstream).
echo ">>> [A] STAR -> GRCh38 (GENCODE v47)"
# --outSAMtype None: we never need the host alignment BAM, only the UNMAPPED
# reads (the microbial pool) and Log.final.out (read-count denominator for RPM).
# STAR still aligns and reports stats; it just skips writing/sorting the BAM,
# which is faster and avoids gigabytes of host reads we would only delete.
STAR_OPTS=()
[ -n "${STAR_GTF:-}" ] && [ -s "${STAR_GTF:-/nonexistent}" ] && STAR_OPTS+=(--sjdbGTFfile "$STAR_GTF")
[ "${STAR_TWOPASS:-false}" = "true" ] && STAR_OPTS+=(--twopassMode Basic)
if [ "${STAR_FFPE_RELAXED:-false}" = "true" ]; then
    STAR_OPTS+=(--outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --outFilterMatchNmin 20)
fi
# Dense seeding (seedSearchStartLmax 25) was the main slowdown; default 50 is ~10x
# faster with negligible host-capture loss. Only added if STAR_SEED_LMAX is set.
[ -n "${STAR_SEED_LMAX:-}" ] && STAR_OPTS+=(--seedSearchStartLmax "$STAR_SEED_LMAX")
STAR --runThreadN "$THREADS" \
     --runMode alignReads \
     --genomeDir "$STAR_GENOME_DIR" \
     --readFilesIn "$r1" "$r2" \
     --readFilesCommand zcat \
     --outSAMtype None \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix "$work" \
     ${STAR_OPTS[@]+"${STAR_OPTS[@]}"}

mate1="${work}Unmapped.out.mate1"
mate2="${work}Unmapped.out.mate2"
[ -s "$mate1" ] || { echo "ERROR: STAR produced no unmapped mate1 for $sample" >&2; exit 1; }

# total input reads (read PAIRS) for RPM normalization later
total_reads=$(grep "Number of input reads" "${work}Log.final.out" | awk '{print $NF}')
star_unmapped=$(($(wc -l < "$mate1") / 4))
echo "    input read pairs=$total_reads ; unmapped pairs=$star_unmapped"

cur1="$mate1"; cur2="$mate2"

# ---- STEP B: optional cleanup of the unmapped pool ----
after_rrna="$star_unmapped"; after_entropy="$star_unmapped"; after_dedup="$star_unmapped"; after_t2t="$star_unmapped"

# B1. rRNA removal (bbduk) -- removes pairs whose reads k-mer-match host rRNA
if have bbduk.sh && [ -s "$HUMAN_RRNA_FASTA" ]; then
    echo ">>> [B1] bbduk rRNA removal"
    bbduk.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/no_rrna_1.fq" out2="$tmp/no_rrna_2.fq" \
        ref="$HUMAN_RRNA_FASTA" k=31 1>>"$work"clean.log 2>&1
    cur1="$tmp/no_rrna_1.fq"; cur2="$tmp/no_rrna_2.fq"
    after_rrna=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B1] SKIP rRNA removal (bbduk or HUMAN_RRNA_FASTA missing)"
fi

# B2. low-complexity filter (bbduk entropy)
if have bbduk.sh && [ "$LOWCOMPLEXITY_ENTROPY" != "0" ]; then
    echo ">>> [B2] bbduk entropy filter (>$LOWCOMPLEXITY_ENTROPY)"
    bbduk.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/hc_1.fq" out2="$tmp/hc_2.fq" \
        entropy="$LOWCOMPLEXITY_ENTROPY" entropywindow=50 entropyk=5 1>>"$work"clean.log 2>&1
    cur1="$tmp/hc_1.fq"; cur2="$tmp/hc_2.fq"
    after_entropy=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B2] SKIP entropy filter"; after_entropy="$after_rrna"
fi

# B3. duplicate removal (clumpify)
if have clumpify.sh && [ "$DEDUP" = "true" ]; then
    echo ">>> [B3] clumpify dedup"
    clumpify.sh -Xmx32g threads="$THREADS" \
        in1="$cur1" in2="$cur2" \
        out1="$tmp/dd_1.fq" out2="$tmp/dd_2.fq" \
        dedupe=t 1>>"$work"clean.log 2>&1
    cur1="$tmp/dd_1.fq"; cur2="$tmp/dd_2.fq"
    after_dedup=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B3] SKIP dedup"; after_dedup="$after_entropy"
fi

# B4. optional second-host removal vs T2T-CHM13 (keep pairs that DON'T map)
if [ -n "$T2T_BT2_INDEX" ] && have bowtie2; then
    echo ">>> [B4] bowtie2 vs T2T-CHM13 (residual human)"
    bowtie2 -x "$T2T_BT2_INDEX" -1 "$cur1" -2 "$cur2" \
        --very-sensitive-local -p "$THREADS" \
        --un-conc "$tmp/noT2T_%.fq" -S /dev/null 2>>"$work"clean.log
    cur1="$tmp/noT2T_1.fq"; cur2="$tmp/noT2T_2.fq"
    after_t2t=$(($(wc -l < "$cur1") / 4))
else
    echo ">>> [B4] SKIP T2T removal (T2T_BT2_INDEX unset or bowtie2 missing)"; after_t2t="$after_dedup"
fi

# ---- STEP C: finalize ----
mv -f "$cur1" "$out_r1"
mv -f "$cur2" "$out_r2"
final=$(($(wc -l < "$out_r1") / 4))

# Per-sample read accounting WITH a header row. final_microbial == after_t2t.
{ printf "sample\tinput_pairs\tstar_unmapped\tafter_rrna\tafter_entropy\tafter_dedup\tafter_t2t\tfinal_microbial\n"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$sample" "$total_reads" "$star_unmapped" "$after_rrna" "$after_entropy" "$after_dedup" "$after_t2t" "$final"
} > "$D_HOST/readcounts/${sample}.readcount.tsv"

# Discard STAR scratch. The single Log.final.out stays in star/ for debugging;
# total reads are taken from readcount.tsv (input_pairs), so no duplicate copy.
rm -rf "${work}Unmapped.out.mate"* "${work}_STARpass1" "${work}_STARtmp" \
       "${work}_STARgenome" "${work}SJ.out.tab" "${work}Log.out" \
       "${work}Log.progress.out" "${work}Aligned"* 2>/dev/null || true

echo ">>> DONE $sample : final candidate-microbial pairs = $final"
