#!/bin/bash
# =============================================================================
# config.sh  --  Central configuration for the NLPHL tumor-microbiome pipeline
# =============================================================================
# Every step script sources this file, so all paths/parameters live in ONE place.
# Edit the values below to match your environment, then submit the steps.
#
# Convention note: this mirrors the lab's existing RNA-seq jobs
# (run_sf_array.sh, TRUST4.sh): SLURM array jobs, partition=emoding,
# mail-user=ssu42, conda envs under brian/conda_envs, tools under brian/tools,
# and a tab-separated sample_list.tsv ( batch <TAB> sample <TAB> r1 <TAB> r2 ).
# =============================================================================

# Resolve the directory this pipeline lives in (so scripts can find siblings).
export PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# ---------------------------------------------------------------------------
# 1. COHORT INPUT / OUTPUT
# ---------------------------------------------------------------------------
export COHORT_DIR="/oak/stanford/groups/emoding/sequencing/pipeline/runs/rnaseq/Binkley_NLPHL"

# Demuxed + fastp-trimmed paired FASTQs (the same inputs used for STAR-RSEM / salmon).
export TRIM_DIR="$COHORT_DIR/NovaSeq_trimmed"

# All pipeline outputs are written under here.
export RUN_DIR="$COHORT_DIR/microbiome"
export LOG_DIR="$RUN_DIR/logs"
export SAMPLE_LIST="$LOG_DIR/sample_list.tsv"

# Authoritative sample list. If set AND present, the sample sheet is built from
# this manifest (one sample name per line; column 1) instead of globbing the
# directory -- this matches the expression pipeline (bs_STAR_align.sh reads the
# same file) and avoids stray/duplicate FASTQs in subfolders. Empty => glob.
export SAMPLE_MANIFEST="$TRIM_DIR/LP_samples.txt"

# Per-step output sub-directories.
export D_HOST="$RUN_DIR/01_host_removed"      # candidate-microbial FASTQs + STAR logs
export D_KRAKEN="$RUN_DIR/02_kraken2"         # .kraken / .kreport
export D_BRACKEN="$RUN_DIR/03_bracken"        # .bracken + merged matrices/biom
export D_KAIJU="$RUN_DIR/04_kaiju"            # .kaiju.out + kaiju2table matrices
export D_RESULTS="$RUN_DIR/05_results"        # filtered, cross-validated tables
export D_FIGS="$RUN_DIR/06_figures"           # plots / PDFs

# FASTQ naming. Trimmed files are R1_<sample>.trimmed.fastq.gz (lab convention).
# The sample-sheet builder also tries a few common alternates automatically.
export R1_TAG="R1_"
export R1_SUFFIX=".trimmed.fastq.gz"

# ---------------------------------------------------------------------------
# 2. COMPUTE (SLURM)
# ---------------------------------------------------------------------------
export PARTITION="emoding"
export MAIL_USER="ssu42"

# ---------------------------------------------------------------------------
# 3. HOST / HUMAN REMOVAL REFERENCES   (step 01)
# ---------------------------------------------------------------------------
# STAR index used for host depletion (GRCh38 / GENCODE v47 — the SAME index your
# bs_STAR_align.sh expression pipeline uses). Built WITH its annotation (sjdb
# baked in), so a GTF is NOT required at mapping time; leave STAR_GTF empty.
export STAR_GENOME_DIR="/oak/stanford/groups/emoding/scripts/rnaseq_scripts/index/STARIndex/gencode47_GRCh38"
export STAR_GTF=""

# STAR binary: use the locally-installed STAR 2.7.11b that your expression pipeline
# uses (and that built the index) — NOT the `star` module (2.7.10b), which would
# shadow it and refuse to load the genome.
export STAR_BIN_DIR="/oak/stanford/groups/emoding/analysis/brian/tools/bin"
export STAR_TWOPASS="false"      # 1-pass is enough for host removal & far faster (2-pass was the slow part)
export STAR_FFPE_RELAXED="true"  # relaxed match/score thresholds -> capture more short FFPE host reads
export STAR_SEED_LMAX=""         # empty = STAR default (50, fast). Set 25 for max sensitivity (much SLOWER).

# BBMap tools (bbduk.sh, clumpify.sh) for step-01 rRNA / low-complexity / dedup
# cleanup. `module load biology bbmap` does NOT expose them on this cluster, so
# install once and point BBMAP_BIN at the env's bin/ (empty = skip those steps):
#   conda create -p /oak/stanford/groups/emoding/analysis/brian/conda_envs/bbmap -c bioconda bbmap
#   then point BBMAP_BIN at the env's bin/ (it also contains java, so PATH alone works).
export BBMAP_BIN="/oak/stanford/groups/emoding/analysis/brian/conda_envs/bbmap/bin"

# Second-pass cleanup (recommended). Each is OPTIONAL and auto-skipped if the
# tool or reference is missing (the step logs what it actually ran).
#   - Human rRNA FASTA (e.g. SILVA human entries / NR_146144 etc.) for bbduk.
#   - Optional second host reference: T2T-CHM13 bowtie2 index prefix, which
#     traps residual human reads that GRCh38 misses.
# Human-only rRNA reference for the step-01 rRNA scrub. Build it once with
# build_human_rrna.sh (writes to this path). If absent, the rRNA sub-step skips.
export HUMAN_RRNA_FASTA="$RUN_DIR/refs/human_rRNA.fa"
export T2T_BT2_INDEX="$RUN_DIR/refs/chm13v2.0"   # build once with build_t2t_index.sh (B4 skips if absent)
export LOWCOMPLEXITY_ENTROPY="0.30"   # bbduk entropy filter (0 disables)
export DEDUP="true"                    # clumpify exact-duplicate (PCR) removal if available
                                       # (optical dedup is NOT enabled; add optical=t dupedist=N
                                       #  in step 01 if you need it for patterned flowcells)

# ---------------------------------------------------------------------------
# 4. KRAKEN2 + BRACKEN   (steps 02, 03)
# ---------------------------------------------------------------------------
export KRAKEN2_SIF="/oak/stanford/groups/emoding/analysis/ajayss/TaxonomyProfiling/kraken2/kraken2.sif"
export BRACKEN_SIF="/oak/stanford/groups/emoding/analysis/ajayss/TaxonomyProfiling/kraken2/bracken.sif"
# DB MUST contain the human genome (so host reads are trapped as Homo sapiens,
# not misassigned to microbes) and must be a complete build (taxo.k2d present).
export KRAKEN2_DB="/oak/stanford/groups/emoding/analysis/ajayss/TaxonomyProfiling/kraken2/kraken2_standard"

# Container runtime for the Kraken2/Bracken .sif images (steps 02 & 03). Detected
# at runtime -- apptainer preferred, singularity as a fallback -- so a node that
# only provides one of them still works. Steps must use "$CONTAINER_RUN", never a
# hardcoded `apptainer`, or preflight could pass on singularity while the job dies.
if command -v apptainer >/dev/null 2>&1;   then export CONTAINER_RUN="apptainer"
elif command -v singularity >/dev/null 2>&1; then export CONTAINER_RUN="singularity"
else export CONTAINER_RUN="apptainer"; fi   # fall back; preflight flags a genuine absence

# False-positive controls (the colleagues' runs used Kraken2 DEFAULTS -> ~99%
# of host-depleted reads "classified", which is a classic over-call).
export KRAKEN2_CONFIDENCE="0.1"       # require 10% of a read's k-mers to agree
export KRAKEN2_MIN_HIT_GROUPS="3"     # >=3 distinct minimizer groups per call
export REPORT_MINIMIZER_DATA="true"   # adds distinct-minimizer cols (KrakenUniq-style)

export READ_LEN="150"                 # NovaSeq 2x150; set to your actual length
export BRACKEN_THRESH="10"            # drop species below this many reads
export BRACKEN_LEVELS="S G"           # species and genus

# ---------------------------------------------------------------------------
# 5. KAIJU   (step 04)
# ---------------------------------------------------------------------------
export KAIJU_BIN="/home/groups/emoding/Ajay_Home/anaconda_ajay/envs/kaiju_env/bin"
export KAIJU_DB_DIR="/oak/stanford/groups/emoding/analysis/ajayss/TaxonomyProfiling/kaiju_databases"
export KAIJU_INDEX="$KAIJU_DB_DIR/kaiju_db_nr_euk.fmi"
export KAIJU_NODES="$KAIJU_DB_DIR/nodes.dmp"
export KAIJU_NAMES="$KAIJU_DB_DIR/names.dmp"
export KAIJU_MODE="greedy"
export KAIJU_MISMATCHES="3"           # -e 3

# ---------------------------------------------------------------------------
# 6. DOWNSTREAM FILTERING / DECONTAMINATION   (step 05)
# ---------------------------------------------------------------------------
export CONTAMINANT_BLOCKLIST="$PIPELINE_DIR/contaminants_blocklist.txt"
export MIN_READS="10"                  # min reads for a taxon to be kept in a sample
export MIN_RPM="1"                     # min reads per million INPUT read pairs (depth-normalized:
                                       # denominator = STAR "Number of input reads", i.e. the full
                                       # pre-host library, NOT the post-host microbial pool)
export MIN_PREVALENCE_FRAC="0.0"       # set >0 to require detection in a fraction of samples
export MIN_MINIMIZER_COVERAGE="0.01"   # Kraken2 distinct-minimizers / total k-mers floor
export MIN_DISTINCT_MINIMIZERS="10"    # min distinct minimizers to keep a taxon (KrakenUniq-style)

# OPTIONAL decontam (R package). Provide a 2-column TSV ( sample <TAB> type )
# where type is "sample" or "control"; if set, step 05 runs decontam against the
# negative controls. Leave empty to skip (a kitome blocklist is used instead).
export CONTROLS_LIST=""

# Concordance: a taxon is "high-confidence" if called by BOTH Kraken2 and Kaiju.
export REQUIRE_BOTH_CLASSIFIERS="true"
