# NLPHL tumor-microbiome pipeline (FFPE RNA-seq)

Detect and quantify microbial reads in the NLPHL FFPE RNA-seq cohort, starting
from the **demuxed + fastp-trimmed FASTQs** you already make for STAR-RSEM /
salmon. It runs **both** classifiers (Kraken2/Bracken and Kaiju), cross-validates
them, and applies the false-positive and contamination controls that the
exploratory scratch run lacked.

This rewrites and hardens the exploratory scripts in `../scratch/` (which were
prototyped on the public sarcoma set GSE213065).

---

## Workflow at a glance

Reads flow top-to-bottom; every `── drop ──>` is counted in the QC tables, so you
can see exactly how many reads or taxa leave at each gate. (`workflow.mermaid`
renders this as a diagram in the Cowork file viewer.)

```
 Trimmed FASTQs (192 samples)
        │
 STEP 01 host removal ── STAR maps human ───────────── drop ─> host reads
        │               ── bbduk host rRNA ──────────── drop ─> rRNA
        │               ── bbduk low-complexity ─────── drop ─> simple repeats
        │               ── clumpify dedup ───────────── drop ─> PCR/optical dups
        ▼
 Candidate-microbial reads ──────────────┐   (counts: QC_read_accounting.csv)
        │                                 │
 STEP 02 Kraken2 (DNA)            STEP 04 Kaiju (protein)
   conf 0.1 + hit-groups 3          greedy, -e 3
   (rest -> unclassified)                 │
        │                                 │
 STEP 03 Bracken (species+genus)   STEP 04b kaiju2table
        └──────────────┬──────────────────┘
                       ▼
 STEP 05 merge + filter + cross-validate
        ── remove human ───────────────── drop
        ── kitome blocklist ───────────── drop
        ── distinct-minimizer filter ──── drop   (counts: filter_funnel_summary.csv)
        ── abundance + prevalence ─────── drop
        ── Kraken2 ∩ Kaiju concordance ── drop
                       ▼
            HIGH-CONFIDENCE TAXA  ->  STEP 06 figures + pipeline_summary.txt
```

Read-level drops (host, rRNA, low-complexity, dups, Kraken2/Kaiju unclassified) land in
`QC_read_accounting.csv`; taxon-level drops (the discovery filters) land in
`filter_funnel_summary.csv`; `pipeline_summary.txt` narrates both in one page.

---

## TL;DR — how to run

```bash
cd script/taxonomic/pipeline
# 1. Edit config.sh  (paths marked "*** SET THESE ***", databases, read length)
nano config.sh
# 2. Submit everything (builds the sample sheet, chains all steps):
bash run_all.sh 8
# 3. Watch:
squeue -u $USER
```

Final tables land in `$RUN_DIR/05_results/`, figures in `$RUN_DIR/06_figures/`.

To run steps by hand instead of `run_all.sh`:

```bash
bash 00_build_sample_sheet.sh                       # makes logs/sample_list.tsv (N lines)
sbatch --array=1-N%8 01_host_removal.sbatch
sbatch --array=1-N%4 02_kraken2.sbatch              # after 01
sbatch            03_bracken_merge.sbatch           # after 02
sbatch --array=1-N%2 04_kaiju.sbatch                # after 01 (low concurrency!)
sbatch            04b_kaiju2table.sbatch            # after 04
sbatch            05_downstream.sbatch              # after 03 and 04b
```

---

## The steps

| step | script | what it does |
|------|--------|--------------|
| 00 | `00_build_sample_sheet.sh` | finds trimmed FASTQs → `sample_list.tsv` (`batch  sample  r1  r2`) |
| 01 | `01_host_removal.sbatch` | STAR → GRCh38; keeps **unmapped** reads; then rRNA / low-complexity / dup / optional T2T cleanup |
| 02 | `02_kraken2.sbatch` | Kraken2 with **confidence + min-hit-groups + minimizer data** |
| 03 | `03_bracken_merge.sbatch` | Bracken (species+genus), BIOM merge, KrakenUniq-style minimizer table, RPM lookup |
| 04 | `04_kaiju.sbatch` | Kaiju protein search vs `nr_euk` (greedy, -e 3) |
| 04b | `04b_kaiju2table.sbatch` | merge Kaiju into cohort matrices |
| 05 | `05_downstream.sbatch` → `05_aggregate_filter.R`, `06_visualize.R` | host/contaminant removal, filters, **Kraken2 ∩ Kaiju** cross-validation, optional decontam, plots |

---

## Why this differs from the scratch run (read this)

The scratch Kraken2 run reported **~99% of host-depleted reads as "classified."**
On tumor RNA that is a red flag, not a result — it is the same failure mode that
got the 2020 TCGA pan-cancer microbiome paper **retracted**. Fixes built in here:

- **Confidence threshold** `--confidence 0.1` and **`--minimum-hit-groups 3`**
  (Kraken2 default is 0 / 1 → maximal false positives).
- **Distinct-minimizer filter** (`--report-minimizer-data`): drop taxa supported
  by too few unique k-mers — a KrakenUniq-style guard against reads piling onto
  one conserved region.
- **Stronger host depletion**: STAR unmapped reads are further cleaned of host
  rRNA, low-complexity sequence, and PCR duplicates (and optionally residual
  human via T2T-CHM13). The Kraken2 DB must include the human genome so leftover
  host is trapped as *Homo sapiens* and then removed.
- **Two orthogonal classifiers**: a taxon is **high-confidence only if both**
  Kraken2 (nucleotide; sees rRNA) and Kaiju (protein; robust to FFPE damage)
  call it. Pico-kit libraries keep microbial rRNA, so Kraken2 will out-detect
  Kaiju for bacteria — concordance is the believable core.
- **Kitome blocklist** (`contaminants_blocklist.txt`, Salter 2014) and an
  optional **`decontam`** hook for true negative controls.

---

## Key knobs in `config.sh`

- **Paths you must set**: `TRIM_DIR`, `STAR_GENOME_DIR` (GRCh38/GENCODE v47 index
  — already filled in; GTF not needed, it is baked into the index), `KRAKEN2_DB`
  (must include human + be a complete build with `taxo.k2d`), `KAIJU_INDEX`,
  `READ_LEN`. Outputs already default to `…/Binkley_NLPHL/microbiome`.
- **STAR (speed)**: step 01 uses the local STAR 2.7.11b at `STAR_BIN_DIR` (not
  the `star` module) and adds `--outReadsUnmapped Fastx`. Defaults are tuned for
  speed: `STAR_TWOPASS=false` (1-pass is enough for host removal) and
  `STAR_SEED_LMAX=""` (STAR's default seed spacing — dense seeding `25` was the
  main slowdown). `STAR_FFPE_RELAXED=true` keeps relaxed match thresholds so
  short FFPE host reads are still captured. Set `STAR_TWOPASS=true` /
  `STAR_SEED_LMAX=25` only if you want maximum sensitivity at ~10x the runtime.
- **T2T second-host pass (B4)**: build once with `bash build_t2t_index.sh`
  (downloads the prebuilt CHM13v2.0 bowtie2 index). With `T2T_BT2_INDEX` set,
  step 01 removes residual human that GRCh38 misses; absent, B4 just skips.
- **FP stringency**: `KRAKEN2_CONFIDENCE` (0.1), `KRAKEN2_MIN_HIT_GROUPS` (3),
  `MIN_MINIMIZER_COVERAGE` (0.01), `MIN_DISTINCT_MINIMIZERS` (10).
- **Abundance/prevalence**: `MIN_READS` (10), `MIN_RPM` (1),
  `MIN_PREVALENCE_FRAC` (0 = off).
- **Decontam**: set `CONTROLS_LIST` to a `sample<TAB>type` TSV (`type` =
  `sample`|`control`) to enable. Leave empty to rely on the blocklist only.
- **Concordance**: `REQUIRE_BOTH_CLASSIFIERS=true` makes "high-confidence"
  require both Kraken2 and Kaiju.

---

## Outputs (`05_results/`)

- `cross_validated_species.csv` — every sample×species call with all metrics and
  pass/fail flags (host, blocklist, minimizer, abundance, prevalence,
  `detected_by`, `high_confidence`).
- `high_confidence_taxa.csv` — the filtered, believable subset.
- `cohort_high_confidence_summary.csv` — per-species prevalence and RPM summary.
- `pipeline_summary.txt` — **one-page plain-language run summary**: median read
  funnel, classifier rates, the taxon funnel, and the top high-confidence species.
  Read this first.
- `QC_read_accounting.csv` — per-sample read funnel: input → STAR-unmapped →
  after rRNA → after low-complexity → after dedup → final, **plus explicit
  `*_removed` counts/percentages per stage**, Kraken2
  `classified_pairs`/`unclassified_pairs`/`pct_unclassified`, Kaiju
  `pct_kaiju_classified`, and `n_high_confidence_species` (believable microbes
  found in that sample).
- `filter_funnel_summary.csv` — cohort-level count of sample-taxon calls
  surviving each gate (host → blocklist → minimizer → abundance → prevalence →
  classifier concordance → high-confidence), with `calls_dropped_here`.
- `decontam_contaminants.csv` — only if controls were supplied.

Figures (`06_figures/`): read-accounting QC (incl. per-stage removals and
Kraken2 unclassified rate), the filter funnel, cohort landscape, and a
Kraken2-vs-Kaiju concordance scatter.

---

## Dependencies

Faithful to the lab setup: STAR (local 2.7.11b at `STAR_BIN_DIR`) / samtools /
BBMap via `module load`; Kraken2 & Bracken via the Apptainer `.sif` images in
`config.sh`; Kaiju via the `kaiju_env` conda bin; **R 4.4.2** with the tidyverse
*core* packages (`readr`, `dplyr`, `tidyr`, `tibble`, `stringr`, `purrr`,
`ggplot2`), `scales`, `phyloseq`, `biomformat` (+ `decontam` if using controls).
The full `tidyverse` meta-package is deliberately avoided — it Imports
`ragg`/`googledrive`/`reprex`, which need system libs (cmake, freetype, libwebp)
absent on the login node. BBMap tools (`bbduk.sh`, `clumpify.sh`) and bowtie2 are
optional in step 01 — if absent, those cleanup sub-steps are skipped and logged.

**One-time R setup** (package libraries are R-version-specific, so 4.4.2 needs
its own install):

```bash
ml R/4.4.2
Rscript install_r_deps.R     # core tidyverse pkgs + phyloseq, biomformat, scales, decontam
```

`05_downstream.sbatch` loads `R/4.4.2` and fails fast with this exact pointer if a
package is missing. **Do not** add `module --force purge` to any step — it breaks
the Sherlock module environment.
