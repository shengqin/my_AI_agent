# Liquid Biopsy & Tumor RNA-seq Pipelines

HPC pipelines for **cell-free DNA (cfDNA) somatic variant calling** and **tumor-microbiome profiling from RNA-seq**, developed for Hodgkin lymphoma translational studies at Stanford University (Radiation Oncology).

| Pipeline | Input | What it produces |
|---|---|---|
| [`script/mutation_calling/`](script/mutation_calling) | Plasma cfDNA + matched PBMC BAMs (hybrid capture) | SNVs, indels, copy-number changes, structural variants, tumor fraction, cohort OncoPrint |
| [`script/taxonomic/pipeline/`](script/taxonomic/pipeline) | Trimmed RNA-seq FASTQs (FFPE tumors) | High-confidence microbial taxa with QC accounting at every filter |

Both pipelines run on SLURM (Stanford Sherlock), are restartable, and log every read or variant they drop.

---

## 1. cfDNA somatic variant calling

A tumor–normal pipeline for low-allele-fraction variants in plasma, built around published CAPP-Seq filtering heuristics.

### Workflow

```
 Matched PBMC BAMs ──► Mutect2 Panel of Normals  ─┐
                  └──► CNVkit pooled reference  ─┤      (built once per cohort)
                                                 ▼
 cfDNA + PBMC BAMs ──► per-sample SLURM array
                        ├─ Mutect2 (tumor–normal, PoN, 0.1% AF floor)
                        │    └─ LearnReadOrientationModel → FilterMutectCalls
                        ├─ bcftools norm → snpEff (canonical transcripts)
                        ├─ CNVkit (hybrid) → call → bintest → genemetrics
                        └─ Manta (somatic SVs)
                                                 ▼
                       cohort summarizer → tidy TSVs
                                                 ▼
                 R filtering + reporting (CAPP-Seq heuristics)
```

### Filtering and reporting (R stage)
- Depth, alt-read support, and VAF thresholds, with strict limits on allele fraction in the matched normal
- Strand support for indels (F1R2 and F2R1) and a homopolymer/STR slippage filter
- Panel-of-normals frequency filter for recurrent background noise
- COSMIC-based rescue of known drivers near the detection threshold
- **CHIP genes flagged** and excluded from tumor-fraction estimates, so clonal hematopoiesis isn't called as tumor
- Tumor fraction estimated from confident SNVs, with tumor-fraction-aware CNV calls
- Cohort **OncoPrint** of SNVs, indels, CNVs, and SVs (ComplexHeatmap)

### Run

```bash
cd script/mutation_calling
# 1. Edit the configuration block at the top of run_full_pipeline.sh
#    (batches, normals, target BED, reference genome, output directories)
bash run_full_pipeline.sh                      # build references + call + summarize
SKIP_PON_BUILD=true bash run_full_pipeline.sh  # reuse existing PoN / CNV reference

# 2. Filter, estimate tumor fraction, and plot
Rscript test.R

# Optional: re-call CNVs with per-sample purity
bash recall_cnvkit_with_purity.sh
```

The orchestrator submits jobs with SLURM dependencies: (1) CNVkit reference and (2) Mutect2 PoN in parallel, then (3) per-sample calling arrays, then (4) the cohort summary. Each step checks for complete outputs and skips them, so you can safely re-run.

### Key scripts

| Script | Purpose |
|---|---|
| `run_full_pipeline.sh` | Single-command SLURM orchestrator |
| `build_mutect2_pon.sh` | Mutect2 Panel of Normals from matched normals |
| `build_cnvkit_reference.sh` | Pooled-normal CNVkit reference |
| `mutation_calling_from_bam.sh` | Per-sample SNV/indel, CNV, and SV calling |
| `summarize_variant_calls.sh` | Flattens per-sample VCF/CNVkit/Manta outputs into cohort TSVs |
| `test.R` | CAPP-Seq-style filtering, tumor fraction, OncoPrint |
| `recall_cnvkit_with_purity.sh` | Optional purity-aware CNV re-calling |
| `build_panel_gene_coords.py` | Gene spans for the capture panel (UCSC refGene) for SV-to-gene mapping |

### Requirements
SLURM · GATK 4.6 · Manta 1.6 · CNVkit · snpEff · bcftools / samtools · R (tidyverse, ComplexHeatmap, circlize) · Python 3

---

## 2. Tumor-microbiome profiling from RNA-seq

Detects and quantifies microbial reads in an FFPE tumor RNA-seq cohort. It runs **two independent classifiers (Kraken2/Bracken and Kaiju)** and keeps only taxa they agree on, after host, contamination, and abundance filters.

- Host removal (STAR, host rRNA, low-complexity, and duplicate filters)
- Kraken2 → Bracken and Kaiju → kaiju2table in parallel
- Cross-validation by NCBI taxid or genus-aware name matching
- Kitome/contaminant blocklist, minimizer, abundance, prevalence, and concordance filters
- Full accounting: `QC_read_accounting.csv`, `filter_funnel_summary.csv`, and `rejected_calls.csv.gz`, where every rejected call carries its reason

See [`script/taxonomic/pipeline/README.md`](script/taxonomic/pipeline/README.md) for the full workflow diagram and run instructions.

---

## Repository layout

```
script/
├── mutation_calling/      # cfDNA somatic pipeline (see above)
├── taxonomic/pipeline/    # tumor-microbiome pipeline (own README)
├── cnv_calling_from_bam.sh   # earlier standalone CNVkit workflow (superseded by mutation_calling/)
├── cHL.R                  # exploratory cohort analysis of SNV calls
├── debug.R                # scratch
└── ConsensusExomePipelineHuman.sh   # third-party reference (see note)
```

Raw data, BAMs, VCFs, and figures are excluded by `.gitignore`. **No patient data is stored in this repository.**

**Note on attribution:** `ConsensusExomePipelineHuman.sh` is a colleague's exome consensus-calling pipeline, kept here for reference only. I did not write it.

**Paths:** scripts reference Stanford Sherlock storage (`/oak/...`). Edit the configuration blocks to run them elsewhere.

---

## Author
**Shengqin (Brian) Su, PhD**, Research Scientist, Stanford University
[shengqin.github.io](https://shengqin.github.io) · [LinkedIn](https://www.linkedin.com/in/shengqinsu/)
