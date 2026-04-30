# cfDNA-PBMC Somatic Mutation Calling Pipeline

This repository contains a high-confidence bioinformatics pipeline designed for calling somatic variants from paired cell-free DNA (cfDNA) and matched normal peripheral blood mononuclear cell (PBMC) samples. 

The pipeline is rigorously optimized to meet the **CAPP-Seq** (Cancer Personalized Profiling by deep Sequencing) publication standards for detecting ultra-low frequency variants (≤1% VAF), while minimizing false positives derived from sequencing and PCR artifacts.

---

## 🧬 Pipeline Architecture

The workflow integrates standard state-of-the-art tools, enhanced with custom artifact suppression modeling:

1. **SNV & Indel Calling (GATK Mutect2)**
   - Utilizes a **Panel of Normals (PoN)** to filter out recurrent germline artifacts and systematic sequencing errors.
   - Leverages **Read Orientation Modeling** (`LearnReadOrientationModel`) to eliminate strand-bias artifacts like 8-oxoguanine oxidation errors, critical for low-VAF cfDNA calling.
   - Restricts calling to the targeted capture regions via BED files to optimize runtime.

2. **Copy Number Variation (CNVkit)**
   - Operates in hybrid-capture mode utilizing a pooled normal reference baseline for stable copy number calculations.
   - Generates segment-level integer copy numbers and explicitly calculates statistical p-values using the binomial test (`cnvkit.py bintest`).

3. **Structural Variants (Manta)**
   - Executes somatic structural variant calling using matched tumor/normal pairs.

4. **Rigorous Filtering (R & CAPP-Seq Heuristics)**
   - Applies strict depth, VAF, and standard subtraction criteria against matched normal samples.
   - Validates multi-strand support for indels (`F1R2` and `F2R1` read parsing) to ensure true biological signal over PCR artifacts.
   - Cross-references COSMIC for biological rescue logic based on PoN frequencies.

---

## 📂 Directory Structure

```
├── README.md
└── script/
    └── mutation_calling/
        ├── run_full_pipeline.sh               # Master SLURM orchestrator
        ├── build_cnvkit_reference.sh          # Generates pooled normal CNV reference
        ├── build_mutect2_pon.sh               # Generates Mutect2 Panel of Normals (PoN)
        ├── mutation_calling_from_bam.sh       # Main per-sample calling script (SNVs, CNVs, SVs)
        ├── summarize_variant_calls.sh         # Aggregates VCF/CNS outputs into tidy TSV formats
        └── test.R                             # R script for final CAPP-Seq statistical filtering
```

> Note: All raw data, intermediate BAM files, VCFs, and PDF analysis outputs are explicitly ignored by Git to keep the repository lightweight.

---

## 🚀 Usage Guide

The pipeline is orchestrated to run on the **Stanford Sherlock HPC cluster** using the SLURM workload manager.

### 1. Configuration
Edit the configuration section at the top of `script/mutation_calling/run_full_pipeline.sh` to define your input batches, reference genome paths, target BED files, and directories.

### 2. Execution
Run the master orchestrator. It will automatically handle dependencies and submit jobs in the correct order:
1. Builds the CNVkit pooled reference
2. Builds the Mutect2 PoN
3. Executes per-sample variant calling arrays
4. Summarizes all cohort results

```bash
cd script/mutation_calling
bash run_full_pipeline.sh
```

### 3. Final Filtering & OncoPrint Generation
Once the summary jobs complete, run the R script to apply final CAPP-Seq heuristics and visualize the results.

```bash
Rscript script/mutation_calling/test.R
```

---

## ⚙️ Requirements
- SLURM Workload Manager
- Conda Environments: `cnvkit`, `snpeff`
- GATK4 (v4.6.0.0+)
- Samtools (v1.16.1+)
- bcftools (v1.16+)
- Python 2.7 (for Manta) / Python 3.x (for CNVkit)
- R (with `dplyr`, `tidyr`, `ComplexHeatmap` packages)
