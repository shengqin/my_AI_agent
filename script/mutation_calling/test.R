# ==============================================================================
# Cohort OncoPrint Generator for Somatic Variant Calls
# Reads from combined_summary/ mutect2, cnvkit, and manta outputs.
# Produces an OncoPrint with SNVs, Indels, CNVs, and Structural Variants.
# ==============================================================================

# 1. Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Define paths (absolute to avoid working directory issues)
PROJECT_DIR <- "/Users/Brian/Library/CloudStorage/Box-Box/BrianAnalysis/Lymphoma/NLPHL"
DATA_DIR <- file.path(PROJECT_DIR, "data")
RESULTS_DIR <- file.path(DATA_DIR, "combined_summary")
MUTECT_FILE <- file.path(RESULTS_DIR, "mutect2_variants.tsv")
CNVKIT_FILE <- file.path(RESULTS_DIR, "cnvkit_aberrations.tsv")
CNV_GENE_FILE <- file.path(RESULTS_DIR, "cnvkit_gene_cnvs.tsv")
SV_FILE <- file.path(RESULTS_DIR, "manta_somatic_svs.tsv")
COHORT_VARIANT_SUMMARY_FILE <- file.path(RESULTS_DIR, "cohort_variant_summary.csv")
OUTPUT_PDF <- file.path(PROJECT_DIR, "analysis", "cohort_oncoprint.pdf")

# PoN frequency summary from build_mutect2_pon.sh (set to NA to skip PoN filtering)
# Lives in combined_summary/ as written by summarize_variant_calls.sh on Sherlock.
PON_SUMMARY_FILE <- file.path(RESULTS_DIR, "pon_site_summary.tsv")

# COSMIC coding mutations VCF for rescue logic (set to NA to skip COSMIC rescue)
# Download from: https://cancer.sanger.ac.uk/cosmic/download (requires free registration)
# Example: CosmicCodingMuts.vcf.gz (hg19/GRCh37)
COSMIC_VCF <- file.path(DATA_DIR, "Cosmic_CompleteTargetedScreensMutant_v103_GRCh37.vcf.gz")

# Aberrant somatic hypermutation (aSHM) hotspots in germinal-center B-cell lymphoma.
# Mutations here are predominantly aSHM by-products (passengers). They are (a) FLAGGED and shown
# distinctly as "aSHM-associated (likely passenger)", and (b) ALLOWED TO BYPASS the panel-of-normals
# recurrence filter — recurrence at these loci reflects aSHM biology, not a technical artifact.
# Per-sample normal subtraction (QC Filter 2) STILL applies, so germline/CHIP leakage is guarded.
# Curated/editable. Bona fide point-mutation drivers (SGK1, DUSP2, SOCS1, STAT6, GNA13, CREBBP,
# EP300, EZH2, TNFAIP3, B2M, CARD11, MYD88, TP53, KMT2D, ...) are intentionally EXCLUDED.
# NOTE: MYC/BCL2/BCL6 are listed because their recurrent *coding point mutations* in GC lymphoma
# are largely aSHM; their driver role is via translocation (captured as SVs), not point mutations.
ASHM_GENES <- c(
  "IGLL5", "BCL2", "BCL6", "MYC", "PAX5", "PIM1", "RHOH", "IRF4", "BCL7A",
  "KLHL6", "SERPINA9", "ST6GAL1", "DTX1", "ZFP36L1", "CD83", "BTG1", "BTG2",
  "H1-2", "H1-3", "H1-4", "H1-5", "H3C2"
)

# Clonal hematopoiesis (CHIP) genes. Because the matched normal is PBL (blood) and the
# tumor input is plasma cfDNA, mutations in these genes may originate from clonal
# hematopoiesis rather than lymphoma. They are KEPT but tagged for manual adjudication.
CHIP_GENES <- c("DNMT3A", "TET2", "ASXL1", "PPM1D", "TP53", "JAK2",
                "SF3B1", "SRSF2", "GNB1", "CBL", "GNAS")

# Canonical NLPHL/GCB-lymphoma drivers eligible for a narrow near-threshold rescue.
# Intentionally excludes aSHM-passenger genes and CHIP genes at filtering time below.
DRIVER_GENES <- c("SGK1", "DUSP2", "JUNB", "SOCS1", "NFKBIE", "STAT6", "GNA13",
                  "CREBBP", "EP300", "EZH2", "TNFAIP3", "B2M", "CIITA", "BCL6",
                  "ITPKB", "MEF2B", "KMT2D", "CARD11", "XPO1", "ARID1A", "ATM",
                  "TP53", "NOTCH1", "NOTCH2", "SPEN", "HIST1H1E", "GNA13",
                  "MYD88", "CD79B", "PIM1")

# Ensure the output directory exists (for OncoPrint PDF, QC plots, and review CSVs below).
dir.create(file.path(PROJECT_DIR, "analysis"), showWarnings = FALSE, recursive = TRUE)

cat("Loading and processing variant data...\n")

# ==============================================================================
# 2. Read and QC SNV/Indel Data (Mutect2 + SnpEff)
# ==============================================================================
# Read data
snv_data <- read_tsv(MUTECT_FILE, col_types = cols(Tumor_F1R2 = col_character(), Tumor_F2R1 = col_character(), Tumor_AF = col_character(), Normal_AF = col_character(), Tumor_AD_Alt = col_character(), Tumor_AD_Ref = col_character()), show_col_types = FALSE)

sample_ids_from_file <- function(file, reader) {
  if (!file.exists(file)) return(character(0))
  tryCatch({
    dat <- reader(file)
    if ("Sample" %in% colnames(dat)) as.character(dat$Sample) else character(0)
  }, error = function(e) character(0))
}

sample_ids_from_object_or_file <- function(object_name, file, reader) {
  if (exists(object_name)) {
    dat <- get(object_name)
    if ("Sample" %in% colnames(dat)) return(as.character(dat$Sample))
  }
  sample_ids_from_file(file, reader)
}

all_cohort_samples <- sort(unique(na.omit(c(
  as.character(snv_data$Sample),
  sample_ids_from_object_or_file("cnv_seg_data", CNVKIT_FILE, function(x) read_tsv(x, show_col_types = FALSE)),
  sample_ids_from_object_or_file("sv_data", SV_FILE, function(x) read_tsv(x, show_col_types = FALSE)),
  sample_ids_from_file(COHORT_VARIANT_SUMMARY_FILE, function(x) read_csv(x, show_col_types = FALSE))
))))

# QC Filter 1: Only keep high-confidence PASS variants
snv_filt <- snv_data %>%
  filter(Filter == "PASS") %>%
  # Remove known false-positive / mapping-artifact / pseudogene genes.
  # NOTE: IGLL5 is intentionally NOT dropped here — it is a real aSHM target in B-cell
  # lymphoma, so it is retained and flagged as likely-passenger below (see ASHM_GENES).
  filter(!Gene %in% c("ACTB", "ABO", "MPEG1", "ELK2AP", "ELK2AP-MIR4507", "MIR8071-1-ELK2AP"))

# QC Filter 2: Enforce the CAPP-Seq Published Strict Heuristics
snv_filt <- snv_filt %>%
  # Handle scientific notation (e.g. 4.5e-03) and multiallelic comma separation
  mutate(
    Tumor_AF_num = as.numeric(sapply(strsplit(as.character(Tumor_AF), ","), `[`, 1)),
    Normal_AF_num = as.numeric(sapply(strsplit(as.character(Normal_AF), ","), `[`, 1)),
    Tumor_AD_Alt_num = as.numeric(sapply(strsplit(as.character(Tumor_AD_Alt), ","), `[`, 1))
  ) %>%
  mutate(Normal_AF_num = replace_na(Normal_AF_num, 0)) %>%
  filter(
    Tumor_DP >= 100, # Deduplicated depth of 100 or more
    Tumor_AD_Alt_num >= 4, # 4 or more variant supporting reads
    Tumor_AF_num >= 0.005, # VAF >= 0.5% (genotyping threshold)
    Normal_DP >= 20, # Normal deduplicated depth of 20 or more
    # Normal VAF strict subtraction rule:
    # < 0.25%, OR < 1% provided sample VAF is > 20-fold normal VAF
    (Normal_AF_num < 0.0025) |
      (Normal_AF_num < 0.01 & Tumor_AF_num > (20 * Normal_AF_num))
  )

# QC Filter 2.5: Indel-specific Strand Bias Filter
# CAPP-Seq requires indels to have read support from both forward and reverse orientations
# to eliminate PCR and sequencing artifacts (especially at <1% VAF)
# Backward compatible: only applies if F1R2/F2R1 columns exist in the data.
has_strand_cols <- all(c("Tumor_F1R2", "Tumor_F2R1") %in% colnames(snv_filt))

if (has_strand_cols) {
  snv_filt <- snv_filt %>%
    mutate(
      F1R2_Alt = as.numeric(sapply(strsplit(as.character(Tumor_F1R2), ","), function(x) if(length(x)>1) x[2] else "0")),
      F2R1_Alt = as.numeric(sapply(strsplit(as.character(Tumor_F2R1), ","), function(x) if(length(x)>1) x[2] else "0")),
      Is_Indel = (nchar(Ref) != nchar(Alt)) | nchar(Ref) > 1 | nchar(Alt) > 1
    ) %>%
    filter(
      # If it's not an indel, keep it. If it is an indel, require both strands > 0
      !Is_Indel | (F1R2_Alt > 0 & F2R1_Alt > 0)
    )
  cat("  Indel strand bias filter: applied (F1R2/F2R1 columns present)\n")
} else {
  cat("  Indel strand bias filter: SKIPPED (F1R2/F2R1 columns not found — re-run pipeline with orientation model)\n")
}

# QC Filter 2b: PoN Frequency Filter (per published CAPP-Seq heuristics)
# Variant must have read support >= 0.1% VAF in:
#   (1) < 10% of PoN samples, OR
#   (2) < 30% of PoN samples, provided sample VAF > 20-fold max PoN VAF
# When variant IS seen in any PoN sample: apply stricter normal VAF (0.1%)
#   unless tumor VAF > 20× max(normal VAF, max PoN VAF)
if (!is.na(PON_SUMMARY_FILE) && file.exists(PON_SUMMARY_FILE)) {
  cat("  Applying PoN frequency filter...\n")
  pon_data <- read_tsv(PON_SUMMARY_FILE, show_col_types = FALSE)

  n_before_pon <- nrow(snv_filt)

  # Join variants with PoN summary on coordinate + alleles
  snv_filt <- snv_filt %>%
    left_join(
      pon_data,
      by = c("Chromosome" = "Chromosome", "Position" = "Position",
             "Ref" = "Ref", "Alt" = "Alt")
    ) %>%
    # Fill NA = variant not seen in any PoN sample (clean)
    mutate(
      N_PON_Samples = replace_na(N_PON_Samples, 0),
      Max_PON_AF = replace_na(Max_PON_AF, 0),
      Fraction_PON_Samples = replace_na(Fraction_PON_Samples, 0),
      In_PON = N_PON_Samples > 0
    )

  # Save the FULL set before PoN filtering — needed for COSMIC rescue below
  snv_pre_pon <- snv_filt

  # Apply PoN frequency filter:
  # Keep if: not in PoN, or in < 10% of PoN, or in < 30% with tumor >> PoN
  snv_filt <- snv_filt %>%
    filter(
      # aSHM hotspots bypass PoN recurrence: recurrence here is aSHM biology, not artifact
      # (per-sample normal subtraction from QC Filter 2 still applied above).
      Gene %in% ASHM_GENES |
      # Not in PoN at all — automatically passes
      Fraction_PON_Samples == 0 |
      # In < 10% of PoN samples
      Fraction_PON_Samples < 0.10 |
      # In < 30% of PoN samples AND tumor VAF > 20× max PoN VAF
      (Fraction_PON_Samples < 0.30 & Tumor_AF_num > (20 * Max_PON_AF))
    )

  # When variant IS seen in any PoN sample, apply stricter normal VAF threshold
  # Normal VAF must be < 0.1% (instead of 0.25%),
  # UNLESS tumor VAF > 20× max(normal VAF, max PoN VAF)
  #   AND read support in < 30% of PoN samples (per paper's additional constraint)
  snv_filt <- snv_filt %>%
    filter(
      # aSHM hotspots are exempt from the stricter in-PoN normal-VAF gate (biology, not artifact);
      # they already passed the 0.25% normal subtraction in QC Filter 2.
      Gene %in% ASHM_GENES |
      !In_PON |
      (Normal_AF_num < 0.001) |
      (Tumor_AF_num > 20 * pmax(Normal_AF_num, Max_PON_AF) &
         Fraction_PON_Samples < 0.30)
    )

  n_pon_removed <- n_before_pon - nrow(snv_filt)
  cat(sprintf("  PoN filter: %d -> %d variants (%d removed)\n",
              n_before_pon, nrow(snv_filt), n_pon_removed))
} else {
  cat("  PoN summary not found — skipping PoN frequency filter.\n")
  cat("  To enable: run build_mutect2_pon.sh and set PON_SUMMARY_FILE path above.\n")
  snv_filt <- snv_filt %>%
    mutate(In_PON = FALSE, Max_PON_AF = 0, Fraction_PON_Samples = 0)
  snv_pre_pon <- snv_filt  # No PoN filter applied, so pre = post
}

# QC Filter 2c: COSMIC Rescue
# Per the paper: "COSMIC-annotated variants with sample VAF > 20-fold normal
# VAF were rescued from PON filtering steps."
# This re-adds PoN-filtered variants that are known oncogenic hotspots,
# preventing the loss of true somatic mutations that happen to appear
# at low levels in normal/PoN samples.
if (!is.na(COSMIC_VCF) && file.exists(COSMIC_VCF)) {
  cat("  Loading COSMIC annotations for rescue logic...\n")

  # Stream the COSMIC VCF instead of readLines() of the whole file: COSMIC coding VCFs are
  # large, and readLines loads every line into memory. We decompress, drop the header, and
  # cut to the first 5 columns in the shell pipe so R only ever sees chr/pos/id/ref/alt.
  cosmic_df <- tryCatch(
    read_tsv(
      pipe(sprintf("gzip -dc %s | grep -v '^#' | cut -f1-5", shQuote(COSMIC_VCF))),
      col_names = FALSE, show_col_types = FALSE, progress = FALSE
    ),
    error = function(e) tibble()
  )

  if (nrow(cosmic_df) > 0) {
    cosmic_df <- cosmic_df %>%
      rename(Chromosome = X1, Position = X2, COSMIC_ID = X3, Ref = X4, Alt = X5) %>%
      mutate(Position = as.numeric(Position)) %>%
      separate_rows(Alt, sep = ",") %>%
      # Normalize chromosome naming to match the pipeline's hg19 "chr"-prefixed contigs.
      # COSMIC GRCh37 VCFs use "1"/"MT" (no "chr"); without this the join silently matches
      # NOTHING and the COSMIC rescue becomes a no-op.
      mutate(
        Chromosome = ifelse(str_detect(Chromosome, "^chr"), Chromosome, paste0("chr", Chromosome)),
        Chromosome = ifelse(Chromosome == "chrMT", "chrM", Chromosome),
        Ref = toupper(Ref),
        Alt = toupper(Alt)
      ) %>%
      distinct(Chromosome, Position, Ref, Alt, .keep_all = TRUE)

    # Sanity check: how many of our variant positions actually match COSMIC after
    # normalization? A zero here means the join keys still disagree (investigate naming).
    n_cosmic_overlap <- snv_filt %>%
      inner_join(cosmic_df %>% select(Chromosome, Position, Ref, Alt),
                 by = c("Chromosome", "Position", "Ref", "Alt")) %>%
      nrow()
    cat(sprintf("  COSMIC join: %d of %d surviving variants match a COSMIC site (post chr-normalization)\n",
                n_cosmic_overlap, nrow(snv_filt)))

    # --- Step 1: Tag surviving variants with COSMIC status ---
    snv_filt <- snv_filt %>%
      left_join(
        cosmic_df %>% select(Chromosome, Position, Ref, Alt, COSMIC_ID),
        by = c("Chromosome", "Position", "Ref", "Alt")
      ) %>%
      mutate(Is_COSMIC = !is.na(COSMIC_ID), COSMIC_Rescued = FALSE)

    # --- Step 2: Identify PoN-removed variants eligible for COSMIC rescue ---
    # These are variants that were in snv_pre_pon but NOT in snv_filt
    snv_pon_removed <- snv_pre_pon %>%
      anti_join(snv_filt, by = c("Chromosome", "Position", "Ref", "Alt", "Sample"))

    if (nrow(snv_pon_removed) > 0) {
      # Annotate removed variants with COSMIC
      snv_pon_removed <- snv_pon_removed %>%
        left_join(
          cosmic_df %>% select(Chromosome, Position, Ref, Alt, COSMIC_ID),
          by = c("Chromosome", "Position", "Ref", "Alt")
        ) %>%
        mutate(Is_COSMIC = !is.na(COSMIC_ID))

      # Rescue condition: COSMIC-annotated AND tumor VAF > 20× normal VAF
      rescued <- snv_pon_removed %>%
        filter(
          Is_COSMIC,
          Tumor_AF_num > (20 * Normal_AF_num)
        ) %>%
        mutate(COSMIC_Rescued = TRUE)

      n_rescued <- nrow(rescued)

      if (n_rescued > 0) {
        # Add rescued variants back to the filtered set
        snv_filt <- bind_rows(snv_filt, rescued)
        cat(sprintf("  COSMIC rescue: %d variants rescued from PoN filter\n", n_rescued))
        if (n_rescued <= 20) {
          cat("  Rescued variants:\n")
          rescued %>%
            select(Sample, Gene, Chromosome, Position, Ref, Alt,
                   Tumor_AF_num, Normal_AF_num, COSMIC_ID) %>%
            print(n = Inf)
        }
      } else {
        cat("  COSMIC rescue: 0 variants eligible for rescue\n")
      }
    } else {
      cat("  COSMIC rescue: no variants were removed by PoN filter\n")
    }

    n_cosmic <- sum(snv_filt$Is_COSMIC)
    cat(sprintf("  COSMIC variants in final set: %d\n", n_cosmic))
  } else {
    snv_filt <- snv_filt %>%
      mutate(Is_COSMIC = FALSE, COSMIC_ID = NA_character_, COSMIC_Rescued = FALSE)
  }
} else {
  cat("  COSMIC VCF not configured — skipping COSMIC rescue.\n")
  cat("  To enable: download CosmicCodingMuts.vcf.gz and set COSMIC_VCF path above.\n")
  snv_filt <- snv_filt %>%
    mutate(Is_COSMIC = FALSE, COSMIC_ID = NA_character_, COSMIC_Rescued = FALSE)
}

# QC Filter 3: Remove synonymous/silent or unannotated variants
# Keep only protein-altering variants
altering_effects <- c(
  "missense_variant", "stop_gained", "frameshift_variant",
  "splice_acceptor_variant", "splice_donor_variant",
  "splice_acceptor_variant&intron_variant", "splice_donor_variant&intron_variant",
  "inframe_deletion", "inframe_insertion", "disruptive_inframe_deletion",
  "disruptive_inframe_insertion",
  "conservative_inframe_deletion", "conservative_inframe_insertion",
  "start_lost",
  "stop_lost", "5_prime_UTR_premature_start_codon_gain_variant",
  "missense_variant&splice_region_variant"
)

snv_qc <- snv_filt %>%
  filter(Variant_Type %in% altering_effects & Gene != ".") %>%
  # Classify as SNV or Indel based on Ref/Alt lengths
  mutate(
    Is_Indel = (nchar(Ref) != nchar(Alt)) | nchar(Ref) > 1 | nchar(Alt) > 1,
    Mutation_Class = case_when(
      # Indels (frameshift, inframe)
      str_detect(Variant_Type, "frameshift") ~ "Frameshift_Indel",
      str_detect(Variant_Type, "inframe") ~ "Inframe_Indel",
      # SNVs
      Variant_Type == "missense_variant" ~ "Missense_SNV",
      Variant_Type == "missense_variant&splice_region_variant" ~ "Missense_SNV",
      Variant_Type == "stop_gained" ~ "Nonsense_SNV",
      str_detect(Variant_Type, "splice_acceptor|splice_donor") ~ "Splice_Site",
      Variant_Type == "5_prime_UTR_premature_start_codon_gain_variant" ~ "UTR_SNV",
      TRUE ~ "Other_SNV"
    )
  ) %>%
  # Compute aSHM/CHIP gene flags for INTERNAL use only (the PoN-recurrence bypass applied
  # upstream, tumor-fraction exclusion, and CHIP review export). Per user request, aSHM
  # variants are shown as ORDINARY mutations (their real SNV/indel class) in the OncoPrint —
  # NO separate aSHM tag/category.
  mutate(
    Is_aSHM = Gene %in% ASHM_GENES,
    Is_CHIP = Gene %in% CHIP_GENES
  )

# Driver-gene near-threshold rescue: strictly limited to PASS, protein-altering,
# non-CHIP/non-aSHM drivers in the 0.25%-<0.5% VAF band. This recovers plausible
# low-VAF somatic drivers without bypassing Mutect2 filters, normal subtraction, or PoN checks.
snv_rescue <- snv_data %>%
  filter(Filter == "PASS") %>%
  filter(
    Gene %in% DRIVER_GENES,
    !(Gene %in% CHIP_GENES),
    !(Gene %in% ASHM_GENES),
    !Gene %in% c("ACTB", "ABO", "MPEG1"),
    !str_detect(Gene, "^ELK2AP")
  ) %>%
  mutate(
    Tumor_AF_num = as.numeric(sapply(strsplit(as.character(Tumor_AF), ","), `[`, 1)),
    Normal_AF_num = as.numeric(sapply(strsplit(as.character(Normal_AF), ","), `[`, 1)),
    Tumor_AD_Alt_num = as.numeric(sapply(strsplit(as.character(Tumor_AD_Alt), ","), `[`, 1))
  ) %>%
  mutate(Normal_AF_num = replace_na(Normal_AF_num, 0)) %>%
  filter(
    Variant_Type %in% altering_effects,
    Tumor_DP >= 100,
    Normal_DP >= 20,
    Tumor_AD_Alt_num >= 4,
    Tumor_AF_num >= 0.0025,
    Tumor_AF_num < 0.005,
    (Normal_AF_num < 0.0025) |
      (Normal_AF_num < 0.01 & Tumor_AF_num > (20 * Normal_AF_num))
  )

if (exists("pon_data")) {
  snv_rescue <- snv_rescue %>%
    left_join(
      pon_data,
      by = c("Chromosome" = "Chromosome", "Position" = "Position",
             "Ref" = "Ref", "Alt" = "Alt")
    ) %>%
    mutate(Fraction_PON_Samples = replace_na(Fraction_PON_Samples, 0)) %>%
    filter(Fraction_PON_Samples < 0.10)
}

if (all(c("Tumor_F1R2", "Tumor_F2R1") %in% colnames(snv_rescue))) {
  snv_rescue <- snv_rescue %>%
    mutate(
      F1R2_Alt = as.numeric(sapply(strsplit(as.character(Tumor_F1R2), ","), function(x) if(length(x)>1) x[2] else "0")),
      F2R1_Alt = as.numeric(sapply(strsplit(as.character(Tumor_F2R1), ","), function(x) if(length(x)>1) x[2] else "0")),
      Is_Indel = (nchar(Ref) != nchar(Alt)) | nchar(Ref) > 1 | nchar(Alt) > 1
    ) %>%
    filter(!Is_Indel | (F1R2_Alt > 0 & F2R1_Alt > 0))
}

snv_rescue <- snv_rescue %>%
  anti_join(snv_qc, by = c("Sample", "Chromosome", "Position", "Ref", "Alt")) %>%
  mutate(
    Is_Indel = (nchar(Ref) != nchar(Alt)) | nchar(Ref) > 1 | nchar(Alt) > 1,
    Is_aSHM = FALSE,
    Is_CHIP = FALSE,
    # Per user request, rescued driver variants are shown as ORDINARY mutations (their real
    # SNV/indel class) — NO separate "rescued" tag/category. The rescue criteria above are
    # unchanged; only the display label differs.
    Mutation_Class = case_when(
      str_detect(Variant_Type, "frameshift") ~ "Frameshift_Indel",
      str_detect(Variant_Type, "inframe") ~ "Inframe_Indel",
      Variant_Type == "missense_variant" ~ "Missense_SNV",
      Variant_Type == "missense_variant&splice_region_variant" ~ "Missense_SNV",
      Variant_Type == "stop_gained" ~ "Nonsense_SNV",
      str_detect(Variant_Type, "splice_acceptor|splice_donor") ~ "Splice_Site",
      Variant_Type == "5_prime_UTR_premature_start_codon_gain_variant" ~ "UTR_SNV",
      TRUE ~ "Other_SNV"
    )
  )

cat(sprintf("  Driver-gene near-threshold rescue: %d variants rescued\n", nrow(snv_rescue)))
if (nrow(snv_rescue) > 0) {
  snv_rescue %>%
    select(Sample, Gene, Tumor_AF = Tumor_AF_num) %>%
    arrange(Sample, Gene) %>%
    print(n = Inf)
  snv_qc <- bind_rows(snv_qc, snv_rescue)
}

cat(sprintf(
  "  Mutect2 after QC: %d protein-altering variants across %d genes\n",
  nrow(snv_qc), n_distinct(snv_qc$Gene)
))
cat(sprintf(
  "    SNVs: %d | Indels: %d | aSHM-gene (shown as normal mutations): %d\n",
  sum(snv_qc$Mutation_Class %in% c("Missense_SNV", "Nonsense_SNV", "Splice_Site", "UTR_SNV", "Other_SNV")),
  sum(snv_qc$Mutation_Class %in% c("Frameshift_Indel", "Inframe_Indel")),
  sum(snv_qc$Is_aSHM)
))

# (4) Tag CHIP-gene calls and export them for manual adjudication. They remain in the
# OncoPrint, but a plasma cfDNA mutation in these genes may be clonal hematopoiesis (blood)
# rather than lymphoma — especially where it is also detectable in the matched PBL normal.
chip_flagged <- snv_qc %>% filter(Is_CHIP)
if (nrow(chip_flagged) > 0) {
  CHIP_CSV <- file.path(PROJECT_DIR, "analysis", "chip_flagged_variants.csv")
  chip_flagged %>%
    select(any_of(c("Sample", "Gene", "Chromosome", "Position", "Ref", "Alt",
                    "Variant_Type", "Mutation_Class", "Tumor_AF_num", "Normal_AF_num",
                    "Tumor_DP", "Normal_DP", "Is_COSMIC", "COSMIC_ID"))) %>%
    arrange(Gene, Sample) %>%
    write_csv(CHIP_CSV)
  cat(sprintf("  CHIP-gene calls flagged for manual review: %d (genes: %s) -> %s\n",
              nrow(chip_flagged), paste(sort(unique(chip_flagged$Gene)), collapse = ", "), CHIP_CSV))
} else {
  cat("  CHIP-gene calls flagged for manual review: 0\n")
}

# First-pass tumor-fraction estimate from confident clonal SNVs.
# Source = snv_filt: the full confident somatic SNV set AFTER depth/alt/VAF QC, normal
# background subtraction, the panel-of-normals filter, and COSMIC rescue — but BEFORE the
# protein-altering functional whitelist. Using ALL confident SNVs (synonymous/non-coding
# passengers included), rather than only driver-level protein-altering calls, gives far more
# data points per sample and less driver-selection bias, so most samples become estimable.
# CHIP/aSHM genes and indels are excluded (computed inline since snv_filt predates those flags).
# Assumptions/limitations: 2 x median SNV VAF under a heterozygous-diploid model, so it is a
# lower-bound-ish estimate; subclonality biases it down and CNV-altered regions can confound it.
tumor_fraction_counts <- snv_filt %>%
  filter(
    !(Gene %in% CHIP_GENES),
    !(Gene %in% ASHM_GENES),
    nchar(Ref) == 1 & nchar(Alt) == 1, # SNVs only (exclude indels)
    !is.na(Tumor_AF_num)
  ) %>%
  group_by(Sample) %>%
  summarize(
    N_SNV_for_TF = n(),
    Median_Tumor_AF = median(Tumor_AF_num),
    .groups = "drop"
  )

# Sample universe = every sample seen in SNV, CNV, SV, or cohort summary inputs.
tumor_fraction <- tibble(Sample = all_cohort_samples) %>%
  left_join(tumor_fraction_counts, by = "Sample") %>%
  mutate(
    N_SNV_for_TF = replace_na(N_SNV_for_TF, 0L),
    Tumor_Fraction = if_else(N_SNV_for_TF >= 3, pmin(1, 2 * Median_Tumor_AF), NA_real_),
    TF_Method = "2x_median_confident_SNV_VAF"
  ) %>%
  select(Sample, N_SNV_for_TF, Tumor_Fraction, TF_Method)

TUMOR_FRACTION_CSV <- file.path(PROJECT_DIR, "analysis", "tumor_fraction.csv")
write_csv(tumor_fraction, TUMOR_FRACTION_CSV)
cat(sprintf(
  "  Tumor-fraction estimates: %d samples (%d with >=3 confident SNVs) -> %s\n",
  nrow(tumor_fraction), sum(!is.na(tumor_fraction$Tumor_Fraction)), TUMOR_FRACTION_CSV
))

# ==============================================================================
# 3. Read Segment-Level CNV Data and Map to Panel Genes
# ==============================================================================
cat("Processing CNV data...\n")

# Gene coordinate lookup for mapping SV breakpoints to genes.
# PREFERRED: authoritative gene spans for panel-targeted genes, built from UCSC refGene
# intersected with the BLYMv2 design BED (see script/mutation_calling/build_panel_gene_coords.py).
# Using true gene boundaries — instead of the min/max of observed mutation positions — means
# genes with zero or a single SNV can still receive SV calls, and breakpoints near (not exactly
# on) a mutation map correctly.
GENE_COORDS_FILE <- file.path(DATA_DIR, "BLYMv2_gene_coords_hg19.tsv")

if (file.exists(GENE_COORDS_FILE)) {
  gene_coords <- read_tsv(GENE_COORDS_FILE, show_col_types = FALSE) %>%
    select(Gene, Chromosome, Gene_Start, Gene_End) %>%
    mutate(Gene_Start = as.numeric(Gene_Start), Gene_End = as.numeric(Gene_End))
  cat(sprintf("  Loaded panel gene coordinates from %s: %d genes\n",
              basename(GENE_COORDS_FILE), nrow(gene_coords)))
} else {
  # FALLBACK (legacy, less accurate): derive approximate gene boundaries from the
  # min/max position of observed Mutect2 variants per gene. Genes without SNVs are absent.
  cat("  WARNING: panel gene-coordinate file not found:\n")
  cat(sprintf("           %s\n", GENE_COORDS_FILE))
  cat("           Falling back to mutation-derived gene boundaries (less accurate for SV->gene mapping).\n")
  cat("           Generate it with: python3 script/mutation_calling/build_panel_gene_coords.py\n")
  gene_coords <- snv_data %>%
    filter(Gene != ".") %>%
    group_by(Gene, Chromosome) %>%
    summarize(
      Gene_Start = min(Position),
      Gene_End = max(Position),
      .groups = "drop"
    )
  cat(sprintf("  Built fallback gene coordinate lookup: %d genes\n", nrow(gene_coords)))
}

# Note: Gene-level CNV file (cnvkit_gene_cnvs.tsv) is not used here.
# All CNV events are mapped to their exact cytoband locations from segment data.

# Read segment-level CNV data and intersect with panel gene coordinates
cnv_seg_mapped <- tibble(Sample = character(), Gene = character(), Mutation_Class = character())
n_sensitive_cnv <- 0L

# Load cytoband data to name large regional CNVs
cytoband_file <- file.path(DATA_DIR, "cytoBand_hg19.txt")
cytobands <- NULL
if (file.exists(cytoband_file)) {
  cytobands <- read_tsv(cytoband_file, col_names = c("chrom", "start", "end", "name", "gieStain"), show_col_types = FALSE)
}

if (file.exists(CNVKIT_FILE)) {
  cnv_seg_data <- read_tsv(CNVKIT_FILE, show_col_types = FALSE)

  if (nrow(cnv_seg_data) > 0) {
    n_before <- nrow(cnv_seg_data)

    # --- CNV Artifact Filters ---
    # 1. Exclude sex chromosomes (chrX/chrY): male patients show constitutional
    #    "deletion" on chrX vs diploid reference — not somatic events
    cnv_seg_data <- cnv_seg_data %>%
      filter(!Chromosome %in% c("chrX", "chrY"))

    # 2. Exclude Immunoglobulin (Ig) Loci: V(D)J recombination in B-cell
    #    lineage causes apparent copy loss/gains — germline/lineage artifact
    cnv_seg_data <- cnv_seg_data %>%
      filter(!(
        # IGH locus (chr14:105.8M-107.3M)
        (Chromosome == "chr14" & Start <= 108000000 & End >= 105000000) |
        # IGK locus (chr2:88.8M-90.3M)
        (Chromosome == "chr2" & Start <= 91000000 & End >= 88000000) |
        # IGL locus (chr22:22M-23.5M)
        (Chromosome == "chr22" & Start <= 24000000 & End >= 21000000)
      ))

    # 2b. Exclude ACTB pseudogene mapping artifact region (chr7:5.56M-5.57M)
    #    Highly fragmented cfDNA mis-maps to pseudogenes, appearing as a focal deletion
    cnv_seg_data <- cnv_seg_data %>%
      filter(!(Chromosome == "chr7" & Start <= 5600000 & End >= 5500000))

    # 3. Require minimum probe support: segments driven by 1-4 probes are
    #    highly prone to technical noise in targeted panels
    cnv_seg_data <- cnv_seg_data %>%
      filter(Probes >= 5)

    # 3. Stricter threshold for arm-level events: large non-SubThreshold segments
    #    (>10Mb) require abs(log2) > 0.5; SubThreshold rows remain eligible for
    #    the tumor-fraction-aware sensitive tier below.
    cnv_seg_data <- cnv_seg_data %>%
      mutate(Seg_Size = End - Start) %>%
      filter(!(Call_Type != "SubThreshold" & Seg_Size > 10e6 & abs(Log2_Ratio) < 0.5)) %>%
      select(-Seg_Size)

    # 4. Bintest evidence filter: require >= 1 bin inside the segment to have
    #    passed CNVkit bintest (p_bintest <= 0.005 by default). Segments called
    #    by CBS but with zero significant bins are likely smoothing artifacts.
    #    Backward compatible: only applies if the column exists (i.e., the
    #    summarize_variant_calls.sh that emits Min_P_Value / N_Significant_Bins
    #    has been run).
    if ("N_Significant_Bins" %in% colnames(cnv_seg_data)) {
      n_pre_bintest <- nrow(cnv_seg_data)
      cnv_seg_data <- cnv_seg_data %>%
        filter(N_Significant_Bins >= 1)
      cat(sprintf(
        "  Bintest evidence filter: %d -> %d segments (removed %d with no significant bins)\n",
        n_pre_bintest, nrow(cnv_seg_data), n_pre_bintest - nrow(cnv_seg_data)
      ))
    } else {
      cat("  Bintest evidence filter: SKIPPED (N_Significant_Bins column not in TSV — re-run summarize_variant_calls.sh)\n")
    }

    cat(sprintf(
      "  CNV segments: %d total -> %d after artifact filters (removed %d)\n",
      n_before, nrow(cnv_seg_data), n_before - nrow(cnv_seg_data)
    ))

    if (nrow(cnv_seg_data) > 0) {
      cnv_seg_data <- cnv_seg_data %>%
        left_join(tumor_fraction %>% select(Sample, Tumor_Fraction), by = "Sample")

      # High-confidence tier: preserve the prior behavior exactly by admitting only
      # strict CNVkit Amplification/Deletion rows and mapping CN=0 deletions separately.
      cnv_high_conf <- cnv_seg_data %>%
        filter(Call_Type %in% c("Amplification", "Deletion")) %>%
        mutate(Mutation_Class = case_when(
          Call_Type == "Amplification" ~ "Amplification",
          Copy_Number == 0 ~ "Deep_Deletion",
          TRUE ~ "Deletion"
        ))

      # Sensitive low-tumor-fraction tier: evaluate only SubThreshold segments from
      # samples with an SNV-derived tumor fraction >=5%. The adaptive threshold is
      # anchored to the expected single-copy gain/loss log2 shift at purity p, with a
      # hard 0.1 log2 noise floor so very small estimated shifts do not create calls.
      NOISE_FLOOR <- 0.1
      cnv_sensitive <- cnv_seg_data %>%
        filter(Call_Type == "SubThreshold", !is.na(Tumor_Fraction), Tumor_Fraction >= 0.05) %>%
        mutate(
          gain_shift = log2((2 + Tumor_Fraction) / 2),
          loss_shift = abs(log2((2 - Tumor_Fraction) / 2)),
          sens_gain = pmax(NOISE_FLOOR, 0.5 * gain_shift),
          sens_loss = pmax(NOISE_FLOOR, 0.5 * loss_shift),
          Mutation_Class = case_when(
            Log2_Ratio >= sens_gain ~ "Amplification_LowTF",
            Log2_Ratio <= -sens_loss ~ "Deletion_LowTF",
            TRUE ~ NA_character_
          )
        ) %>%
        filter(!is.na(Mutation_Class)) %>%
        select(-gain_shift, -loss_shift, -sens_gain, -sens_loss)

      n_sensitive_cnv <- nrow(cnv_sensitive)
      cnv_seg_to_map <- bind_rows(cnv_high_conf, cnv_sensitive)

      if (nrow(cnv_seg_to_map) > 0) {
        # We map EVERY tiered CNV segment (focal or large) to its exact cytoband.
        cnv_hits <- list()

        for (j in seq_len(nrow(cnv_seg_to_map))) {
          seg_j <- cnv_seg_to_map[j, ]

          feature_name <- paste0(seg_j$Chromosome, " CNV") # Fallback

          if (!is.null(cytobands)) {
            midpoint <- seg_j$Start + (seg_j$End - seg_j$Start) / 2
            band_hit <- cytobands %>%
              filter(chrom == seg_j$Chromosome, start <= midpoint, end >= midpoint)

            if (nrow(band_hit) > 0) {
              clean_chr <- str_remove(seg_j$Chromosome, "^chr")
              # Exact cytoband (e.g., "6q25.1", "2p16.1")
              feature_name <- paste0(clean_chr, band_hit$name[1])
            }
          }

          cnv_hits[[length(cnv_hits) + 1]] <- tibble(
            Sample = seg_j$Sample,
            Gene = feature_name, # Storing location replacing the 'Gene' row name
            Log2_Ratio = seg_j$Log2_Ratio,
            Copy_Number = seg_j$Copy_Number,
            Call_Type = seg_j$Call_Type,
            Mutation_Class = seg_j$Mutation_Class
          )
        }

        cnv_seg_mapped <- bind_rows(cnv_hits) %>%
          select(Sample, Gene, Mutation_Class) %>%
          distinct()
      }

      cat(sprintf("  Processed %d CNVs exactly to their locations.\n", nrow(cnv_seg_mapped)))
    }
  }
}

cat(sprintf("  Sensitive low tumor-fraction CNVs added: %d\n", n_sensitive_cnv))

# The user explicitly requested to show ALL deletions and amplifications at their exact locations.
# Therefore, we discard the gene-level definitions for CNVs and exclusively use the segment locations.
cnv_qc <- cnv_seg_mapped

cat(sprintf(
  "  Total CNV events for OncoPrint: %d across %d genes\n",
  nrow(cnv_qc), n_distinct(cnv_qc$Gene)
))

# ==============================================================================
# 4. Read Structural Variant Data (Manta)
# ==============================================================================
cat("Processing structural variant data...\n")

sv_qc <- tibble(Sample = character(), Gene = character(), Mutation_Class = character())

if (file.exists(SV_FILE)) {
  sv_data <- read_tsv(SV_FILE, show_col_types = FALSE)

  # Keep only PASS structural variants. Manta flags low-confidence calls
  # (e.g. MinSomaticScore); without this they would leak into the OncoPrint.
  if ("Filter" %in% colnames(sv_data)) {
    n_sv_before <- nrow(sv_data)
    sv_data <- sv_data %>% filter(Filter == "PASS")
    cat(sprintf("  SV PASS filter: %d -> %d structural variants\n", n_sv_before, nrow(sv_data)))
  }

  if (nrow(sv_data) > 0) {
    # Helper: panel genes overlapping a breakpoint window chr:[min(p1,p2)-flank, max(p1,p2)+flank]
    genes_near <- function(chrom, p1, p2, flank = 50000) {
      if (is.na(p1)) return(character(0))
      if (is.na(p2)) p2 <- p1
      gene_coords %>%
        filter(
          Chromosome == chrom,
          Gene_Start <= max(p1, p2) + flank,
          Gene_End   >= min(p1, p2) - flank
        ) %>%
        pull(Gene)
    }

    # Map SV breakpoints to nearby panel genes
    sv_gene_hits <- list()

    for (i in seq_len(nrow(sv_data))) {
      sv <- sv_data[i, ]
      sv_pos   <- as.numeric(sv$Position)
      mate_pos <- as.numeric(sv$Mate_Pos)
      mate_chr <- if ("Mate_Chrom" %in% colnames(sv_data)) sv$Mate_Chrom else NA_character_
      sv_type  <- as.character(sv$SV_Type)
      # Effective mate chromosome for focal mapping: fall back to the event chromosome when
      # the mate chromosome is unknown ("." / NA) — e.g. Manta INV/BND records that encode the
      # partner via END rather than CHR2 — so the mate breakpoint still maps to real genes.
      mate_chr_eff <- if (is.na(mate_chr) || mate_chr == ".") sv$Chromosome else mate_chr

      if ((!is.na(mate_chr) && mate_chr != "." && mate_chr != sv$Chromosome) ||
          sv_type %in% c("BND", "INV")) {
        # Breakpoint events (including same-chromosome BND/INV) map EACH breakpoint
        # focally; do NOT span the interval.
        hit_genes <- unique(c(
          genes_near(sv$Chromosome, sv_pos,  sv_pos),
          genes_near(mate_chr_eff,  mate_pos, mate_pos)
        ))
      } else if (sv_type %in% c("DEL", "DUP")) {
        # Intrachromosomal DEL/DUP: map genes across the spanned region.
        hit_genes <- genes_near(sv$Chromosome, sv_pos, mate_pos)
      } else {
        hit_genes <- unique(c(
          genes_near(sv$Chromosome, sv_pos,  sv_pos),
          genes_near(mate_chr_eff,  mate_pos, mate_pos)
        ))
      }

      if (length(hit_genes) > 0) {
        sv_gene_hits[[length(sv_gene_hits) + 1]] <- tibble(
          Sample = sv$Sample,
          Gene = hit_genes,
          SV_Type = sv$SV_Type,
          Mutation_Class = paste0("SV_", sv$SV_Type)
        )
      }
    }

    if (length(sv_gene_hits) > 0) {
      sv_qc <- bind_rows(sv_gene_hits) %>%
        filter(!str_detect(Gene, "^LINC|^LOC|^MIR\\d|^OR\\d|-AS\\d|ACTB|IGLL5|ABO|MPEG1|ELK2AP")) %>%
        select(Sample, Gene, Mutation_Class) %>%
        distinct()

      cat(sprintf("  Structural variants mapped to %d gene-sample pairs\n", nrow(sv_qc)))
    } else {
      cat("  No SVs mapped to panel genes.\n")
    }
  } else {
    cat("  No somatic SVs detected.\n")
  }
}

# ==============================================================================
# 5. Build Combined OncoPrint Matrix
# ==============================================================================
cat("Building OncoPrint matrix...\n")

# Collapse SNV/Indel mutations per gene/sample
snv_onco <- snv_qc %>%
  group_by(Sample, Gene) %>%
  summarize(Mutations = paste(unique(Mutation_Class), collapse = ";"), .groups = "drop")

# Collapse CNV calls per gene/sample
cnv_onco <- cnv_qc %>%
  group_by(Sample, Gene) %>%
  summarize(Mutations = paste(unique(Mutation_Class), collapse = ";"), .groups = "drop")

# Collapse SV calls per gene/sample
sv_onco <- sv_qc %>%
  group_by(Sample, Gene) %>%
  summarize(Mutations = paste(unique(Mutation_Class), collapse = ";"), .groups = "drop")

# Merge all alteration types
onco_df <- bind_rows(snv_onco, cnv_onco, sv_onco) %>%
  group_by(Sample, Gene) %>%
  summarize(Mutations = paste(unique(Mutations), collapse = ";"), .groups = "drop")

# Clean patient IDs: strip "-T1" suffix for cleaner column names
onco_df <- onco_df %>%
  mutate(Sample = str_remove(Sample, "-T1$"))

# Reshape to wide format matrix (Rows = Genes, Cols = Samples)
onco_mat <- onco_df %>%
  pivot_wider(names_from = Sample, values_from = Mutations, values_fill = "") %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Separate mutation rows (gene names) from CNV rows (cytoband locations)
all_cnv_features <- unique(cnv_qc$Gene)
is_cnv <- rownames(onco_mat) %in% all_cnv_features

# For mutations: keep top 40 most frequently altered genes
mut_rows <- onco_mat[!is_cnv, , drop = FALSE]
mut_counts <- rowSums(mut_rows != "")
# Use seq_len() to avoid the R `1:0` foot-gun: when there are no mutation rows,
# `1:min(40, 0)` would yield c(1, 0) and mis-index; seq_len(0) is integer(0).
top_mut <- names(sort(mut_counts, decreasing = TRUE)[seq_len(min(40, length(mut_counts)))])

# For CNVs: keep ALL cytoband rows (don't cut any — they're already QC-filtered)
cnv_rows <- onco_mat[is_cnv, , drop = FALSE]
cnv_rows_nonempty <- cnv_rows[rowSums(cnv_rows != "") > 0, , drop = FALSE]

# Combine: mutations on top, CNVs on bottom
onco_mat_plot <- rbind(
  onco_mat[top_mut, , drop = FALSE],
  cnv_rows_nonempty
)

# Add cohort samples with no post-QC alterations as empty OncoPrint columns.
cohort_samples <- sort(unique(str_remove(all_cohort_samples, "-T1$")))
missing_cohort_samples <- setdiff(cohort_samples, colnames(onco_mat_plot))
if (length(missing_cohort_samples) > 0) {
  plot_rownames <- rownames(onco_mat_plot)
  if (is.null(plot_rownames)) {
    plot_rownames <- character(0)
  }
  empty_cols <- matrix(
    "",
    nrow = nrow(onco_mat_plot),
    ncol = length(missing_cohort_samples),
    dimnames = list(plot_rownames, missing_cohort_samples)
  )
  onco_mat_plot <- cbind(onco_mat_plot, empty_cols)
}
cat(sprintf(
  "  Added %d alteration-free cohort samples to OncoPrint: %s\n",
  length(missing_cohort_samples),
  ifelse(length(missing_cohort_samples) > 0, paste(missing_cohort_samples, collapse = ", "), "none")
))

# Build row_split vector matching the combined matrix
row_splits <- c(
  rep("Mutations (SNV/Indel/SV)", length(top_mut)),
  rep("Copy Number Alterations", nrow(cnv_rows_nonempty))
)
# Ensure the order is "Mutations" then "Copy Number Alterations" in the plot
row_splits <- factor(row_splits, levels = c("Mutations (SNV/Indel/SV)", "Copy Number Alterations"))

cat(sprintf("  Final OncoPrint: %d rows x %d samples\n", nrow(onco_mat_plot), ncol(onco_mat_plot)))

# Print alteration type breakdown
all_types <- unlist(str_split(onco_mat_plot[onco_mat_plot != ""], ";"))
cat("  Alteration type counts in OncoPrint:\n")
print(sort(table(all_types), decreasing = TRUE))

# ==============================================================================
# 6. Define Colors and Draw OncoPrint
# ==============================================================================

# Standard cBioPortal-style color palette for OncoPrint
col <- c(
  # SNV types (cBioPortal convention: green for missense, black for truncating)
  "Missense_SNV"     = "#008000", # Green (cBioPortal standard)
  "Nonsense_SNV"     = "#000000", # Black (truncating)
  "Splice_Site"      = "#9B59B6", # Purple
  "UTR_SNV"          = "#F39C12", # Gold/Amber
  "Other_SNV"        = "#7F8C8D", # Gray
  "aSHM_Variant"     = "#85C1AE", # Teal — aSHM-associated, likely passenger (e.g. IGLL5)
  "Rescued_Driver_LowVAF" = "#4DBD74", # Light green with hatch marks in alter_fun
  # Indel types
  "Frameshift_Indel" = "#2C3E50", # Dark Navy Blue (differentiated from Nonsense)
  "Inframe_Indel"    = "#993404", # Brown (cBioPortal standard)
  # CNV types (cBioPortal convention: red for amp, blue for del)
  "Amplification"    = "#FF0000", # Red (cBioPortal standard)
  "Amplification_LowTF" = "#E79A9A", # Desaturated red for tumor-fraction-aware sensitive tier
  "Deletion"         = "#56B4E9", # Light Blue (heterozygous loss)
  "Deletion_LowTF"   = "#9CC7E6", # Desaturated blue for tumor-fraction-aware sensitive tier
  "Deep_Deletion"    = "#0000FF", # Blue (cBioPortal standard for homozygous del)
  # Structural Variant types
  "SV_DEL"           = "#FDB462", # Light Orange
  "SV_DUP"           = "#FB8072", # Salmon
  "SV_INV"           = "#BEBADA", # Lavender
  "SV_BND"           = "#BC80BD" # Mauve
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = "#ECECEC", col = NA))
  },
  # --- SNV types: full height, full width ---
  "Missense_SNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Missense_SNV"], col = NA))
  },
  "Nonsense_SNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Nonsense_SNV"], col = NA))
  },
  "Splice_Site" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  "UTR_SNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["UTR_SNV"], col = NA))
  },
  "Other_SNV" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.9, gp = gpar(fill = col["Other_SNV"], col = NA))
  },
  # aSHM-associated (likely passenger): half-height to visually denote lower confidence
  "aSHM_Variant" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.5, gp = gpar(fill = col["aSHM_Variant"], col = NA))
  },
  # Near-threshold driver rescue: shorter light-green rectangle with cross-hatch overlay.
  "Rescued_Driver_LowVAF" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.55, gp = gpar(fill = col["Rescued_Driver_LowVAF"], col = NA))
    grid.segments(x - w * 0.35, y - h * 0.20, x + w * 0.35, y + h * 0.20,
                  gp = gpar(col = "#1B5E20", lwd = 0.8))
    grid.segments(x - w * 0.35, y + h * 0.20, x + w * 0.35, y - h * 0.20,
                  gp = gpar(col = "#1B5E20", lwd = 0.8))
  },
  # --- Indel types: slightly shorter (cBioPortal uses ~2/3 height for inframe) ---
  "Frameshift_Indel" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.66, gp = gpar(fill = col["Frameshift_Indel"], col = NA))
  },
  "Inframe_Indel" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.66, gp = gpar(fill = col["Inframe_Indel"], col = NA))
  },
  # --- CNV types: full height in their own section ---
  "Amplification" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Amplification"], col = NA))
  },
  "Amplification_LowTF" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.5, gp = gpar(fill = col["Amplification_LowTF"], col = NA))
  },
  "Deletion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Deletion"], col = NA))
  },
  "Deletion_LowTF" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.5, gp = gpar(fill = col["Deletion_LowTF"], col = NA))
  },
  "Deep_Deletion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Deep_Deletion"], col = NA))
  },
  # --- SV types: triangle marker ---
  "SV_DEL" = function(x, y, w, h) {
    grid.points(x, y, pch = 24, size = unit(0.6, "char"), gp = gpar(fill = col["SV_DEL"], col = "black"))
  },
  "SV_DUP" = function(x, y, w, h) {
    grid.points(x, y, pch = 24, size = unit(0.6, "char"), gp = gpar(fill = col["SV_DUP"], col = "black"))
  },
  "SV_INV" = function(x, y, w, h) {
    grid.points(x, y, pch = 24, size = unit(0.6, "char"), gp = gpar(fill = col["SV_INV"], col = "black"))
  },
  "SV_BND" = function(x, y, w, h) {
    grid.points(x, y, pch = 24, size = unit(0.6, "char"), gp = gpar(fill = col["SV_BND"], col = "black"))
  }
)

# Only include colors/legend entries for alteration types actually present in the data
types_present <- unique(unlist(str_split(onco_mat_plot[onco_mat_plot != ""], ";")))
col_used <- col[names(col) %in% types_present]
alter_fun_used <- alter_fun[c("background", names(col_used))]

# Create human-readable legend labels
legend_labels <- names(col_used) %>%
  str_replace("^aSHM_Variant$", "aSHM-associated (likely passenger)") %>%
  str_replace("^Rescued_Driver_LowVAF$", "Driver rescue (0.25-0.5% VAF, low-confidence)") %>%
  str_replace("^Amplification_LowTF$", "Amplification (low tumor-fraction)") %>%
  str_replace("^Deletion_LowTF$", "Deletion (low tumor-fraction)") %>%
  str_replace("_SNV", " (SNV)") %>%
  str_replace("_Indel", " (Indel)") %>%
  str_replace("Deep_Deletion", "Deep Deletion (CN=0)") %>%
  str_replace("^Deletion$", "Deletion (CN=1)") %>%
  str_replace("SV_DEL", "SV Deletion") %>%
  str_replace("SV_DUP", "SV Duplication") %>%
  str_replace("SV_INV", "SV Inversion") %>%
  str_replace("SV_BND", "SV Breakend")

heatmap_title <- "CAPP-Seq Somatic Alterations"
heatmap_subtitle <- "SNVs | Indels | Copy Number Variants | Structural Variants"

# Generate PDF
pdf(OUTPUT_PDF, width = max(14, ncol(onco_mat_plot) * 0.6 + 4), height = max(10, nrow(onco_mat_plot) * 0.35 + 3))

if (nrow(onco_mat_plot) == 0 || ncol(onco_mat_plot) == 0) {
  cat("  No variants left to plot after filtering. Creating an empty PDF placeholder.\n")
  plot.new()
  text(0.5, 0.5, "No somatic variants passed filtering across all samples.", cex = 1.2)
} else {
  p <- oncoPrint(
    onco_mat_plot,
    alter_fun = alter_fun_used,
    col = col_used,
    remove_empty_columns = FALSE,
    remove_empty_rows = TRUE,
    row_split = row_splits,
    cluster_row_slices = FALSE, # Prevents flipping the order of the two sections randomly
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = gpar(fontsize = 9, fontface = "italic"),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_names_side = "bottom",
    pct_gp = gpar(fontsize = 8),
    column_title = paste0(heatmap_title, "\n", heatmap_subtitle),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Alterations",
      at = names(col_used),
      labels = legend_labels,
      nrow = 3
    ),
    top_annotation = HeatmapAnnotation(
      cbar = anno_oncoprint_barplot(
        border = TRUE
      ),
      annotation_height = unit(2, "cm")
    ),
    right_annotation = rowAnnotation(
      rbar = anno_oncoprint_barplot(
        border = TRUE
      ),
      annotation_width = unit(2, "cm")
    )
  )
  draw(p)
}

dev.off()

cat(sprintf("\nOncoPrint generated successfully: %s\n", OUTPUT_PDF))

# ==============================================================================
# 7. Sequencing Depth QC Plots
# ==============================================================================
cat("Generating sequencing depth QC plots...\n")

QC_PDF <- file.path(PROJECT_DIR, "analysis", "sequencing_depth_qc.pdf")

# Parse numeric depths from raw data (all variants, not just PASS)
depth_data <- snv_data %>%
  mutate(
    Tumor_DP_num  = as.numeric(Tumor_DP),
    Normal_DP_num = as.numeric(Normal_DP)
  ) %>%
  filter(!is.na(Tumor_DP_num), !is.na(Normal_DP_num))

# Strip "-T1" suffix to get clean patient IDs for labeling
depth_data <- depth_data %>%
  mutate(Patient = str_remove(Sample, "-T1$"))

# --- Callable-depth / effective-sensitivity report ---
# The genotyping floor is governed by BOTH Tumor_AF >= 0.5% AND Tumor_AD_Alt >= 4, so the
# minimum DETECTABLE VAF at a site is ~ 4 / Tumor_DP. Thus a 0.5% variant needs DP ~ 800 and
# a 1% variant needs DP ~ 400. This table reports, per sample, the median depth, the implied
# minimum detectable VAF at that median, and the fraction of sites reaching the 1% / 0.5%
# depths. CAVEAT: computed over CALLED candidate sites in mutect2_variants.tsv — a proxy for
# panel-wide coverage. True callable sensitivity needs per-base depth across the panel.
callable_report <- depth_data %>%
  group_by(Patient, Batch) %>%
  summarize(
    N_sites = n(),
    Median_Tumor_DP = median(Tumor_DP_num),
    Median_Normal_DP = median(Normal_DP_num),
    Implied_min_VAF_at_median = round(4 / median(Tumor_DP_num), 4),
    Frac_sites_DP_ge_400 = round(mean(Tumor_DP_num >= 400), 3), # ~1.0% LOD achievable
    Frac_sites_DP_ge_800 = round(mean(Tumor_DP_num >= 800), 3), # ~0.5% LOD achievable
    .groups = "drop"
  ) %>%
  arrange(Median_Tumor_DP)

CALLABLE_CSV <- file.path(PROJECT_DIR, "analysis", "callable_depth_sensitivity.csv")
write_csv(callable_report, CALLABLE_CSV)
cat(sprintf("  Callable-depth/sensitivity report -> %s\n", CALLABLE_CSV))
cat("  (Proxy over called sites; min detectable VAF ~= 4 / Tumor_DP given the AD_Alt>=4 floor.)\n")
print(callable_report, n = Inf)

pdf(QC_PDF, width = 14, height = 10)

# --- Plot 1: Per-sample Tumor Depth (Boxplot) ---
# Calculate per-sample median for labeling
sample_medians <- depth_data %>%
  group_by(Patient, Batch) %>%
  summarize(Median_Tumor_DP = median(Tumor_DP_num), .groups = "drop")

p1 <- ggplot(depth_data, aes(x = reorder(Patient, Tumor_DP_num, FUN = median), y = Tumor_DP_num, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = 110, label = "Depth ≥ 100 threshold", hjust = 0, color = "red", size = 3) +
  coord_cartesian(ylim = c(0, quantile(depth_data$Tumor_DP_num, 0.99))) +
  labs(
    title = "Tumor (cfDNA) Sequencing Depth per Sample",
    subtitle = "Dashed red line = minimum depth threshold (100x) from CAPP-Seq method",
    x = "Patient", y = "Tumor Depth (deduplicated reads)",
    fill = "Batch"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# --- Plot 2: Per-sample Normal Depth (Boxplot) ---
p2 <- ggplot(depth_data, aes(x = reorder(Patient, Normal_DP_num, FUN = median), y = Normal_DP_num, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = 25, label = "Depth ≥ 20 threshold", hjust = 0, color = "red", size = 3) +
  coord_cartesian(ylim = c(0, quantile(depth_data$Normal_DP_num, 0.99))) +
  labs(
    title = "Normal (PBL) Sequencing Depth per Sample",
    subtitle = "Dashed red line = minimum depth threshold (20x) from CAPP-Seq method",
    x = "Patient", y = "Normal Depth (deduplicated reads)",
    fill = "Batch"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)

# --- Plot 3: Depth Distribution Density (Tumor vs Normal overlay) ---
depth_long <- depth_data %>%
  select(Patient, Tumor_DP_num, Normal_DP_num) %>%
  pivot_longer(
    cols = c(Tumor_DP_num, Normal_DP_num),
    names_to = "Type", values_to = "Depth"
  ) %>%
  mutate(Type = ifelse(Type == "Tumor_DP_num", "Tumor (cfDNA)", "Normal (PBL)"))

p3 <- ggplot(depth_long, aes(x = Depth, fill = Type)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "#E31A1C", linewidth = 0.6) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "#1F78B4", linewidth = 0.6) +
  scale_x_continuous(limits = c(0, quantile(depth_long$Depth, 0.99))) +
  labs(
    title = "Sequencing Depth Distribution: Tumor vs Normal",
    subtitle = "Red dashed = Tumor threshold (100x) | Blue dashed = Normal threshold (20x)",
    x = "Sequencing Depth", y = "Density",
    fill = "Sample Type"
  ) +
  theme_minimal(base_size = 12)
print(p3)

# --- Plot 4: Per-Sample Median Summary Table (as a bar chart) ---
summary_stats <- depth_data %>%
  group_by(Patient, Batch) %>%
  summarize(
    Median_Tumor_DP = median(Tumor_DP_num),
    Median_Normal_DP = median(Normal_DP_num),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Median_Tumor_DP, Median_Normal_DP),
    names_to = "Type", values_to = "Median_Depth"
  ) %>%
  mutate(Type = ifelse(Type == "Median_Tumor_DP", "Tumor (cfDNA)", "Normal (PBL)"))

p4 <- ggplot(summary_stats, aes(x = reorder(Patient, Median_Depth), y = Median_Depth, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(
    title = "Median Sequencing Depth per Patient",
    x = "Patient", y = "Median Depth",
    fill = "Sample Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p4)

dev.off()

cat(sprintf("Depth QC plots generated successfully: %s\n", QC_PDF))
