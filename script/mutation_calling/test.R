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
OUTPUT_PDF <- file.path(PROJECT_DIR, "analysis", "cohort_oncoprint.pdf")

# PoN frequency summary from build_mutect2_pon.sh (set to NA to skip PoN filtering)
# Lives in combined_summary/ as written by summarize_variant_calls.sh on Sherlock.
PON_SUMMARY_FILE <- file.path(RESULTS_DIR, "pon_site_summary.tsv")

# COSMIC coding mutations VCF for rescue logic (set to NA to skip COSMIC rescue)
# Download from: https://cancer.sanger.ac.uk/cosmic/download (requires free registration)
# Example: CosmicCodingMuts.vcf.gz (hg19/GRCh37)
COSMIC_VCF <- file.path(DATA_DIR, "Cosmic_CompleteTargetedScreensMutant_v103_GRCh37.vcf.gz")

cat("Loading and processing variant data...\n")

# ==============================================================================
# 2. Read and QC SNV/Indel Data (Mutect2 + SnpEff)
# ==============================================================================
# Read data
snv_data <- read_tsv(MUTECT_FILE, show_col_types = FALSE)

# QC Filter 1: Only keep high-confidence PASS variants
snv_filt <- snv_data %>%
  filter(Filter == "PASS")

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

  # Read COSMIC VCF (just chr, pos, ref, alt, ID columns)
  cosmic_lines <- readLines(gzfile(COSMIC_VCF))
  cosmic_lines <- cosmic_lines[!startsWith(cosmic_lines, "#")]

  if (length(cosmic_lines) > 0) {
    cosmic_df <- read_tsv(
      I(cosmic_lines), col_names = FALSE, show_col_types = FALSE,
      col_select = c(1, 2, 3, 4, 5)
    ) %>%
      rename(Chromosome = X1, Position = X2, COSMIC_ID = X3, Ref = X4, Alt = X5) %>%
      mutate(Position = as.numeric(Position)) %>%
      distinct(Chromosome, Position, Ref, Alt, .keep_all = TRUE)

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
  )

cat(sprintf(
  "  Mutect2 after QC: %d protein-altering variants across %d genes\n",
  nrow(snv_qc), n_distinct(snv_qc$Gene)
))
cat(sprintf(
  "    SNVs: %d | Indels: %d\n",
  sum(snv_qc$Mutation_Class %in% c("Missense_SNV", "Nonsense_SNV", "Splice_Site", "UTR_SNV", "Other_SNV")),
  sum(snv_qc$Mutation_Class %in% c("Frameshift_Indel", "Inframe_Indel"))
))

# ==============================================================================
# 3. Read Segment-Level CNV Data and Map to Panel Genes
# ==============================================================================
cat("Processing CNV data...\n")

# Build a gene coordinate lookup from the Mutect2 data (all variants, not just PASS)
# Extract min/max position per gene per chromosome to define gene boundaries
gene_coords <- snv_data %>%
  filter(Gene != ".") %>%
  group_by(Gene, Chromosome) %>%
  summarize(
    Gene_Start = min(Position),
    Gene_End = max(Position),
    .groups = "drop"
  )

cat(sprintf("  Built gene coordinate lookup: %d genes\n", nrow(gene_coords)))

# Note: Gene-level CNV file (cnvkit_gene_cnvs.tsv) is not used here.
# All CNV events are mapped to their exact cytoband locations from segment data.

# Read segment-level CNV data and intersect with panel gene coordinates
cnv_seg_mapped <- tibble(Sample = character(), Gene = character(), Mutation_Class = character())

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

    # 2. Exclude IGH locus region (chr14:35.8M-69.3M): V(D)J recombination in
    #    B-cell lineage causes apparent copy loss — germline/lineage artifact
    cnv_seg_data <- cnv_seg_data %>%
      filter(!(Chromosome == "chr14" & Start >= 35000000 & End <= 70000000))

    # 3. Require minimum probe support: segments driven by 1-4 probes are
    #    highly prone to technical noise in targeted panels
    cnv_seg_data <- cnv_seg_data %>%
      filter(Probes >= 5)

    # 3. Stricter threshold for arm-level events: large segments (>10Mb) require
    #    abs(log2) > 0.5 to reduce noise from GC bias and coverage unevenness
    cnv_seg_data <- cnv_seg_data %>%
      mutate(Seg_Size = End - Start) %>%
      filter(!(Seg_Size > 10e6 & abs(Log2_Ratio) < 0.5)) %>%
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
      # We map EVERY CNV segment (focal or large) to its exact cytoband
      cnv_hits <- list()

      for (j in seq_len(nrow(cnv_seg_data))) {
        seg_j <- cnv_seg_data[j, ]

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
          Call_Type = seg_j$Call_Type
        )
      }

      cnv_seg_mapped <- bind_rows(cnv_hits) %>%
        mutate(Mutation_Class = case_when(
          Call_Type == "Amplification" ~ "Amplification",
          Copy_Number == 0 ~ "Deep_Deletion",
          TRUE ~ "Deletion"
        )) %>%
        select(Sample, Gene, Mutation_Class) %>%
        distinct()

      cat(sprintf("  Processed %d CNVs exactly to their locations.\n", nrow(cnv_seg_mapped)))
    }
  }
}

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

  if (nrow(sv_data) > 0) {
    # Map SV breakpoints to nearby panel genes
    sv_gene_hits <- list()

    for (i in seq_len(nrow(sv_data))) {
      sv <- sv_data[i, ]
      sv_pos <- as.numeric(sv$Position)
      sv_end <- as.numeric(sv$Mate_Pos)
      if (is.na(sv_end)) sv_end <- sv_pos

      # Find genes on the same chromosome near the breakpoint (within the SV span)
      hits <- gene_coords %>%
        filter(
          Chromosome == sv$Chromosome,
          Gene_Start <= max(sv_pos, sv_end) + 50000, # 50kb flanking
          Gene_End >= min(sv_pos, sv_end) - 50000
        )

      if (nrow(hits) > 0) {
        sv_gene_hits[[length(sv_gene_hits) + 1]] <- tibble(
          Sample = sv$Sample,
          Gene = hits$Gene,
          SV_Type = sv$SV_Type,
          Mutation_Class = paste0("SV_", sv$SV_Type)
        )
      }
    }

    if (length(sv_gene_hits) > 0) {
      sv_qc <- bind_rows(sv_gene_hits) %>%
        filter(!str_detect(Gene, "^LINC|^LOC|^MIR\\d|^OR\\d|-AS\\d")) %>%
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
top_mut <- names(sort(mut_counts, decreasing = TRUE)[1:min(40, length(mut_counts))])

# For CNVs: keep ALL cytoband rows (don't cut any — they're already QC-filtered)
cnv_rows <- onco_mat[is_cnv, , drop = FALSE]
cnv_rows_nonempty <- cnv_rows[rowSums(cnv_rows != "") > 0, , drop = FALSE]

# Combine: mutations on top, CNVs on bottom
onco_mat_plot <- rbind(
  onco_mat[top_mut, , drop = FALSE],
  cnv_rows_nonempty
)

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
  # Indel types
  "Frameshift_Indel" = "#000000", # Black (truncating, same as nonsense per cBioPortal)
  "Inframe_Indel"    = "#993404", # Brown (cBioPortal standard)
  # CNV types (cBioPortal convention: red for amp, blue for del)
  "Amplification"    = "#FF0000", # Red (cBioPortal standard)
  "Deletion"         = "#56B4E9", # Light Blue (heterozygous loss)
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
  # --- Indel types: slightly shorter (cBioPortal uses ~2/3 height for inframe) ---
  "Frameshift_Indel" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Frameshift_Indel"], col = NA))
  },
  "Inframe_Indel" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h * 0.66, gp = gpar(fill = col["Inframe_Indel"], col = NA))
  },
  # --- CNV types: full height in their own section ---
  "Amplification" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Amplification"], col = NA))
  },
  "Deletion" = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "pt"), h - unit(0.5, "pt"), gp = gpar(fill = col["Deletion"], col = NA))
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
