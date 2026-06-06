#!/usr/bin/env Rscript
# =============================================================================
# 06_visualize.R  --  QC + cohort landscape + Kraken2-vs-Kaiju concordance
# Reads step-05 outputs and writes PDFs to FIGS_DIR.
#   Args: 1 RESULTS_DIR   2 FIGS_DIR
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(scales)
})
op <- function(p) tryCatch(p, error = function(e) { message("plot skipped: ", conditionMessage(e)); NULL })

a <- commandArgs(trailingOnly = TRUE)
RESULTS_DIR <- a[1]; FIGS_DIR <- ifelse(length(a) >= 2, a[2], file.path(RESULTS_DIR, "figures"))
dir.create(FIGS_DIR, showWarnings = FALSE, recursive = TRUE)

xv_path <- file.path(RESULTS_DIR, "cross_validated_species.csv")
if (!file.exists(xv_path)) stop("Missing cross_validated_species.csv; run step 05 first.")
xv <- read_csv(xv_path, show_col_types = FALSE)

# ---- 1. QC: read accounting + classification rate ---------------------------
qc_path <- file.path(RESULTS_DIR, "QC_read_accounting.csv")
if (file.exists(qc_path)) {
  qc <- read_csv(qc_path, show_col_types = FALSE)
  pdf(file.path(FIGS_DIR, "QC_read_accounting.pdf"), width = 12, height = 8)
  op(print(
    qc %>% select(sample, star_unmapped, final_microbial) %>%
      pivot_longer(-sample) %>%
      ggplot(aes(reorder(sample, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_y_continuous(labels = comma) + theme_minimal(base_size = 8) +
      labs(title = "Read accounting per sample",
           subtitle = "STAR-unmapped vs final candidate-microbial pairs",
           x = "", y = "read pairs", fill = "")))
  op(print(
    qc %>% ggplot(aes(pct_microbial_final)) +
      geom_histogram(bins = 40, fill = "#4682B4") + theme_minimal() +
      labs(title = "Final microbial fraction of input reads",
           subtitle = "Very high values (>>1-5%) often indicate residual host / artifact",
           x = "% of input read pairs", y = "samples")))
  # reads removed at each cleanup stage (step 01)
  rmcols <- c("host_removed_STAR","rRNA_removed","lowcomplexity_removed","duplicates_removed","t2t_removed")
  if (all(rmcols %in% names(qc))) op(print(
    qc %>% select(sample, all_of(rmcols)) %>%
      pivot_longer(-sample, names_to = "stage", values_to = "reads") %>%
      mutate(stage = factor(stage, levels = rmcols)) %>%
      ggplot(aes(stage, reads)) + geom_boxplot(fill = "#4682B4", alpha = .6) +
      scale_y_continuous(labels = comma) + theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
      labs(title = "Read pairs removed per cleanup stage (across samples)",
           x = "", y = "read pairs removed")))
  # Kraken2 unclassified rate (reads not passing confidence / min-hit-groups, or no DB hit)
  if ("pct_unclassified" %in% names(qc)) op(print(
    qc %>% ggplot(aes(as.numeric(pct_unclassified))) +
      geom_histogram(bins = 40, fill = "#6A5ACD") + theme_minimal() +
      labs(title = "Kraken2 unclassified rate per sample",
           subtitle = "fraction of read pairs left unclassified after confidence/min-hit-groups",
           x = "% unclassified read pairs", y = "samples")))
  dev.off()
}

# ---- 1b. Filter funnel ------------------------------------------------------
ff_path <- file.path(RESULTS_DIR, "filter_funnel_summary.csv")
if (file.exists(ff_path)) {
  ff <- read_csv(ff_path, show_col_types = FALSE)
  pdf(file.path(FIGS_DIR, "filter_funnel.pdf"), width = 11, height = 6)
  op(print(
    ff %>% mutate(stage = factor(stage, levels = rev(stage))) %>%
      ggplot(aes(stage, calls)) + geom_col(fill = "#B22222", alpha = .8) +
      coord_flip() + scale_y_continuous(labels = comma) + theme_minimal(base_size = 11) +
      labs(title = "Filter funnel: sample-taxon calls surviving each gate",
           subtitle = "how many Kraken2/Bracken calls remain after each filter",
           x = "", y = "calls")))
  dev.off()
}

# ---- 2. Cohort high-confidence landscape ------------------------------------
hc <- xv %>% filter(high_confidence)
landscape <- function(df, kingdom) {
  d <- df %>% filter(superkingdom == kingdom)
  if (nrow(d) == 0) return(ggplot() + theme_void() + labs(title = paste("No", kingdom)))
  top <- d %>% group_by(species) %>%
    summarise(n = n_distinct(sample), tot = sum(RPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n), desc(tot)) %>% slice_head(n = 15)
  d %>% filter(species %in% top$species) %>%
    mutate(species = factor(species, levels = rev(top$species))) %>%
    ggplot(aes(species, RPM, fill = sample)) + geom_col() + coord_flip() +
    theme_minimal(base_size = 10) + theme(legend.position = "none") +
    labs(title = paste(kingdom, "- high-confidence taxa"),
         subtitle = "top 15 by prevalence then total RPM", x = "", y = "RPM")
}
if (nrow(hc) > 0) {
  pdf(file.path(FIGS_DIR, "COHORT_high_confidence_landscape.pdf"), width = 12, height = 14)
  for (k in c("Bacteria","Viruses","Eukaryota","Archaea")) op(print(landscape(hc, k)))
  dev.off()
} else message("No high-confidence taxa to plot.")

# ---- 3. Kraken2 vs Kaiju concordance ----------------------------------------
conc <- xv %>% filter(reads > 0 | kaiju_reads > 0)
if (nrow(conc) > 0 && "kaiju_reads" %in% names(conc)) {
  pdf(file.path(FIGS_DIR, "kraken2_vs_kaiju_concordance.pdf"), width = 9, height = 8)
  op(print(
    conc %>% ggplot(aes(reads + 1, kaiju_reads + 1, colour = detected_by)) +
      geom_point(alpha = 0.4, size = 0.8) +
      scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      theme_minimal(base_size = 11) +
      labs(title = "Per sample-species: Kraken2 vs Kaiju read support",
           subtitle = "Points on the diagonal are concordant; Kraken2-only (off-axis) are suspect",
           x = "Kraken2/Bracken reads (+1)", y = "Kaiju reads (+1)", colour = "")))
  det <- conc %>% count(detected_by)
  op(print(
    det %>% ggplot(aes(detected_by, n, fill = detected_by)) + geom_col() +
      theme_minimal(base_size = 12) + theme(legend.position = "none") +
      labs(title = "Detections by classifier agreement", x = "", y = "sample-species calls")))
  dev.off()
}

message("Figures written to: ", FIGS_DIR)
