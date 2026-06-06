#!/usr/bin/env Rscript
# =============================================================================
# 05_aggregate_filter.R
# Merge Kraken2/Bracken + Kaiju, remove host & contaminants, apply confidence
# and abundance filters, cross-validate the two classifiers, and (optionally)
# run decontam against negative controls.
#
# Called by 05_downstream.sbatch, which passes config values as positional args:
#  1 RESULTS_DIR          2 BRACKEN_BIOM(species)   3 KAIJU_MATRIX(species)
#  4 MINIMIZER_STATS      5 TOTAL_READS_LOOKUP      6 READCOUNTS_DIR
#  7 BLOCKLIST            8 MIN_READS               9 MIN_RPM
# 10 MIN_PREVALENCE_FRAC 11 MIN_MINIMIZER_COVERAGE 12 MIN_DISTINCT_MINIMIZERS
# 13 REQUIRE_BOTH        14 CONTROLS_LIST (path or "NONE")
# =============================================================================
suppressPackageStartupMessages({
  # explicit core packages (not the tidyverse meta-package)
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(stringr); library(purrr)
  library(phyloseq); library(biomformat)
})

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 14) stop("Expected 14 arguments; see header.")
RESULTS_DIR <- a[1]; BIOM <- a[2]; KAIJU <- a[3]; MMZ <- a[4]; TOTALS <- a[5]
READCOUNTS <- a[6]; BLOCKLIST <- a[7]
MIN_READS <- as.numeric(a[8]); MIN_RPM <- as.numeric(a[9])
MIN_PREV <- as.numeric(a[10]); MIN_COV <- as.numeric(a[11])
MIN_DISTINCT <- as.numeric(a[12]); REQUIRE_BOTH <- tolower(a[13]) == "true"
CONTROLS <- a[14]

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
log <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# ---- helpers ----------------------------------------------------------------
clean_sample <- function(x) {
  # strip extensions outside-in so e.g. "S1.species_bracken.kreport" -> "S1"
  x %>% basename() %>%
    str_remove("\\.names\\.out$") %>%
    str_remove("\\.kaiju\\.out$") %>%
    str_remove("\\.kreport\\.minimizer$") %>%
    str_remove("\\.kreport$") %>%
    str_remove("\\.species_bracken$") %>% str_remove("\\.genus_bracken$") %>%
    str_remove("_bracken$") %>% str_remove("\\.bracken$") %>%
    str_remove("\\.species$") %>% str_remove("\\.genus$") %>%
    str_remove("_microbe$") %>%
    str_trim()
}
clean_tax <- function(x) x %>% str_remove_all("^[a-z]__") %>%
  str_replace_all(";", "") %>% str_trim()
norm_name <- function(x) x %>% tolower() %>% str_trim()

# ---- 1. Kraken2/Bracken (species) from BIOM ---------------------------------
log("Loading Bracken BIOM:", BIOM)
kb <- tryCatch({
  ps <- import_biom(BIOM)
  m  <- psmelt(ps)
  rk <- paste0("Rank", 1:7)
  m %>%
    rename(reads = Abundance, sample = Sample) %>%
    mutate(taxid = as.character(OTU),
           superkingdom = clean_tax(.data[[rk[1]]]),
           phylum = clean_tax(.data[[rk[2]]]), class = clean_tax(.data[[rk[3]]]),
           order = clean_tax(.data[[rk[4]]]), family = clean_tax(.data[[rk[5]]]),
           genus = clean_tax(.data[[rk[6]]]), species = clean_tax(.data[[rk[7]]]),
           sample = clean_sample(sample)) %>%
    filter(reads > 0) %>%
    select(sample, taxid, superkingdom, phylum, class, order, family, genus, species, reads)
}, error = function(e) { log("WARN: could not read BIOM:", conditionMessage(e)); NULL })

if (is.null(kb) || nrow(kb) == 0) stop("No Kraken2/Bracken data parsed; cannot continue.")

# ---- 2. Total reads (host-input) for RPM ------------------------------------
totals <- tryCatch(
  read_tsv(TOTALS, col_names = c("sample","total_reads"), show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample)),
  error = function(e) tibble(sample = character(), total_reads = numeric()))

kb <- kb %>% left_join(totals, by = "sample") %>%
  group_by(sample) %>%
  mutate(classified_reads = sum(reads),
         total_reads = ifelse(is.na(total_reads), classified_reads, total_reads)) %>%
  ungroup() %>%
  mutate(RPM = reads / total_reads * 1e6,
         percent = reads / classified_reads * 100)

# ---- 3. Distinct-minimizer (KrakenUniq-style) FP filter ---------------------
if (file.exists(MMZ)) {
  mmz <- read_tsv(MMZ, show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample), taxid = as.character(taxid),
           coverage = ifelse(total_minimizers > 0, distinct_minimizers/total_minimizers, 0))
  kb <- kb %>% left_join(select(mmz, sample, taxid, distinct_minimizers, coverage),
                         by = c("sample","taxid"))
  kb <- kb %>% mutate(
    distinct_minimizers = ifelse(is.na(distinct_minimizers), NA_real_, distinct_minimizers),
    pass_minimizer = is.na(distinct_minimizers) |
      (distinct_minimizers >= MIN_DISTINCT & coverage >= MIN_COV))
  n_na_mmz <- sum(is.na(kb$distinct_minimizers))
  log("Applied minimizer filter (distinct>=", MIN_DISTINCT, ", coverage>=", MIN_COV, ")")
  if (n_na_mmz > 0)
    log("NOTE:", n_na_mmz, "taxon-calls had no minimizer row and auto-passed (fail-open);",
        "verify step 02 wrote a *.kreport.minimizer for every sample.")
} else {
  kb <- kb %>% mutate(distinct_minimizers = NA_real_, coverage = NA_real_, pass_minimizer = TRUE)
  log("WARN: no minimizer_stats.tsv; skipping KrakenUniq-style filter.")
}

# ---- 4. Host + contaminant (kitome) flags -----------------------------------
block <- character(0)
if (file.exists(BLOCKLIST)) {
  block <- read_lines(BLOCKLIST) %>% str_remove("#.*$") %>% str_trim()
  block <- norm_name(block[block != ""])
}
is_host <- function(g, s) norm_name(g) == "homo" | norm_name(s) == "homo sapiens"
kb <- kb %>% mutate(
  is_host = is_host(genus, species),
  is_blocklisted = norm_name(genus) %in% block | norm_name(species) %in% block,
  pass_abundance = reads >= MIN_READS & RPM >= MIN_RPM)

# ---- 5. Kaiju (species, full lineage) ---------------------------------------
kj <- tryCatch({
  # kaiju2table columns: file  percent  reads  taxon_id  taxon_name(lineage)
  kj_raw <- read_tsv(KAIJU, show_col_types = FALSE)
  if (!"taxon_id" %in% names(kj_raw)) kj_raw$taxon_id <- NA   # tolerate older kaiju2table output
  kj_raw %>%
    mutate(sample = clean_sample(file), taxid = as.character(taxon_id)) %>%
    separate(taxon_name,
             into = c("superkingdom","phylum","class","order","family","genus","species"),
             sep = ";", extra = "drop", fill = "right") %>%
    mutate(across(c(superkingdom,phylum,class,order,family,genus,species), clean_tax)) %>%
    filter(!is.na(species), species != "", species != "NA", reads > 0) %>%
    group_by(sample, taxid, superkingdom, phylum, class, order, family, genus, species) %>%
    summarise(kaiju_reads = sum(reads), .groups = "drop")
}, error = function(e) { log("WARN: could not read Kaiju matrix:", conditionMessage(e)); NULL })

# ---- 6. Cross-validate Kraken2 vs Kaiju -------------------------------------
# Match on NCBI taxid (robust) OR normalized species name (covers taxonomy-
# version skew between the Kraken2 and Kaiju DBs). A call is concordant ("both")
# if EITHER key matches. Name-only matching silently lost real concordances
# whenever the two DBs spelled a species differently, or Kaiju emitted a
# truncated lineage that pushed the species name out of the last field.
kb2 <- kb %>%
  mutate(sp_key = norm_name(species), ge_key = norm_name(genus)) %>%
  filter(sp_key != "", sp_key != "na")

if (!is.null(kj) && nrow(kj) > 0) {
  kj_name <- kj %>% mutate(sp_key = norm_name(species)) %>%
    filter(sp_key != "", sp_key != "na") %>%
    group_by(sample, sp_key) %>% summarise(kaiju_reads_name = sum(kaiju_reads), .groups = "drop")
  kj_taxid <- kj %>% filter(!is.na(taxid), taxid != "", taxid != "NA") %>%
    group_by(sample, taxid) %>% summarise(kaiju_reads_taxid = sum(kaiju_reads), .groups = "drop")
  xval <- kb2 %>%
    left_join(kj_name,  by = c("sample","sp_key")) %>%
    left_join(kj_taxid, by = c("sample","taxid")) %>%
    mutate(kaiju_reads_name  = replace_na(kaiju_reads_name, 0),
           kaiju_reads_taxid = replace_na(kaiju_reads_taxid, 0),
           kaiju_reads = pmax(kaiju_reads_name, kaiju_reads_taxid),
           detected_by = if_else(kaiju_reads > 0, "both", "kraken2_only")) %>%
    select(-kaiju_reads_name, -kaiju_reads_taxid)
} else {
  xval <- kb2 %>% mutate(kaiju_reads = 0, detected_by = "kraken2_only")
  log("WARN: no Kaiju data; all calls flagged kraken2_only.")
}

# prevalence across samples (Kraken2 detections)
prev <- xval %>% filter(pass_abundance, !is_host, !is_blocklisted, pass_minimizer) %>%
  group_by(sp_key) %>% summarise(n_samples = n_distinct(sample), .groups = "drop")
# Cohort denominator = every sample that entered classification (one
# readcount.tsv per host-removed sample), NOT just samples that happened to
# yield a surviving Bracken call -- otherwise prevalence is inflated and every
# "X / Y samples" figure reports the wrong Y.
cohort_n <- if (dir.exists(READCOUNTS))
  length(list.files(READCOUNTS, pattern = "\\.readcount\\.tsv$")) else 0L
n_total <- max(cohort_n, n_distinct(xval$sample))
xval <- xval %>% left_join(prev, by = "sp_key") %>%
  mutate(n_samples = replace_na(n_samples, 0),
         pass_prevalence = (n_samples / pmax(n_total,1)) >= MIN_PREV)

# final high-confidence flag
xval <- xval %>% mutate(
  high_confidence = pass_abundance & pass_minimizer & pass_prevalence &
    !is_host & !is_blocklisted &
    (if (REQUIRE_BOTH) detected_by == "both" else TRUE))

# ---- 7. Optional decontam against negative controls -------------------------
decontam_tab <- NULL
if (CONTROLS != "NONE" && file.exists(CONTROLS) && requireNamespace("decontam", quietly = TRUE)) {
  log("Running decontam (prevalence) against controls in:", CONTROLS)
  ctl <- read_tsv(CONTROLS, col_names = c("sample","type"), show_col_types = FALSE) %>%
    mutate(sample = clean_sample(sample))
  mat <- xval %>% group_by(sample, taxid) %>% summarise(reads = sum(reads), .groups="drop") %>%
    pivot_wider(names_from = sample, values_from = reads, values_fill = 0) %>%
    column_to_rownames("taxid") %>% as.matrix()
  smp <- intersect(colnames(mat), ctl$sample)
  if (length(smp) > 2 && any(ctl$type[match(smp, ctl$sample)] == "control")) {
    is_ctrl <- ctl$type[match(smp, ctl$sample)] == "control"
    dc <- decontam::isContaminant(t(mat[, smp, drop=FALSE]), neg = is_ctrl, method = "prevalence")
    decontam_tab <- dc %>% rownames_to_column("taxid")
    contam_ids <- decontam_tab %>% filter(contaminant) %>% pull(taxid)
    xval <- xval %>% mutate(is_decontam_contaminant = taxid %in% contam_ids,
                            high_confidence = high_confidence & !is_decontam_contaminant)
    write_csv(decontam_tab, file.path(RESULTS_DIR, "decontam_contaminants.csv"))
    log("decontam flagged", length(contam_ids), "contaminant taxa.")
  } else log("WARN: controls list lacks usable control samples; skipping decontam.")
} else if (CONTROLS != "NONE") {
  log("NOTE: decontam not run (no controls file or 'decontam' package missing).")
}

# ---- 7b. Kaiju-only taxa (detected by Kaiju, NOT by Kraken2/Bracken) --------
# Append them as their own rows so cross_validated_species.csv shows the full
# kraken2_only / kaiju_only / both breakdown. Kraken2/Bracken-derived metrics
# (reads, RPM, percent, minimizers, prevalence) are NA by definition; abundance
# and host/blocklist flags use the Kaiju read count. A Kaiju-only call is
# high-confidence ONLY under REQUIRE_BOTH=false (it can never be "both" and has
# no Kraken2 FP gates), so with the default REQUIRE_BOTH=true these appear in the
# cross-validated table but never in high_confidence_taxa.csv.
n_kaiju_only <- 0L
if (!is.null(kj) && nrow(kj) > 0) {
  kj_sp <- kj %>% mutate(sp_key = norm_name(species), ge_key = norm_name(genus)) %>%
    filter(sp_key != "", sp_key != "na")
  kaiju_only <- kj_sp %>%
    anti_join(distinct(kb2, sample, sp_key), by = c("sample","sp_key")) %>%
    anti_join(distinct(kb2, sample, taxid),  by = c("sample","taxid")) %>%
    mutate(
      reads = NA_real_, classified_reads = NA_real_, total_reads = NA_real_,
      RPM = NA_real_, percent = NA_real_,
      distinct_minimizers = NA_real_, coverage = NA_real_, pass_minimizer = NA,
      is_host        = is_host(genus, species),
      is_blocklisted = norm_name(genus) %in% block | norm_name(species) %in% block,
      pass_abundance = kaiju_reads >= MIN_READS,
      n_samples = NA_integer_, pass_prevalence = NA,
      detected_by = "kaiju_only",
      high_confidence = if (REQUIRE_BOTH) FALSE
                        else (pass_abundance & !is_host & !is_blocklisted))
  n_kaiju_only <- nrow(kaiju_only)
  if (n_kaiju_only > 0)
    xval <- bind_rows(xval, kaiju_only %>% select(any_of(names(xval))))
}
log("Classifier breakdown -- both:", sum(xval$detected_by == "both"),
    " kraken2_only:", sum(xval$detected_by == "kraken2_only"),
    " kaiju_only:", n_kaiju_only)

# Kraken2/Bracken-backbone view (excludes the appended Kaiju-only rows). The
# filter funnel and the prevalence cascade below describe the Kraken2 path only.
xk <- xval %>% filter(detected_by != "kaiju_only")

# ---- 8. Write outputs -------------------------------------------------------
out_cols <- c("sample","taxid","superkingdom","phylum","class","order","family",
              "genus","species","reads","kaiju_reads","total_reads","RPM","percent",
              "distinct_minimizers","coverage","n_samples","detected_by",
              "is_host","is_blocklisted","pass_minimizer","pass_abundance",
              "pass_prevalence","high_confidence")
out_cols <- intersect(out_cols, names(xval))
write_csv(xval %>% select(all_of(out_cols)),
          file.path(RESULTS_DIR, "cross_validated_species.csv"))

hc <- xval %>% filter(high_confidence)
write_csv(hc %>% select(all_of(out_cols)),
          file.path(RESULTS_DIR, "high_confidence_taxa.csv"))

# RPM is Kraken2/Bracken-derived (NA for any kaiju_only member of a group), so
# aggregate it NA-safely; kaiju_reads carries the orthogonal Kaiju support.
rpm_agg <- function(f, x) if (all(is.na(x))) NA_real_ else f(x[!is.na(x)])
cohort <- hc %>% group_by(superkingdom, genus, species) %>%
  summarise(n_samples = n_distinct(sample),
            mean_RPM   = rpm_agg(mean, RPM),
            median_RPM = rpm_agg(stats::median, RPM),
            max_RPM    = rpm_agg(max, RPM),
            kaiju_reads = sum(kaiju_reads, na.rm = TRUE),
            detected_by = paste(sort(unique(detected_by)), collapse=","),
            .groups = "drop") %>% arrange(desc(n_samples), desc(mean_RPM))
write_csv(cohort, file.path(RESULTS_DIR, "cohort_high_confidence_summary.csv"))

# --- Filter funnel: how many sample-taxon calls survive each successive gate ---
# Kraken2/Bracken backbone only (xk); Kaiju-only rows are reported separately.
g_host  <- !xk$is_host
g_block <- g_host  & !xk$is_blocklisted
g_minz  <- g_block & xk$pass_minimizer
g_abund <- g_minz  & xk$pass_abundance
g_prev  <- g_abund & xk$pass_prevalence
# Stage 6 mirrors the FINAL concordance gate: it only narrows when REQUIRE_BOTH
# is on, so the funnel stays monotonic (calls_dropped_here never goes negative).
g_both  <- if (REQUIRE_BOTH) g_prev & (xk$detected_by == "both") else g_prev
hc_k <- xk %>% filter(high_confidence)   # Kraken2-backbone high-confidence (funnel endpoint)
funnel <- tibble(
  stage = c("0. raw sample-taxon calls (Kraken2/Bracken)",
            "1. after host (Homo sapiens) removal",
            "2. after kitome blocklist",
            "3. after distinct-minimizer filter",
            "4. after abundance filter (reads/RPM)",
            "5. after prevalence filter",
            "6. concordant with Kaiju (both classifiers)",
            "7. HIGH-CONFIDENCE (final)"),
  calls = c(nrow(xk), sum(g_host), sum(g_block), sum(g_minz),
            sum(g_abund), sum(g_prev), sum(g_both), nrow(hc_k)),
  distinct_species = c(n_distinct(xk$sp_key), n_distinct(xk$sp_key[g_host]),
            n_distinct(xk$sp_key[g_block]), n_distinct(xk$sp_key[g_minz]),
            n_distinct(xk$sp_key[g_abund]), n_distinct(xk$sp_key[g_prev]),
            n_distinct(xk$sp_key[g_both]), n_distinct(hc_k$sp_key)))
funnel <- funnel %>% mutate(calls_dropped_here = dplyr::lag(calls) - calls)
write_csv(funnel, file.path(RESULTS_DIR, "filter_funnel_summary.csv"))

# --- QC read accounting (per-stage removals + Kraken2 classified/unclassified) ---
qc <- tryCatch({
  files <- list.files(READCOUNTS, pattern="\\.readcount\\.tsv$", full.names=TRUE)
  if (length(files)==0) NULL else map_dfr(files, ~read_tsv(.x, show_col_types = FALSE))
}, error=function(e) NULL)
if (!is.null(qc)) {
  if (!"after_t2t" %in% names(qc)) qc$after_t2t <- qc$after_dedup   # back-compat
  qc <- qc %>% mutate(
    host_removed_STAR         = input_pairs - star_unmapped,
    rRNA_removed              = star_unmapped - after_rrna,
    lowcomplexity_removed     = after_rrna - after_entropy,
    duplicates_removed        = after_entropy - after_dedup,
    t2t_removed               = after_dedup - after_t2t,
    pct_host_STAR             = round(100*(input_pairs - star_unmapped)/pmax(input_pairs,1), 2),
    pct_rRNA_removed          = round(100*(star_unmapped - after_rrna)/pmax(star_unmapped,1), 2),
    pct_lowcomplexity_removed = round(100*(after_rrna - after_entropy)/pmax(after_rrna,1), 2),
    pct_duplicates_removed    = round(100*(after_entropy - after_dedup)/pmax(after_entropy,1), 2),
    pct_t2t_removed           = round(100*(after_dedup - after_t2t)/pmax(after_dedup,1), 2),
    pct_microbial_final       = round(100*final_microbial/pmax(input_pairs,1), 4))
  # fold in Kraken2 per-sample classified/unclassified (written by step 03)
  kqc_path <- file.path(dirname(TOTALS), "kraken2_classification_QC.tsv")
  if (file.exists(kqc_path)) {
    kqc <- tryCatch(read_tsv(kqc_path, show_col_types = FALSE) %>%
                      mutate(sample = clean_sample(sample)), error = function(e) NULL)
    if (!is.null(kqc)) qc <- qc %>% left_join(kqc, by = "sample")
  }
  # fold in Kaiju per-sample classified/unclassified (written by step 04b)
  kjqc_path <- file.path(dirname(KAIJU), "kaiju_classification_QC.tsv")
  if (file.exists(kjqc_path)) {
    kjqc <- tryCatch(read_tsv(kjqc_path, show_col_types = FALSE) %>%
                       mutate(sample = clean_sample(sample)), error = function(e) NULL)
    if (!is.null(kjqc)) qc <- qc %>% left_join(kjqc, by = "sample")
  }
  # per-sample count of high-confidence species (believable microbes per sample)
  qc <- qc %>%
    left_join(hc %>% count(sample, name = "n_high_confidence_species"), by = "sample") %>%
    mutate(n_high_confidence_species = tidyr::replace_na(n_high_confidence_species, 0L))
  write_csv(qc, file.path(RESULTS_DIR, "QC_read_accounting.csv"))
}

# --- One-glance narrative summary (pipeline_summary.txt) ---------------------
med <- function(x) if (length(x)) round(stats::median(suppressWarnings(as.numeric(x)), na.rm = TRUE), 2) else NA
n_samp <- n_total   # cohort size (samples that entered classification), matches the prevalence denominator
con <- file(file.path(RESULTS_DIR, "pipeline_summary.txt"), "w")
wl <- function(...) writeLines(paste0(...), con)
wl("NLPHL tumor-microbiome pipeline — run summary")
wl("=============================================")
wl("Samples analyzed: ", n_samp)
wl("")
wl("READ FUNNEL  (median % of input read pairs per sample):")
if (exists("qc") && !is.null(qc)) {
  wl(sprintf("  host removed (STAR) ........ %s%%", med(qc$pct_host_STAR)))
  wl(sprintf("  rRNA removed ............... %s%%", med(qc$pct_rRNA_removed)))
  wl(sprintf("  low-complexity removed ..... %s%%", med(qc$pct_lowcomplexity_removed)))
  wl(sprintf("  duplicates removed ......... %s%%", med(qc$pct_duplicates_removed)))
  wl(sprintf("  T2T residual host removed .. %s%%", med(qc$pct_t2t_removed)))
  wl(sprintf("  -> final candidate-microbial %s%% of input", med(qc$pct_microbial_final)))
  if ("pct_unclassified" %in% names(qc))     wl(sprintf("  Kraken2 unclassified (median) %s%%", med(qc$pct_unclassified)))
  if ("pct_kaiju_classified" %in% names(qc)) wl(sprintf("  Kaiju classified (median) .. %s%%", med(qc$pct_kaiju_classified)))
} else wl("  (no read-accounting files found)")
wl("")
wl("TAXON-CALL FUNNEL  (Kraken2/Bracken backbone; one call = one species in one sample):")
for (i in seq_len(nrow(funnel)))
  wl(sprintf("  %-46s calls=%-7d species=%d", funnel$stage[i], funnel$calls[i], funnel$distinct_species[i]))
wl("")
wl("CLASSIFIER BREAKDOWN  (all candidate calls, before high-confidence filtering):")
wl(sprintf("  both (Kraken2 & Kaiju) ..... %d", sum(xval$detected_by == "both")))
wl(sprintf("  Kraken2 only ............... %d", sum(xval$detected_by == "kraken2_only")))
wl(sprintf("  Kaiju only ................. %d", sum(xval$detected_by == "kaiju_only")))
wl("")
wl("TOP HIGH-CONFIDENCE SPECIES (by prevalence):")
top <- utils::head(cohort, 15)
if (nrow(top) > 0) for (i in seq_len(nrow(top)))
  wl(sprintf("  %2d. %-32s %d/%d samples  median RPM %.1f  [%s]",
             i, top$species[i], top$n_samples[i], n_samp, top$median_RPM[i], top$detected_by[i]))
wl("")
wl("Key files: high_confidence_taxa.csv | cohort_high_confidence_summary.csv |")
wl("           cross_validated_species.csv | QC_read_accounting.csv | filter_funnel_summary.csv")
close(con)

log("Done. Wrote results to:", RESULTS_DIR)
log("Wrote pipeline_summary.txt (human-readable run summary).")
log(sprintf("High-confidence: %d sample-taxon calls, %d distinct species.",
            nrow(hc), n_distinct(hc$species)))
