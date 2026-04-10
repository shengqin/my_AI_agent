library(data.table)
library(dplyr)

# Define data directory
data_dir <- "/Users/Brian/Library/CloudStorage/Box-Box/BinkleyLab/cHL/cHL_genotyping/data"

# Define the 4 target files
files_to_read <- c(
  "2022-05-20_Stefan1_SNV_detail.txt",
  "2022-09-20_Stefan2_SNV_detail.txt",
  "2023-03-30_Binkley_SNV_detail.txt",
  "2025-11-07_France_SNV_detail.txt"
)

# Safe reading function to avoid memory issues and column mismatches
read_snv_safely <- function(file_name) {
  file_path <- file.path(data_dir, file_name)
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Load the entire file using fast data.table
  # ADDED quote="" string. Otherwise fread can "swallow" thousands of rows if there's a rogue quote in a gene description!
  df <- fread(file_path, data.table = FALSE, quote = "")
  
  # Harmonize variant naming conventions (Looking at you, France Dataset)
  if ("Variant class" %in% colnames(df)) df$Variantclass <- df$`Variant class`
  if ("Variant type" %in% colnames(df)) df$Varianttype <- df$`Variant type`
  
  # Extract only the needed columns to save memory and avoid binding type issues
  needed_cols <- c("Sample", "Variant_Classification", "Varianttype", "VAF")
  available_cols <- needed_cols[needed_cols %in% colnames(df)]
  df_sub <- df[, available_cols, drop = FALSE]
  
  # Inject NAs if any specific naming standard wasn't matched
  for (col in needed_cols) {
    if (!col %in% colnames(df_sub)) {
      df_sub[[col]] <- NA
    }
  }
  
  df_sub$Dataset <- file_name
  return(df_sub)
}

# Bind the data from all files
snv_data <- bind_rows(lapply(files_to_read, read_snv_safely))

# We actively search for terms describing HARMFUL coding mutations
harmful_keywords <- c("missense", "nonsense", "stop", "splice", "startloss", "frameshift")

# Filter out non-harmful mutations and compute VAF stats per sample
sample_vaf_stats <- snv_data %>%
  mutate(
    # Merge classification types for a robust search text
    search_text = paste(tolower(Variant_Classification), tolower(Varianttype), sep = " ")
  ) %>%
  # KEEP only those matching harmful keywords
  filter(grepl(paste(harmful_keywords, collapse = "|"), search_text)) %>%
  # STRIP out any "near-splice" variants that are still synonymous or intronic
  filter(!grepl("synonymous|intron", search_text)) %>%
  # Force VAF as numeric for safe means
  mutate(VAF = as.numeric(VAF)) %>%
  filter(!is.na(VAF)) %>%
  # Group by sample ID and summarize
  group_by(Sample) %>%
  summarise(
    Harmful_Mutation_Count = n(),
    Mean_VAF = mean(VAF),
    Median_VAF = median(VAF),
    .groups = "drop"
  )

library(tidyr)

# Parse Patient mapping and Source type from Sample ID
sample_vaf_stats <- sample_vaf_stats %>%
  mutate(
    Patient = gsub("-.*", "", Sample),
    Source = case_when(
      grepl("-P", Sample) ~ "Plasma",
      grepl("-T", Sample) ~ "Tumor",
      TRUE ~ "Other"
    )
  )

# Aggregate multiple samples (e.g. P1 and P1v2) into a Patient-Level summary
patient_vaf_stats <- sample_vaf_stats %>%
  group_by(Patient, Source) %>%
  summarise(
    # Average the mean VAFs if there are multiple sample replicate runs
    Aggregated_Mean_VAF = mean(Mean_VAF, na.rm = TRUE),
    Aggregated_Median_VAF = median(Median_VAF, na.rm = TRUE),
    # Combine total unique harmful mutations (or sum them if considered distinct)
    Total_Harmful_Mutations = sum(Harmful_Mutation_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot wide so Plasma and Tumor stats are strictly side-by-side per patient
  pivot_wider(
    names_from = Source,
    values_from = c(Aggregated_Mean_VAF, Aggregated_Median_VAF, Total_Harmful_Mutations),
    names_glue = "{Source}_{.value}"
  ) %>%
  # Clean up column names for neatness
  rename_with(~ gsub("Aggregated_", "", .x))

print("\n--- PATIENT LEVEL VAF SUMMARY ---")
print(head(patient_vaf_stats))

library(ggplot2)
library(ggpubr)

fractions_file <- "/Users/Brian/Library/CloudStorage/Box-Box/BrianAnalysis/Lymphoma/cybersort/data/2025-Q3Q4/2025-07-12_cHL_Fractions_Norm.txt"
nomenclature_file <- "/Users/Brian/Library/CloudStorage/Box-Box/BrianAnalysis/Lymphoma/02-2023_Collaboration/binkley_nomenclature_HL_numbers.txt"

if (file.exists(fractions_file) & file.exists(nomenclature_file)) {
  # Load mapping and cell fractions
  fractions <- fread(fractions_file, data.table = FALSE)
  mapping <- fread(nomenclature_file, data.table = FALSE)

  # Harmonize the names using the binkley nomenclature
  fractions <- fractions %>%
    left_join(mapping, by = c("Mixture" = "Binkley_id")) %>%
    mutate(Target_Patient = ifelse(!is.na(AALab_id), AALab_id, Mixture))

  # Join with patient VAFs and determine definitive VAF dynamically
  plot_df <- patient_vaf_stats %>%
    inner_join(fractions, by = c("Patient" = "Target_Patient")) %>%
    # If a person has both plasma or tumor, use plasma data automatically utilizing dplyr::coalesce
    mutate(Final_VAF_Median = coalesce(Plasma_Median_VAF, Tumor_Median_VAF),
           Final_VAF_Mean = coalesce(Plasma_Mean_VAF, Tumor_Mean_VAF),) %>%
    # Filter valid pairs
    filter(!is.na(Final_VAF_Median) & !is.na(HRS))

  # Generate intuitive ggplot with Pearson correlation overlays
  p1 <- ggplot(plot_df, aes(x = Final_VAF_Median, y = HRS)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Fraction vs Patient Median VAF",
      subtitle = "Automatically Prioritizing Plasma VAF for matched collections",
      x = "Median VAF",
      y = "Cybersort HRS Cell Fraction"
    )
  
  p2 <- ggplot(plot_df, aes(x = Final_VAF_Mean, y = HRS)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Fraction vs Patient Median VAF",
      subtitle = "Automatically Prioritizing Plasma VAF for matched collections",
      x = "Mean VAF",
      y = "Cybersort HRS Cell Fraction"
    )
  p3 <- ggplot(plot_df, aes(x = Final_VAF_Median, y = HRS + 0.001)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Fraction vs Patient Median VAF (Log10)",
      subtitle = "Automatically Prioritizing Plasma VAF for matched collections",
      x = "Median VAF (log10)",
      y = "Cybersort HRS Cell Fraction (log10)"
    )

  p4 <- ggplot(plot_df, aes(x = Final_VAF_Mean, y = HRS + 0.001)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    scale_x_log10() +
    scale_y_log10() +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Fraction vs Patient Mean VAF (Log10)",
      subtitle = "Automatically Prioritizing Plasma VAF for matched collections",
      x = "Mean VAF (log10)",
      y = "Cybersort HRS Cell Fraction (log10)"
    )

  plot_df_tumor <- plot_df %>% filter(!is.na(Tumor_Median_VAF))

  p5 <- ggplot(plot_df_tumor, aes(x = Tumor_Median_VAF, y = HRS)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Percentage vs Patient Median VAF",
      subtitle = "Strictly Tumor Samples Only (Linear Scale)",
      x = "Tumor Median VAF",
      y = "Cybersort HRS Cell Fraction (%)"
    )

  p6 <- ggplot(plot_df_tumor, aes(x = Tumor_Mean_VAF, y = HRS)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "HRS Cell Percentage vs Patient Mean VAF",
      subtitle = "Strictly Tumor Samples Only (Linear Scale)",
      x = "Tumor Mean VAF",
      y = "Cybersort HRS Cell Fraction (%)"
    )
  # Pre-transform explicit log variables into the actual dataframe to force pearson correlation calculation upon logarithmic outputs
  plot_df_log_tumor <- plot_df %>%
    filter(!is.na(Tumor_Median_VAF)) %>%
    mutate(
      Log_Tumor_VAF_Median = log10(Tumor_Median_VAF),
      Log_Tumor_VAF_Mean = log10(Tumor_Mean_VAF),
      Log_HRS = log10(HRS + 0.001)  # safely applying pseudocount inherently within vector computation space
    )

  p7 <- ggplot(plot_df_log_tumor, aes(x = Log_Tumor_VAF_Median, y = Log_HRS)) +
    geom_point(alpha = 0.7, color = "darkblue", size = 3) +
    geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Log10(HRS Fraction) vs Log10(Tumor Median VAF)",
      subtitle = "Pearson correlation mapping strictly against Tumor biopsies",
      x = "Log10(Tumor Median VAF)",
      y = "Log10(Cybersort HRS Cell Fraction)"
    )

  p8 <- ggplot(plot_df_log_tumor, aes(x = Log_Tumor_VAF_Mean, y = Log_HRS)) +
    geom_point(alpha = 0.7, color = "darkgrey", size = 3) +
    geom_smooth(method = "lm", color = "blue", , se=FALSE) +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Log10(HRS) vs Log10(Tumor Mean VAF)",
      subtitle = "Pearson correlation mapping strictly against Tumor biopsies",
      x = "Log10(Tumor Mean VAF)",
      y = "Log10(Cybersort HRS Cell Fraction)"
    )

  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  print(p7)
  print(p8)
}
