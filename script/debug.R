library(dplyr)
library(readr)
library(stringr)

cytobands <- read_tsv("cytoBand_hg19.txt", col_names = c("chrom", "start", "end", "name", "gieStain"), show_col_types = FALSE)
cnv_seg_data <- read_tsv("combined_summary/cnvkit_aberrations.tsv", show_col_types=FALSE) %>% filter(Chromosome=="chr6")

for (i in 1:nrow(cnv_seg_data)) {
  seg_j <- cnv_seg_data[i, ]
  midpoint <- seg_j$Start + (seg_j$End - seg_j$Start) / 2
  band_hit <- cytobands %>% filter(chrom == seg_j$Chromosome, start <= midpoint, end >= midpoint)
  clean_chr <- str_remove(seg_j$Chromosome, "^chr")
  exact_band <- paste0(clean_chr, band_hit$name[1])
  arm_only <- paste0(clean_chr, substr(band_hit$name[1], 1, 1))
  cat(sprintf("Sample: %s, Range: %d-%d, Midpoint: %d, Exact: %s, Arm: %s\n", 
              seg_j$Sample, seg_j$Start, seg_j$End, as.integer(midpoint), exact_band, arm_only))
}
