#!/usr/bin/env Rscript
# =============================================================================
# install_r_deps.R  --  one-time install of the pipeline's R packages
# =============================================================================
# Package libraries are R-VERSION-SPECIFIC: anything you installed under an older
# R (e.g. 4.1.2) is NOT visible to R 4.4.2. Run this once after loading 4.4.2:
#
#     ml R/4.4.2
#     Rscript install_r_deps.R
#
# Installs into your personal library (R_LIBS_USER, ~/R/.../4.4). Re-running is
# safe; it only installs what is missing.
# =============================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Core tidyverse component packages ONLY -- NOT the `tidyverse` meta-package,
# which Imports ragg/googledrive/reprex/etc. that need cmake + image system libs
# (freetype, libwebp, libuv) absent on the login node and not used by this pipeline.
cran <- c("readr", "dplyr", "tidyr", "tibble", "stringr", "purrr",
          "ggplot2", "scales", "patchwork", "ggrepel")     # CRAN deps
bioc <- c("phyloseq", "biomformat", "decontam")            # Bioconductor deps

lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(lib) && !dir.exists(lib)) dir.create(lib, recursive = TRUE, showWarnings = FALSE)
cat("R", as.character(getRversion()), "| user library:", lib, "\n\n")

missing <- function(p) !requireNamespace(p, quietly = TRUE)

to_cran <- Filter(missing, cran)
if (length(to_cran)) { cat("Installing CRAN:", paste(to_cran, collapse = ", "), "\n")
  install.packages(to_cran) } else cat("CRAN packages already present.\n")

if (missing("BiocManager")) install.packages("BiocManager")
to_bioc <- Filter(missing, bioc)
if (length(to_bioc)) { cat("Installing Bioconductor:", paste(to_bioc, collapse = ", "), "\n")
  BiocManager::install(to_bioc, update = FALSE, ask = FALSE) } else cat("Bioconductor packages already present.\n")

# ---- final report ----
all_pkgs <- c(cran, bioc)
status <- vapply(all_pkgs, function(p) !missing(p), logical(1))
cat("\nDependency status (R", as.character(getRversion()), "):\n")
for (p in all_pkgs) cat(sprintf("  [%s] %s\n", if (status[[p]]) "OK " else "MISS", p))
if (any(!status)) {
  cat("\nSome packages failed to install. Re-run, or check compiler/module errors above.\n")
  quit(status = 1)
}
cat("\nAll dependencies present. You can now submit 05_downstream.sbatch.\n")
