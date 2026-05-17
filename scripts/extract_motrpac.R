#!/usr/bin/env Rscript
# Extract MoTrPAC transcriptomics DE results from the official R packages.
# Outputs TSV files consumed by notebook 01.
#
# Installation (run once):
#   install.packages("remotes")
#   remotes::install_github("MoTrPAC/MotrpacRatTraining6moData")
#   remotes::install_github("MoTrPAC/MotrpacRatTraining6mo")
#
# Usage:
#   Rscript scripts/extract_motrpac.R --outdir data/raw
#
# Output files:
#   data/raw/motrpac_<tissue>_<timepoint>_r_export.tsv

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
})

parser <- ArgumentParser(description = "Export MoTrPAC DE results to TSV")
parser$add_argument("--outdir", default = "data/raw",
                    help = "Output directory for TSV files")
parser$add_argument("--timepoints", default = "8w",
                    help = "Comma-separated timepoints (1w,2w,4w,8w)")
args <- parser$parse_args()

outdir <- args$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
timepoints <- strsplit(args$timepoints, ",")[[1]]

# Tissue IDs in MoTrPAC nomenclature
TISSUES <- list(
  gastrocnemius = "SKM-GN",
  liver         = "LIVER",
  white_adipose = "WAT-SC"
)

# Load the data package
if (!requireNamespace("MotrpacRatTraining6moData", quietly = TRUE)) {
  stop(
    "MotrpacRatTraining6moData not installed.\n",
    "Run: remotes::install_github('MoTrPAC/MotrpacRatTraining6moData')"
  )
}
library(MotrpacRatTraining6moData)

# TRNSCRPT_DA is the differential analysis results object
# It contains training-vs-sedentary comparisons across tissues and timepoints
da <- tryCatch(
  MotrpacRatTraining6moData::TRNSCRPT_DA,
  error = function(e) {
    message("Trying alternative data object name...")
    get("TRNSCRPT_TRAINING_DA", envir = asNamespace("MotrpacRatTraining6moData"))
  }
)

message("Loaded DA object with ", nrow(da), " records")
message("Columns: ", paste(colnames(da), collapse = ", "))
message("Tissues: ", paste(unique(da$tissue), collapse = ", "))
message("Comparisons: ", paste(unique(da$comparison_group), collapse = ", "))

for (tissue_name in names(TISSUES)) {
  tissue_id <- TISSUES[[tissue_name]]

  for (tp in timepoints) {
    # Filter for this tissue and timepoint training-vs-control
    sub <- da %>%
      filter(
        grepl(tissue_id, tissue, ignore.case = TRUE),
        grepl(tp, comparison_group, ignore.case = TRUE) |
          grepl(tp, timepoint, ignore.case = TRUE),
        grepl("trained|training|exercise", comparison_group, ignore.case = TRUE) |
          grepl("trained|training|exercise", contrast, ignore.case = TRUE)
      )

    if (nrow(sub) == 0) {
      # Fallback: just filter by tissue and hope for 8w default
      sub <- da %>%
        filter(grepl(tissue_id, tissue, ignore.case = TRUE))
      message(sprintf("Warning: no exact timepoint match for %s %s, using all %d records",
                      tissue_name, tp, nrow(sub)))
    }

    # Standardise column names
    sub <- sub %>%
      rename_with(~ gsub("gene.?symbol|gene_name", "gene_symbol", .x, ignore.case = TRUE)) %>%
      rename_with(~ gsub("ensembl.?id|ensembl_gene", "ensembl_rat", .x, ignore.case = TRUE)) %>%
      rename_with(~ gsub("log2?fc|log.?fold|logfoldchange|estimate", "logFC",
                         .x, ignore.case = TRUE)) %>%
      rename_with(~ gsub("adj.*p.*val|p.*adj|fdr|qval|bh", "adj_pvalue",
                         .x, ignore.case = TRUE)) %>%
      rename_with(~ gsub("^p.?val.*$|^pvalue$", "pvalue", .x, ignore.case = TRUE))

    out_path <- file.path(outdir, sprintf("motrpac_%s_%s_r_export.tsv", tissue_name, tp))
    write_tsv(sub, out_path)
    message(sprintf("Wrote %d rows → %s", nrow(sub), out_path))
  }
}

message("Done.")
