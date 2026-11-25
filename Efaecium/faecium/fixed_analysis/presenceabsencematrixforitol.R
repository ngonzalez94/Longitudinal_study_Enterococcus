# =============================================================================
# AMR Gene Presence / Absence Matrix
# Values: 1 = present, -1 = absent
# Rows = isolates (derived from the filename prefix in column 1)
# Columns = unique Element symbols across all files
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)

# -----------------------------------------------------------------------------
# 1. Load data
#    The file has a compound first column: "filename,Protein id"
#    e.g. "NG0016bc_contigs_amr.tsv,NA"
#    We read everything as raw text, then parse.
# -----------------------------------------------------------------------------
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/fixed_analysis")
# ---- Option A: single concatenated file (as shown in your example) ----------
raw <- read.delim("amr_faecium.tsv",
                  sep       = "\t",
                  header    = TRUE,
                  stringsAsFactors = FALSE,
                  check.names = FALSE)

# The first column contains "filename,ProteinID" — split it
first_col <- names(raw)[1]   # e.g. "NG0016bc_contigs_amr.tsv,Protein id"

# Extract isolate name from the filename prefix (everything before "_contigs")
raw <- raw %>%
  mutate(
    isolate = str_extract(.[[1]], "^[^,]+"),          # grab filename portion
    isolate = tools::file_path_sans_ext(isolate),     # drop .tsv extension
    isolate = str_remove(isolate, "_contigs_amr$")    # clean suffix
  )

# Keep only AMR rows and relevant columns
df <- raw %>%
  filter(Type == "AMR") %>%
  select(isolate, gene = `Element symbol`) %>%
  distinct()

# -----------------------------------------------------------------------------
# ---- Option B: read multiple TSV files from a folder ------------------------
# Uncomment this block (and comment out Option A above) if you have separate
# per-isolate files in a directory.
# -----------------------------------------------------------------------------
# file_dir <- "path/to/your/amr/files"
# tsv_files <- list.files(file_dir, pattern = "_contigs_amr\\.tsv$",
#                         full.names = TRUE)
#
# df <- lapply(tsv_files, function(f) {
#   d <- read.delim(f, sep = "\t", header = TRUE,
#                   stringsAsFactors = FALSE, check.names = FALSE)
#   # skip the compound first-col issue if files are clean per-sample
#   d$isolate <- tools::file_path_sans_ext(basename(f))
#   d$isolate <- sub("_contigs_amr$", "", d$isolate)
#   d
# }) %>%
#   bind_rows() %>%
#   filter(Type == "AMR") %>%
#   select(isolate, gene = `Element symbol`) %>%
#   distinct()

# -----------------------------------------------------------------------------
# 2. Build presence / absence matrix
#    1  = gene present in that isolate
#    -1 = gene absent  in that isolate
# -----------------------------------------------------------------------------
pa_matrix <- df %>%
  mutate(present = 1L) %>%
  pivot_wider(
    names_from  = gene,
    values_from = present,
    values_fill = -1L          # fill missing = absent
  ) %>%
  arrange(isolate)

# -----------------------------------------------------------------------------
# 3. Print
# -----------------------------------------------------------------------------
cat("\n=== AMR Gene Presence / Absence Matrix (1 = present, -1 = absent) ===\n")
cat("Isolates:", nrow(pa_matrix), "\n")
cat("Genes   :", ncol(pa_matrix) - 1, "\n\n")
print(pa_matrix)

# -----------------------------------------------------------------------------
# 4. Export
# -----------------------------------------------------------------------------
write.table(pa_matrix,
            file      = "amr_presence_absence_matrix.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

cat("\nSaved to: amr_presence_absence_matrix.tsv\n")
