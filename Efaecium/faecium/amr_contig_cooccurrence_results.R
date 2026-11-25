# =============================================================================
# AMR Gene Co-occurrence Analysis
# Goal: Count isolates sharing a specific gene combination on the SAME contig,
#       report antibiotic classes, and express results as % and fraction n/N.
# =============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# setwd
setwd("")
df <- read_xlsx("amrfinder_faecium.xlsx")

# Normalise column names (trim whitespace)
names(df) <- trimws(names(df))

# Keep only AMR rows (ignore STRESS, VIRULENCE, etc. if present)
df <- df %>% filter(Type == "AMR")

# -----------------------------------------------------------------------------
# 2. Build a per-contig gene table
#    Key columns: Sample (isolate), ContigID, Gene symbol, Class
# -----------------------------------------------------------------------------
contig_genes <- df %>%
  select(Sample, ContigID, gene = `Gene symbol`, class = Class) %>%
  distinct()

# -----------------------------------------------------------------------------
# 3. For every contig that carries ≥2 genes, collapse into a gene-combination
#    row, keeping track of which AMR classes those genes belong to.
# -----------------------------------------------------------------------------
contig_combos <- contig_genes %>%
  group_by(Sample, ContigID) %>%
  summarise(
    n_genes        = n_distinct(gene),
    gene_combo     = paste(sort(unique(gene)),  collapse = " | "),
    class_combo    = paste(sort(unique(class)), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_genes >= 2)   # only contigs that co-carry multiple AMR genes

# -----------------------------------------------------------------------------
# 4. Count isolates per gene combination
# -----------------------------------------------------------------------------
total_isolates <- n_distinct(df$Sample)

summary_table <- contig_combos %>%
  group_by(gene_combo, class_combo) %>%
  summarise(
    n_isolates   = n_distinct(Sample),
    .groups      = "drop"
  ) %>%
  mutate(
    pct          = round(n_isolates / total_isolates * 100, 1),
    fraction     = paste0(n_isolates, "/", total_isolates)
  ) %>%
  arrange(desc(n_isolates)) %>%
  rename(
    `Gene combination`  = gene_combo,
    `AMR class(es)`     = class_combo,
    `N isolates`        = n_isolates,
    `% isolates`        = pct,
    `Fraction (n/N)`    = fraction
  )

# -----------------------------------------------------------------------------
# 5. Print results
# -----------------------------------------------------------------------------
cat("\n=== AMR gene co-occurrence on the same contig ===\n")
cat("Total isolates in dataset:", total_isolates, "\n\n")
print(summary_table, n = Inf)

# -----------------------------------------------------------------------------
# 6. (Optional) Export to TSV
# -----------------------------------------------------------------------------
write.table(summary_table,
            file      = "amr_contig_cooccurrence_results.tsv",
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)
cat("\nResults saved to: amr_contig_cooccurrence_results.tsv\n")
