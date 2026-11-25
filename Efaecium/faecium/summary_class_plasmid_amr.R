# summarise_plasmid_AMR.R
# Filters plasmid contigs (prediction == "p"), then for each unique contigID:
#   - lists all Gene symbols found (Genes_contained)
#   - pivots each unique Class into its own numbered column (Class_1, Class_2, ...)
#   - counts total distinct classes (Total_classes)

library(dplyr)
library(readxl)
library(tidyr)
library(readr)

# ── File paths ────────────────────────────────────────────────────────────────
input_path  <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/merged_RFPlasmid_amrfinder_Efaecium.xlsx"   # <- update this
output_path <- "summarised_plasmid_AMR_Efaecium.tsv"

# ── Read & filter ─────────────────────────────────────────────────────────────
df <- read_excel(input_path, col_types = "text")
names(df) <- trimws(names(df))

df_plasmid <- df %>%
  filter(prediction == "p")

# ── Summarise per Sample + contigID ───────────────────────────────────────────
summarised <- df_plasmid %>%
  group_by(Sample, contigID) %>%
  summarise(
    prediction          = first(prediction),
    `votes chromosomal` = first(`votes chromosomal`),
    `votes plasmid`     = first(`votes plasmid`),
    Genes_contained     = paste(unique(`Gene symbol`), collapse = ", "),
    classes             = list(unique(Class)),
    Total_classes       = n_distinct(Class),
    .groups = "drop"
  )

# ── Pivot classes into Class_1, Class_2, ... columns ─────────────────────────
max_classes <- max(lengths(summarised$classes))

class_cols <- summarised %>%
  mutate(class_list = lapply(classes, function(x) {
    length(x) <- max_classes   # pad shorter lists with NA
    x
  })) %>%
  pull(class_list) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  setNames(paste0("Class_", seq_len(max_classes)))

# ── Combine & write ───────────────────────────────────────────────────────────
result <- summarised %>%
  select(-classes) %>%
  bind_cols(class_cols) %>%
  select(
    Sample,
    prediction,
    `votes chromosomal`,
    `votes plasmid`,
    contigID,
    Genes_contained,
    starts_with("Class_"),
    Total_classes
  )

write_tsv(result, output_path)
cat(sprintf("Done. %d rows written to '%s'.\n", nrow(result), output_path))
