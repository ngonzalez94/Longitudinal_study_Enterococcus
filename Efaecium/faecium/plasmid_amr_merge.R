library(dplyr)
library(readr)
library(readxl)

# ── File paths ────────────────────────────────────────────────────────────────
# Update these paths to point to your actual files
file1_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/amrfinder_faecium.xlsx"   # PlasmidFinder / AMRFinder results (TSV)
file2_path <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/Efaecium_prediction.csv"   # Plasmid prediction results (TSV)
output_path <- "merged_RFPlasmid_amrfinder_Efaecium.tsv"

# ── Read files ────────────────────────────────────────────────────────────────
# File 1: tab-separated, header starts with '#FILE'
df1 <- read_excel(
  file1_path)

# Rename '#FILE' to 'FILE' so dplyr handles it cleanly
names(df1)[names(df1) == "#FILE"] <- "FILE"

# File 2: tab-separated
df2 <- read_csv(
  file2_path,
  col_types = cols(.default = "c")
)

# ── Join on SEQUENCE == contigID ──────────────────────────────────────────────
merged <- inner_join(
  df1 %>% select(Sample, `Gene symbol`, `Element subtype`, Class, ContigID),
  df2 %>% select(prediction, `votes chromosomal`, `votes plasmid`, contigID),
  by = c("ContigID" = "contigID")
) %>%
  select(
    Sample              = Sample,
    prediction,
    `votes chromosomal`,
    `votes plasmid`,
    contigID             = ContigID, 
    `Gene symbol`, 
    `Element subtype`, 
    Class
    
  )

# ── Write output ──────────────────────────────────────────────────────────────
write_tsv(merged, output_path)

cat(sprintf(
  "Done. %d matching rows written to '%s'.\n",
  nrow(merged), output_path
))
