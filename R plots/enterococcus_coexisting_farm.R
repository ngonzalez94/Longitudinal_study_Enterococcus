#####################################################################
# Enterococcus species co-occurrence per sample - UpSet plot (R version)
#
# Requires: readxl, dplyr, tidyr, UpSetR
#   install.packages(c("readxl", "dplyr", "tidyr", "UpSetR"))
#####################################################################

library(readxl)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(patchwork)

# ------------------------------------------------------------------
# 1. Load data and keep only Enterococcus isolates
# ------------------------------------------------------------------
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/")
data <- read_excel("NGIsolates.xlsx")

data_small <- data %>%
  select(FKSampleID, bact_genus, bact_species)

# Load your sample metadata file (update filename/sheet name to match yours)
sample_meta <- read_excel("NGSamples.xlsx") %>%
  select(FKSampleID, FarmID, SamplingTimepoint)

#------------------------------------------------------------------
# 2. Build presence/absence matrix per sample, with metadata attached
# ------------------------------------------------------------------
entero <- data_small %>%
  filter(bact_genus == "Enterococcus") %>%
  distinct(FKSampleID, bact_species) %>%
  left_join(sample_meta, by = "FKSampleID")

pres <- entero %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    id_cols = c(FKSampleID, FarmID, SamplingTimepoint),
    names_from = bact_species,   # <- takes each VALUE in bact_species (durans, faecalis...)
    values_from = present,       #    and makes it its own COLUMN of TRUE/FALSE
    values_fill = FALSE
  )

species_cols <- setdiff(names(pres), c("FKSampleID", "FarmID", "SamplingTimepoint"))

# ------------------------------------------------------------------
# 3b. UpSet plot coloured by FarmID
# ------------------------------------------------------------------
upset_by_farm <- upset(
  pres,
  species_cols,
  name = "Species",
  min_size = 1,
  width_ratio = 0.15,
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = FarmID)
    ) + labs(fill = "Farm")
  )
) + plot_annotation(title = "Enterococcus species co-occurrence, coloured by farm")

print(upset_by_farm)
ggsave("enterococcus_upset_by_farm.png", upset_by_farm,
       width = 12, height = 7, dpi = 300)
