#####################################################################
# Enterococcus plate distribution + AMR resistance analysis
#####################################################################

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/")

data <- read_excel("NGIsolates.xlsx")

# ------------------------------------------------------------------
# Species colour palette (matches requested colour scheme)
# ------------------------------------------------------------------
species_colors <- c(
  "faecalis" = "#AEE1F0",  # light blue
  "faecium"  = "salmon"   # light pink
)

# ------------------------------------------------------------------
# 2. Derive Plate from IsolateID suffix, keep Enterococcus only,
#    join Farm in
# ------------------------------------------------------------------
entero <- data %>%
  mutate(
    suffix = str_sub(IsolateID, -2, -1),
    Plate = case_when(
      suffix %in% c("ba", "bb") ~ "No selection",
      suffix %in% c("bc", "bd") ~ "Tylosin",
      suffix %in% c("be", "bf") ~ "Tetracycline",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(bact_genus == "Enterococcus", !is.na(Plate)) %>%
  left_join(sample_meta, by = "FKSampleID")

entero$Plate <- factor(entero$Plate,
                       levels = c("No selection", "Tylosin", "Tetracycline"))

# ------------------------------------------------------------------
# 3. Count isolates per species, farm, and plate
# ------------------------------------------------------------------
plate_farm_counts <- entero %>%
  count(FarmID, Plate, bact_species)

# ------------------------------------------------------------------
# 4. 3x3 grid: rows = Farm, columns = Plate
# ------------------------------------------------------------------
p_plate_farm <- ggplot(plate_farm_counts, aes(x = bact_species, y = n, fill = bact_species)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.3, size = 3) +
  facet_grid(FarmID ~ Plate, scales = "free_x") +
  labs(
    title = "Enterococcus isolates per species, by farm and plate",
    x = "Species",
    y = "Number of isolates",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_plate_farm)

ggsave("enterococcus_by_farm_and_plate.png", p_plate_farm,
       width = 14, height = 10, dpi = 300)

# ------------------------------------------------------------------
# 4. Resistance calculation using ECOFF cutoffs
#    (E. faecalis and E. faecium only, milk samples excluded)
# ------------------------------------------------------------------

# --- Species-specific TYL cutoff -----------------------------------
# ECOFF: faecium TYL <= 11, faecalis TYL <= 17
# Other antibiotics use the same cutoff for both species.

amr <- data %>%
  filter(bact_genus == "Enterococcus",
         bact_species %in% c("faecalis", "faecium")) %>%
  mutate(
    ERY = ERYdiff,
    TYL = TYLdiff,
    VAN = VANdiff,
    CIP = CIPdiff,
    AMP = AMPdiff,
    TET = TETdiff
  )

# ---- TODO: join sample metadata (CalfID / farm / timepoint) here ----
# Once you send the metadata file, this becomes, e.g.:
#
 sample_meta <- read_excel("NGSamples.xlsx")   # has FKSampleID, CalfID, Farm, Timepoint

 amr <- amr %>%
   left_join(sample_meta, by = "FKSampleID")

 amr <- amr %>%
   filter(!str_detect(str_to_lower(CalfID), "milk"))   # excludes "milk", "MILK", "milk and pellet" etc.
# ----------------------------------------------------------------------

amr <- amr %>%
  mutate(
    ERY_R = ifelse(ERY <= 13.01, 1, 0),
    TYL_R = case_when(
      bact_species == "faecium"  & TYL <= 11.01 ~ 1,
      bact_species == "faecium"  & TYL >  11.01 ~ 0,
      bact_species == "faecalis" & TYL <= 17.01 ~ 1,
      bact_species == "faecalis" & TYL >  17.01 ~ 0,
      TRUE ~ NA_real_
    ),
    VAN_R = ifelse(VAN <= 14.01, 1, 0),
    CIP_R = ifelse(CIP <= 15.01, 1, 0),
    AMP_R = ifelse(AMP <= 16.01, 1, 0),
    TET_R = ifelse(TET <= 14.01, 1, 0)
  )

# Count resistant isolates per species per antimicrobial
resistance_summary <- amr %>%
  group_by(bact_species) %>%
  summarise(
    n_isolates   = n(),
    ERY_resistant = sum(ERY_R, na.rm = TRUE),
    TYL_resistant = sum(TYL_R, na.rm = TRUE),
    VAN_resistant = sum(VAN_R, na.rm = TRUE),
    CIP_resistant = sum(CIP_R, na.rm = TRUE),
    AMP_resistant = sum(AMP_R, na.rm = TRUE),
    TET_resistant = sum(TET_R, na.rm = TRUE)
  )

print(resistance_summary)

write.csv(resistance_summary, "faecalis_faecium_resistance_summary.csv", row.names = FALSE)

# Optional: long format for plotting resistance rates
resistance_long <- amr %>%
  select(bact_species, ERY_R, TYL_R, VAN_R, CIP_R, AMP_R, TET_R) %>%
  pivot_longer(cols = -bact_species, names_to = "Antimicrobial", values_to = "Resistant") %>%
  mutate(Antimicrobial = str_remove(Antimicrobial, "_R")) %>%
  group_by(bact_species, Antimicrobial) %>%
  summarise(
    n_tested   = sum(!is.na(Resistant)),
    n_resistant = sum(Resistant, na.rm = TRUE),
    pct_resistant = 100 * n_resistant / n_tested,
    .groups = "drop"
  )

p_resistance <- ggplot(resistance_long, aes(x = Antimicrobial, y = pct_resistant, fill = bact_species)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n_resistant_samples),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5) +
  labs(
    title = "Percentage of resistant isolates by antimicrobial",
    x = "Antimicrobial",
    y = "% resistant",
    fill = "Species"
  ) +
  theme_bw()

print(p_resistance)
ggsave("resistance_by_antimicrobial.png", p_resistance, width = 8, height = 5, dpi = 300)



#now as sample level
# ---------------------------------------------------------------
# Collapse from ISOLATE level to SAMPLE level (FKSampleID):
# a sample counts as "resistant" to a drug if AT LEAST ONE of its
# isolates (of that species) is resistant to that drug.
# ---------------------------------------------------------------
sample_level <- amr %>%
  group_by(FKSampleID, bact_species) %>%
  summarise(
    n_isolates   = n(),
    ERY_R = as.integer(any(ERY_R == 1, na.rm = TRUE)),
    TYL_R = as.integer(any(TYL_R == 1, na.rm = TRUE)),
    VAN_R = as.integer(any(VAN_R == 1, na.rm = TRUE)),
    CIP_R = as.integer(any(CIP_R == 1, na.rm = TRUE)),
    AMP_R = as.integer(any(AMP_R == 1, na.rm = TRUE)),
    TET_R = as.integer(any(TET_R == 1, na.rm = TRUE)),
    .groups = "drop"
  )

# Number of SAMPLES (not isolates) positive for each antimicrobial, per species
resistance_summary <- sample_level %>%
  group_by(bact_species) %>%
  summarise(
    n_samples     = n(),
    ERY_resistant_samples = sum(ERY_R, na.rm = TRUE),
    TYL_resistant_samples = sum(TYL_R, na.rm = TRUE),
    VAN_resistant_samples = sum(VAN_R, na.rm = TRUE),
    CIP_resistant_samples = sum(CIP_R, na.rm = TRUE),
    AMP_resistant_samples = sum(AMP_R, na.rm = TRUE),
    TET_resistant_samples = sum(TET_R, na.rm = TRUE)
  )

print(resistance_summary)

write.csv(resistance_summary, "faecalis_faecium_resistance_summary.csv", row.names = FALSE)

# Optional: long format for plotting % of samples resistant
resistance_long <- sample_level %>%
  select(bact_species, ERY_R, TYL_R, VAN_R, CIP_R, AMP_R, TET_R) %>%
  pivot_longer(cols = -bact_species, names_to = "Antimicrobial", values_to = "Resistant") %>%
  mutate(Antimicrobial = str_remove(Antimicrobial, "_R")) %>%
  group_by(bact_species, Antimicrobial) %>%
  summarise(
    n_samples_tested   = sum(!is.na(Resistant)),
    n_resistant_samples = sum(Resistant, na.rm = TRUE),
    pct_resistant = 100 * n_resistant_samples / n_samples_tested,
    .groups = "drop"
  )

p_resistance <- ggplot(resistance_long, aes(x = Antimicrobial, y = n_resistant_samples, fill = bact_species)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n_resistant_samples),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(
    values = species_colors,
    labels = c("faecalis" = "E. faecalis", "faecium" = "E. faecium")
  ) +
  labs(
    title = "Number of samples resistant by antimicrobial",
    x = "Antimicrobial",
    y = "Number of resistant samples",
    fill = "Species"
  ) +
  theme_bw()

print(p_resistance)
ggsave("resistance_by_antimicrobial.png", p_resistance, width = 8, height = 5, dpi = 300)


###NOW IT MUST BE SPLIT BY FARM
# 1. Bring FarmID into the sample-level resistance data
#    (assumes `amr` already has FarmID joined in from sample_meta -
#    if not, join it here first: amr <- amr %>% left_join(sample_meta, by = "FKSampleID"))
# ------------------------------------------------------------------
sample_level_farm <- amr %>%
  group_by(FKSampleID, FarmID, bact_species) %>%
  summarise(
    n_isolates = n(),
    ERY_R = as.integer(any(ERY_R == 1, na.rm = TRUE)),
    TYL_R = as.integer(any(TYL_R == 1, na.rm = TRUE)),
    VAN_R = as.integer(any(VAN_R == 1, na.rm = TRUE)),
    CIP_R = as.integer(any(CIP_R == 1, na.rm = TRUE)),
    AMP_R = as.integer(any(AMP_R == 1, na.rm = TRUE)),
    TET_R = as.integer(any(TET_R == 1, na.rm = TRUE)),
    .groups = "drop"
  )

# ------------------------------------------------------------------
# 2. Helper: build a "Farm x Antimicrobial" resistance table for one species
#    Cell format: "XX.XX% (n/N)"
# ------------------------------------------------------------------
build_farm_resistance_table <- function(sp) {
  
  df <- sample_level_farm %>% filter(bact_species == sp)
  
  # Per-farm counts
  by_farm <- df %>%
    group_by(FarmID) %>%
    summarise(
      n_samples = n(),
      ERY = sprintf("%.2f%% (%d/%d)", 100 * sum(ERY_R) / n(), sum(ERY_R), n()),
      TYL = sprintf("%.2f%% (%d/%d)", 100 * sum(TYL_R) / n(), sum(TYL_R), n()),
      VAN = sprintf("%.2f%% (%d/%d)", 100 * sum(VAN_R) / n(), sum(VAN_R), n()),
      CIP = sprintf("%.2f%% (%d/%d)", 100 * sum(CIP_R) / n(), sum(CIP_R), n()),
      AMP = sprintf("%.2f%% (%d/%d)", 100 * sum(AMP_R) / n(), sum(AMP_R), n()),
      TET = sprintf("%.2f%% (%d/%d)", 100 * sum(TET_R) / n(), sum(TET_R), n()),
      .groups = "drop"
    ) %>%
    rename(Farm = FarmID, Samples = n_samples)
  
  # Total row (all farms combined, for this species)
  total_row <- df %>%
    summarise(
      Farm = "Total",
      Samples = n(),
      ERY = sprintf("%.2f%% (%d/%d)", 100 * sum(ERY_R) / n(), sum(ERY_R), n()),
      TYL = sprintf("%.2f%% (%d/%d)", 100 * sum(TYL_R) / n(), sum(TYL_R), n()),
      VAN = sprintf("%.2f%% (%d/%d)", 100 * sum(VAN_R) / n(), sum(VAN_R), n()),
      CIP = sprintf("%.2f%% (%d/%d)", 100 * sum(CIP_R) / n(), sum(CIP_R), n()),
      AMP = sprintf("%.2f%% (%d/%d)", 100 * sum(AMP_R) / n(), sum(AMP_R), n()),
      TET = sprintf("%.2f%% (%d/%d)", 100 * sum(TET_R) / n(), sum(TET_R), n())
    )
  
  bind_rows(by_farm, total_row)
}

# ------------------------------------------------------------------
# 3. Build the two tables
# ------------------------------------------------------------------
faecalis_by_farm <- build_farm_resistance_table("faecalis")
faecium_by_farm  <- build_farm_resistance_table("faecium")

print(faecalis_by_farm)
print(faecium_by_farm)

write.csv(faecalis_by_farm, "faecalis_resistance_by_farm.csv", row.names = FALSE)
write.csv(faecium_by_farm, "faecium_resistance_by_farm.csv", row.names = FALSE)

# ------------------------------------------------------------------
# 4. LaTeX output (booktabs style, matching your existing table)
# ------------------------------------------------------------------
kable(faecalis_by_farm, format = "latex", booktabs = TRUE,
      caption = "Antibiotic resistance in \\textit{E. faecalis} isolates by farm",
      col.names = c("Farm", "Samples", "ERY", "TYL", "VAN", "CIP", "AMP", "TET"),
      align = "lccccccc", escape = FALSE)

kable(faecium_by_farm, format = "latex", booktabs = TRUE,
      caption = "Antibiotic resistance in \\textit{E. faecium} isolates by farm",
      col.names = c("Farm", "Samples", "ERY", "TYL", "VAN", "CIP", "AMP", "TET"),
      align = "lccccccc", escape = FALSE)




# ------------------------------------------------------------------
# 5. Milk samples: which isolates, from which farm, at what timepoint,
#    and their diffusion zones
#    ---> requires the sample metadata file (CalfID, Farm, Timepoint)
# ------------------------------------------------------------------
# 1. Rebuild amr from scratch (this brings milk rows back in)
amr <- data %>%
  filter(bact_genus == "Enterococcus",
         bact_species %in% c("faecalis", "faecium")) %>%
  mutate(
    ERY = ERYdiff, TYL = TYLdiff, VAN = VANdiff,
    CIP = CIPdiff, AMP = AMPdiff, TET = TETdiff
  )

# 2. Join metadata (only once!)
amr <- amr %>%
  left_join(sample_meta, by = "FKSampleID")

# 3. Compute resistance flags
amr <- amr %>%
  mutate(
    ERY_R = ifelse(ERY <= 13.01, 1, 0),
    TYL_R = case_when(
      bact_species == "faecium"  & TYL <= 11.01 ~ 1,
      bact_species == "faecium"  & TYL >  11.01 ~ 0,
      bact_species == "faecalis" & TYL <= 17.01 ~ 1,
      bact_species == "faecalis" & TYL >  17.01 ~ 0,
      TRUE ~ NA_real_
    ),
    VAN_R = ifelse(VAN <= 14.01, 1, 0),
    CIP_R = ifelse(CIP <= 15.01, 1, 0),
    AMP_R = ifelse(AMP <= 16.01, 1, 0),
    TET_R = ifelse(TET <= 14.01, 1, 0)
  )

# 4. Pull milk isolates out NOW, with correct _result columns
milk_isolates <- amr %>%
  filter(str_detect(str_to_lower(CalfID), "milk")) %>%
  mutate(
    ERY_result = ifelse(ERY_R == 1, "Resistant", "Not resistant"),
    TYL_result = ifelse(TYL_R == 1, "Resistant", "Not resistant"),
    VAN_result = ifelse(VAN_R == 1, "Resistant", "Not resistant"),
    CIP_result = ifelse(CIP_R == 1, "Resistant", "Not resistant"),
    AMP_result = ifelse(AMP_R == 1, "Resistant", "Not resistant"),
    TET_result = ifelse(TET_R == 1, "Resistant", "Not resistant")
  ) %>%
  select(IsolateID, FKSampleID, CalfID, FarmID, SamplingTimepoint,
         bact_genus, bact_species,
         ERY_result, TYL_result, VAN_result, CIP_result, AMP_result, TET_result)

print(milk_isolates)

# 5. THEN exclude milk for your main summary
amr <- amr %>%
  filter(!str_detect(str_to_lower(CalfID), "milk"))

amr <- amr %>%
  filter(!str_detect(str_to_lower(CalfID), "milk"))

milk_isolates <- amr %>%
     filter(str_detect(str_to_lower(CalfID), "milk")) %>%
  mutate(
    ERY_result = ifelse(ERY_R == 1, "Resistant", "Not resistant"),
    TYL_result = ifelse(TYL_R == 1, "Resistant", "Not resistant"),
    VAN_result = ifelse(VAN_R == 1, "Resistant", "Not resistant"),
    CIP_result = ifelse(CIP_R == 1, "Resistant", "Not resistant"),
    AMP_result = ifelse(AMP_R == 1, "Resistant", "Not resistant"),
    TET_result = ifelse(TET_R == 1, "Resistant", "Not resistant")
  ) %>%
   select(IsolateID, FKSampleID, CalfID, FarmID, SamplingTimepoint,
          bact_genus, bact_species,
          ERY, TYL, VAN, CIP, AMP, TET)
#
 write.csv(milk_isolates, "milk_isolates_summary.csv", row.names = FALSE)
 print(milk_isolates)
# ---------------------------------------------------------------------
