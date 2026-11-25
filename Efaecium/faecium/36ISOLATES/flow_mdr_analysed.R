# ============================================================
# MDR isolate attrition Sankey/alluvial plot
# Total MDR isolates -> Sequenced -> Analysed
# Author: N. Gonzalez | Chapter 3
# ============================================================

library(tidyverse)
library(ggalluvial)

# ------------------------------------------------------------
# 1. EDIT THIS: one row per isolate, with its resistance profile
#    and the stage it reached (NA / "Excluded_..." if dropped out)
# ------------------------------------------------------------
# Build isolate-level table. Replace n= values with your real counts.
# profile = resistance phenotype string (as in your sketch)

profiles <- tribble(
  ~profile,            ~n,
  "CIP+ERY+TET+TYL",    37,
  "CIP+TET+TYL",        1,
  "CIP+ERY+TET",       1,  
  "CIP+TET+AMP",        3,
  "CIP+ERY+AMP",        1
)
# should sum to your starting N (43 here) - check:
sum(profiles$n)

# Expand to isolate-level rows
isolates <- profiles %>%
  uncount(n) %>%
  mutate(id = row_number())

# ------------------------------------------------------------
# 2. Assign exclusions at each stage
#    - lost_seq: isolates that failed to meet sequencing standard (lost between MDR -> Sequenced)
#    - lost_assembly: isolates that failed assembly (lost between Sequenced -> Analysed)
#    EDIT: put the actual isolate ids you excluded, or just set n_lost1 / n_lost2
#    and this will randomly (reproducibly) pick which ones drop out.
# ------------------------------------------------------------
set.seed(43)
n_total   <- nrow(isolates)
n_lost1   <- 3   # didn't meet sequencing standard
n_lost2   <- 3   # failed assembly (drawn from those that WERE sequenced)

lost1_ids <- sample(isolates$id, n_lost1)
remaining_after_seq <- setdiff(isolates$id, lost1_ids)
lost2_ids <- sample(remaining_after_seq, n_lost2)

isolates <- isolates %>%
  mutate(
    stage_MDR      = profile,   # everyone starts here with their profile
    stage_Sequenced = case_when(
      id %in% lost1_ids ~ "Excluded (below standard)",
      TRUE               ~ profile
    ),
    stage_Analysed = case_when(
      id %in% lost1_ids ~ "Excluded (below standard)",
      id %in% lost2_ids ~ "Excluded (no assembly)",
      TRUE               ~ profile
    )
  )

# ------------------------------------------------------------
# 3. Convert to lodes (long) format for ggalluvial
# ------------------------------------------------------------
alluvial_data <- isolates %>%
  select(id, profile, stage_MDR, stage_Sequenced, stage_Analysed) %>%
  to_lodes_form(axes = 2:4, id = "id") %>%
  rename(stage = x, value = stratum)

# order stages left to right
alluvial_data$stage <- factor(
  alluvial_data$stage,
  levels = c("stage_MDR", "stage_Sequenced", "stage_Analysed"),
  labels = c("Total MDR\nisolates", "Sequenced", "Analysed")
)

# make sure "Excluded" categories sit at bottom of stack, greyed out
excl_levels <- unique(alluvial_data$value[grepl("Excluded", alluvial_data$value)])
profile_levels <- setdiff(unique(alluvial_data$value), excl_levels)
alluvial_data$value <- factor(alluvial_data$value,
                              levels = c(profile_levels, excl_levels))

# ------------------------------------------------------------
# 4. Colour palette: profiles get colour, exclusions get grey
# ------------------------------------------------------------
n_profiles <- length(profile_levels)
profile_cols <- setNames(
  scales::hue_pal()(n_profiles),
  profile_levels
)
excl_cols <- setNames(rep("grey80", length(excl_levels)), excl_levels)
fill_cols <- c(profile_cols, excl_cols)

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------
p <- ggplot(alluvial_data,
            aes(x = stage, stratum = value, alluvium = id,
                y = 1, fill = value, label = value)) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "grey40", linewidth = 0.15, alpha = 0.85) +
  geom_stratum(width = 0.35, color = "grey20") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.6, check_overlap = TRUE) +
  scale_fill_manual(values = fill_cols, guide = "none") +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  labs(x = NULL, y = NULL,
       title = "MDR E. faecalis isolate attrition through sequencing and analysis") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 11, face = "bold")
  )

# ------------------------------------------------------------
# 6. Annotate stage totals + losses along the top
# ------------------------------------------------------------
stage_totals <- alluvial_data %>%
  distinct(id, stage, value) %>%
  filter(!grepl("Excluded", value)) %>%
  count(stage, name = "n_remaining")

p <- p + annotate("text",
                  x = stage_totals$stage,
                  y = n_total + 2,
                  label = paste0("n = ", stage_totals$n_remaining),
                  fontface = "bold", size = 3.5)

p

ggsave("MDR_attrition_alluvial.pdf", p, width = 7, height = 5, units = "in")
ggsave("MDR_attrition_alluvial.png", p, width = 7, height = 5, units = "in", dpi = 300)



###LET'S TRY A FUNNEL BETTER: 
# ------------------------------------------------------------
# 2. Assign which isolates were lost at each stage
#    EDIT: replace sample() with explicit isolate ids if you know
#    exactly which ones failed sequencing / assembly
# ------------------------------------------------------------
set.seed(43)
n_lost1 <- 3   # didn't meet sequencing standard
n_lost2 <- 3   # failed assembly

lost1_ids <- sample(isolates$id, n_lost1)
remaining_after_seq <- setdiff(isolates$id, lost1_ids)
lost2_ids <- sample(remaining_after_seq, n_lost2)

isolates <- isolates %>%
  mutate(
    in_MDR        = TRUE,
    in_Sequenced  = !id %in% lost1_ids,
    in_Analysed   = !id %in% lost1_ids & !id %in% lost2_ids
  )

# ------------------------------------------------------------
# 3. Build long-format counts per stage x profile (only isolates present at that stage)
# ------------------------------------------------------------
stage_counts <- bind_rows(
  isolates %>% filter(in_MDR)       %>% count(profile) %>% mutate(stage = "Total MDR\nisolates"),
  isolates %>% filter(in_Sequenced) %>% count(profile) %>% mutate(stage = "Sequenced"),
  isolates %>% filter(in_Analysed)  %>% count(profile) %>% mutate(stage = "Analysed")
) %>%
  mutate(stage = factor(stage, levels = c("Total MDR\nisolates", "Sequenced", "Analysed")))

stage_totals <- stage_counts %>%
  group_by(stage) %>%
  summarise(total = sum(n), .groups = "drop")

# keep profile order consistent across bars (largest at bottom)
profile_order <- profiles %>% arrange(desc(n)) %>% pull(profile)
stage_counts$profile <- factor(stage_counts$profile, levels = profile_order)

# ------------------------------------------------------------
# 4. Loss annotation data (drawn between bars)
# ------------------------------------------------------------
loss_labels <- tibble(
  x     = c(1.5, 2.5),
  y     = c(max(stage_totals$total) * 1.05, max(stage_totals$total) * 1.05),
  label = c(paste0("- ", n_lost1, " (below\nsequencing standard)"),
            paste0("- ", n_lost2, " (failed\nassembly)"))
)

# ------------------------------------------------------------
# 5. Colour palette
# ------------------------------------------------------------
fill_cols <- setNames(scales::hue_pal()(length(profile_order)), profile_order)

# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------
p <- ggplot(stage_counts, aes(x = stage, y = n, fill = profile)) +
  geom_col(width = 0.55, colour = "grey20", linewidth = 0.2) +
  geom_text(data = stage_totals, aes(x = stage, y = total, label = paste0("n=", total)),
            inherit.aes = FALSE, vjust = -0.6, fontface = "bold", size = 4) +
  # funnel connector arrows between bars
    annotate("segment", x = 1.28, xend = 1.72, y = -2, yend = -2,
           arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  annotate("segment", x = 2.28, xend = 2.72, y = -2, yend = -2,
           arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_text(data = loss_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3, colour = "grey30", fontface = "italic") +
  scale_fill_manual(values = fill_cols, name = "Resistance profile") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(x = NULL, y = "Number of isolates",
       title = "MDR E. faecalis isolate attrition through sequencing and analysis") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

p
