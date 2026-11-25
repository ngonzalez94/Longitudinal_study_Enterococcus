# ============================================================
# MDR isolate attrition - Faecalis - stacked funnel bar chart
# Total MDR isolates -> Sequenced -> Analysed
# Author: N. Gonzalez | Chapter 3
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. STARTING MDR isolates (n = 31)
# ------------------------------------------------------------
profiles_mdr <- tribble(
  ~profile,          ~n,
  "TET+TYL+CIP",     25,
  "TET+TYL+CIP+ERY",  5,
  "TYL+CIP+VAN",      1
)
sum(profiles_mdr$n)  # sanity check = 31

isolates_mdr <- profiles_mdr %>%
  uncount(n) %>%
  mutate(id = row_number(),
         source = "MDR_batch")

# ------------------------------------------------------------
# 2. Isolates ADDED at sequencing by mistake (non-MDR, n = 2)
# ------------------------------------------------------------
isolates_added <- tibble(
  profile = rep("ERY+TYL+CIP (non-MDR)", 2),
  id      = max(isolates_mdr$id) + seq_len(2),
  source  = "added_in_error"
)

# ------------------------------------------------------------
# 3. Isolates REMOVED before/at sequencing (milk isolates, TET+TYL+CIP, n = 2)
#    EDIT: replace with the actual isolate ids if you know which 2
# ------------------------------------------------------------
set.seed(31)
removed_ids <- isolates_mdr %>%
  filter(profile == "TET+TYL+CIP") %>%
  slice_sample(n = 2) %>%
  pull(id)

# ------------------------------------------------------------
# 4. Build stage membership flags
# ------------------------------------------------------------
all_isolates <- bind_rows(isolates_mdr, isolates_added)

all_isolates <- all_isolates %>%
  mutate(
    in_MDR       = source == "MDR_batch",                       # only original 31 count as MDR total
    in_Sequenced = source == "added_in_error" | !(id %in% removed_ids)  # add 2, drop the 2 milk isolates
  )

# ------------------------------------------------------------
# 5. Analysed stage: drop 2 (no reads) + 6 (bad quality) = 8 lost
#    from whoever made it into Sequenced
#    EDIT: replace sample() with explicit ids if known
# ------------------------------------------------------------
sequenced_ids <- all_isolates %>% filter(in_Sequenced) %>% pull(id)

set.seed(310)
lost_noreads   <- sample(sequenced_ids, 2)
remaining      <- setdiff(sequenced_ids, lost_noreads)
lost_badqual   <- sample(remaining, 6)

all_isolates <- all_isolates %>%
  mutate(
    in_Analysed = in_Sequenced & !id %in% c(lost_noreads, lost_badqual)
  )

# ------------------------------------------------------------
# 6. Long-format counts per stage x profile
# ------------------------------------------------------------
stage_counts <- bind_rows(
  all_isolates %>% filter(in_MDR)       %>% count(profile) %>% mutate(stage = "Total MDR\nisolates"),
  all_isolates %>% filter(in_Sequenced) %>% count(profile) %>% mutate(stage = "Sequenced"),
  all_isolates %>% filter(in_Analysed)  %>% count(profile) %>% mutate(stage = "Analysed")
) %>%
  mutate(stage = factor(stage, levels = c("Total MDR\nisolates", "Sequenced", "Analysed")))

stage_totals <- stage_counts %>%
  group_by(stage) %>%
  summarise(total = sum(n), .groups = "drop")

profile_order <- c(profiles_mdr$profile, "(non-MDR)")
stage_counts$profile <- factor(stage_counts$profile, levels = profile_order)

# ------------------------------------------------------------
# 7. Loss/gain annotation labels between bars
# ------------------------------------------------------------
loss_labels <- tibble(
  x     = c(1.5, 2.5),
  y     = c(max(stage_totals$total) * 1.1, max(stage_totals$total) * 1.1),
  label = c("+2 added in error\n-2 milk isolates removed",
            "-2 no reads\n-6 bad quality")
)

# ------------------------------------------------------------
# 8. Colours
# ------------------------------------------------------------
fill_cols <- setNames(scales::hue_pal()(length(profile_order)), profile_order)
fill_cols["ERY+TYL+CIP (non-MDR)"] <- "grey60"

# ------------------------------------------------------------
# 9. Plot
# ------------------------------------------------------------
p <- ggplot(stage_counts, aes(x = stage, y = n, fill = profile)) +
  geom_col(width = 0.55, colour = "grey20", linewidth = 0.2) +
  geom_text(data = stage_totals, aes(x = stage, y = total, label = paste0("n=", total)),
            inherit.aes = FALSE, vjust = -0.6, fontface = "bold", size = 4) +
  annotate("segment", x = 1.28, xend = 1.72, y = -2, yend = -2,
           arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  annotate("segment", x = 2.28, xend = 2.72, y = -2, yend = -2,
           arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_text(data = loss_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3, colour = "grey30", fontface = "italic") +
  scale_fill_manual(values = fill_cols, name = "Resistance profile") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, y = "Number of isolates",
       title = "Isolate attrition through sequencing and analysis (Species B)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  )

p

ggsave("MDR_attrition_funnel_speciesB.pdf", p, width = 7.5, height = 5.5, units = "in")
ggsave("MDR_attrition_funnel_speciesB.png", p, width = 7.5, height = 5.5, units = "in", dpi = 300)