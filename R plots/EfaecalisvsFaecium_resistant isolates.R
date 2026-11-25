library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── 1. Load data ──────────────────────────────────────────────────────────────
# Update the path and sheet names to match your file
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus")
faecalis <- read_excel("Ecocc_temp.xlsx", sheet = 1)
faecium  <- read_excel("Ecocc_temp.xlsx", sheet = 3)

# ── 2. Classify resistance (1 = resistant, 0 = susceptible) ──────────────────
classify_resistance <- function(df, tyl_cutoff) {
  df %>%
    mutate(
      AMP_R = ifelse(AMP <= 16.01, 1, 0),
      CIP_R = ifelse(CIP <= 15.01, 1, 0),
      ERY_R = ifelse(ERY <= 13.01, 1, 0),
      TET_R = ifelse(TET <= 14.01, 1, 0),
      TYL_R = ifelse(TYL <= tyl_cutoff, 1, 0),
      VAN_R = ifelse(VAN <= 14.01, 1, 0)
    )
}

faecalis <- classify_resistance(faecalis, tyl_cutoff = 17.01)
faecium  <- classify_resistance(faecium,  tyl_cutoff = 11.01)
# ── 3. Summarise: % resistant, n resistant, n total per antimicrobial ─────────
summarise_resistance <- function(df, species_label) {
  drugs <- c("AMP_R", "CIP_R", "ERY_R", "TET_R", "TYL_R", "VAN_R")
  
  data.frame(
    Species       = species_label,
    Antimicrobial = sub("_R$", "", drugs),
    pct           = sapply(drugs, function(d) mean(df[[d]], na.rm = TRUE) * 100),
    n_res         = sapply(drugs, function(d) sum(df[[d]], na.rm = TRUE)),
    n_tot         = sapply(drugs, function(d) sum(!is.na(df[[d]]))),
    row.names     = NULL
  ) %>%
    mutate(
      label_top = sprintf("%.1f%%", pct),
      label_bot = sprintf("(n=%d/%d)", n_res, n_tot)
    )
}

summary_df <- bind_rows(
  summarise_resistance(faecalis, "E. faecalis"),
  summarise_resistance(faecium,  "E. faecium")
)

# Fix antimicrobial order
summary_df$Antimicrobial <- factor(summary_df$Antimicrobial,
                                   levels = c("AMP", "CIP", "ERY", "TET", "TYL", "VAN"))

# ── 4. Plot ───────────────────────────────────────────────────────────────────
p <- ggplot(summary_df, aes(x = Antimicrobial, y = pct, fill = Species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, colour = "white") +
  # percentage label just above the bar
  geom_text(
    aes(label = label_top),
    position = position_dodge(width = 0.7),
    vjust = -3, size = 4.5, colour = "grey20"
  ) +
  # (n=x/total) label just below the percentage
  geom_text(
    aes(label = label_bot),
    position = position_dodge(width = 0.7),
    vjust = -1.6, size = 3.5, colour = "grey40"
  ) +
  scale_fill_manual(values = c("E. faecalis" = "#2196A6", "E. faecium" = "#E07B54")) +
  scale_y_continuous(
    limits = c(0, 115),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(
    #title    = "Antimicrobial Resistance by Species",
   # subtitle = "TYL cut-off: E. faecalis ≤ 17 mm, E. faecium ≤ 11 mm",
    x        = "Antimicrobial",
    y        = "% Resistant isolates",
    fill     = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(face = "italic"),
    plot.title        = element_text(face = "bold"),
    axis.text         = element_text(colour = "black", size = 12),
    panel.grid.major.y = element_line(colour = "grey90")
  )

print(p)

# ── 5. Save ───────────────────────────────────────────────────────────────────
ggsave("resistance_plot.png", plot = p, width = 8, height = 5, dpi = 300)
message("Plot saved to resistance_plot.png")