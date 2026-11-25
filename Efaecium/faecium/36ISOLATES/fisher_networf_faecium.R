# =============================================================================
# AMR Gene Co-occurrence Network — Fisher's Exact Test + BH Correction
# Input : AMRFinderPlus results (all isolates in one file)
# Output: network plot (.png / .pdf) + full results table (.csv)
# =============================================================================

# install.packages(c("tidyverse", "igraph", "ggraph", "ggrepel", "readxl"))

library(tidyverse)
library(tibble)
library(igraph)
library(ggraph)
library(ggrepel)
library(readxl)   

# =============================================================================
# 1. LOAD DATA
# =============================================================================
# Adjust the path and format to your file:
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/36ISOLATES/")
df <- read_tsv("merged_output.tsv", locale = locale(encoding = "UTF-8"))   # tab-separated
#df <- read_excel("amrfinder_combined_2.xlsx", sheet = 1)  # Excel

# Rename to safe column names (adjust if your headers differ)
df <- dplyr::rename(df,
                    sample    = Isolate,
                    gene      = `Element symbol`,
                    amr_class = Class,
                    type      = Type
)

# =============================================================================
# 2. CLEAN & PREPARE
# =============================================================================
df <- df %>%
  # Keep only AMR entries with a named gene (drop STRESS, POINT mutations
  # without a gene symbol, etc. — adjust the filter to suit your needs)
  filter(type == "AMR", !is.na(gene)) %>%
  # Multi-class genes (e.g. LINCOSAMIDE/MACROLIDE) — keep first class only.
  # Change str_split_i(amr_class, "/", 1) to keep a different token if needed.
  mutate(amr_class = str_split_i(amr_class, "/", 1)) %>%
  distinct(sample, gene, amr_class)

n_isolates <- n_distinct(df$sample)
cat("Isolates:", n_isolates, "| Unique AMR genes:", n_distinct(df$gene), "\n")

# =============================================================================
# 3. BUILD PRESENCE / ABSENCE MATRIX
# =============================================================================
pa_matrix <- df %>%
  distinct(sample, gene) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0L) %>%
  column_to_rownames("sample")

mat   <- as.matrix(pa_matrix)
genes <- colnames(mat)

# =============================================================================
# 4. FISHER'S EXACT TEST FOR EVERY GENE PAIR
# =============================================================================
# For each pair we build a 2x2 contingency table:
#
#              Gene B present  Gene B absent
#  Gene A present      a              b
#  Gene A absent       c              d
#
# alternative = "greater" → one-sided test: are genes co-occurring MORE
# than expected by chance? (odds ratio > 1)

run_fisher <- function(g1, g2, mat) {
  a <- sum(mat[, g1] == 1 & mat[, g2] == 1)   # both present
  b <- sum(mat[, g1] == 1 & mat[, g2] == 0)   # g1 only
  c <- sum(mat[, g1] == 0 & mat[, g2] == 1)   # g2 only
  d <- sum(mat[, g1] == 0 & mat[, g2] == 0)   # neither
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  data.frame(
    gene1      = g1,
    gene2      = g2,
    both       = a,       # isolates carrying both
    only_g1    = b,
    only_g2    = c,
    neither    = d,
    p_value    = ft$p.value,
    odds_ratio = as.numeric(ft$estimate)
  )
}

gene_pairs    <- combn(genes, 2, simplify = FALSE)
fisher_results <- map_dfr(gene_pairs, ~run_fisher(.x[1], .x[2], mat))

# =============================================================================
# 5. BENJAMINI-HOCHBERG CORRECTION FOR MULTIPLE TESTING
# =============================================================================
fisher_results <- fisher_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

cat("Total pairs tested    :", nrow(fisher_results), "\n")
cat("Significant (p_adj<0.05):", sum(fisher_results$p_adj < 0.05), "\n\n")

# Save full results table
write_csv(fisher_results, "amr_fisher_results.csv")
cat("Full results saved to amr_fisher_results.csv\n\n")

# Print significant pairs
fisher_results %>%
  filter(p_adj < 0.05) %>%
  select(gene1, gene2, both, p_value, p_adj, odds_ratio) %>%
  print(n = 50)

# =============================================================================
# 6. FILTER TO SIGNIFICANT POSITIVE CO-OCCURRENCES ONLY
# =============================================================================
sig_pairs <- fisher_results %>%
  filter(p_adj < 0.05, odds_ratio > 1)   # positive co-occurrence only

if (nrow(sig_pairs) == 0) {
  stop("No significant pairs — try relaxing the threshold or check isolate count.")
}

# =============================================================================
# 7. GENE METADATA (frequency + class)
# =============================================================================
gene_meta <- df %>%
  group_by(gene, amr_class) %>%
  summarise(frequency = n_distinct(sample), .groups = "drop")

# =============================================================================
# 8. BUILD IGRAPH OBJECT
# =============================================================================
g_net <- graph_from_data_frame(
  sig_pairs %>% select(gene1, gene2, both, p_adj, odds_ratio),
  directed = FALSE,
  vertices = gene_meta %>% filter(gene %in% c(sig_pairs$gene1, sig_pairs$gene2))
)

# =============================================================================
# 9. COLOUR PALETTE — EDIT TO MATCH YOUR OTHER PLOTS
# =============================================================================
# Run unique(gene_meta$amr_class) to see all classes in your data,
# then add/change hex values below as needed.

class_colours <- c(
  "TETRACYCLINE"   = "#E6814A",
  "AMINOGLYCOSIDE" = "#4A90C4",
  "PLEUROMUTILIN"   = "#50B87A",
  "LINCOSAMIDE"    = "#9B59B6",
  "LINCOSAMIDE/STREPTOGRAMIN"    = "#E8435A",
  "PHENICOL"       = "#7B8FA6",
  "FOSFOMYCIN"     = "#BA7517",
  "COLISTIN"       = "#D4537E",
  "NITROFURAN"     = "#95A5A6",
  "FOSMIDOMYCIN"   = "#BDC3C7",
  "QUINOLONE"      = "#C0392B",
  "LIPOPEPTIDE"   = "#6C3483",
  "RIFAMYCIN"      = "#A04000",
  "TRIMETHOPRIM"   = "#1E8449"
)

# =============================================================================
# 10. NETWORK PLOT
# =============================================================================
set.seed(42)   # change seed if you want a different layout arrangement

p <- ggraph(g_net, layout = "fr") +   # try "stress", "kk", "nicely" for alternatives
  
  # Edges: width = number of co-carrying isolates; opacity = statistical strength
  geom_edge_link(
    aes(width = both, alpha = -log10(p_adj)),
    colour      = "grey50",
    show.legend = FALSE
  ) +
  scale_edge_width(range = c(0.8, 6)) +
  scale_edge_alpha(range = c(0.35, 0.95)) +
  
  # Nodes: size = how many isolates carry the gene; colour = AMR class
  geom_node_point(
    aes(size = frequency, colour = amr_class),
    alpha = 0.92
  ) +
  scale_size_continuous(name = "Frequency\n(isolates)", range = c(4, 14)) +
  scale_colour_manual(name = "AMR class", values = class_colours, na.value = "grey50") +
  
  # Gene name labels with auto-repel to avoid overlap
  geom_node_text(
    aes(label = name),
    repel        = TRUE,
    size         = 5,
    colour       = "grey10",
    bg.colour    = "white",
    bg.r         = 0.12,
    max.overlaps = 30
  ) +
  
  theme_graph(base_family = "sans") +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 13, face = "bold"),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 8,  colour = "grey40", hjust = 0.5),
    plot.margin     = margin(10, 10, 10, 10)
  ) +
  labs(
    title    = "AMR Gene Co-occurrence Network",
    subtitle = paste0(
      "Fisher's exact test (one-sided), BH-adjusted p < 0.05  |  n = ",
      n_isolates,
      " isolates  |  Edge width = shared isolate count  |  Edge opacity = -log\u2081\u2080(p_adj)"
    )
  )

# =============================================================================
# 11. SAVE
# =============================================================================
ggsave("amr_network_fisher.pdf", plot = p, width = 12, height = 8,
       device = cairo_pdf)   # use "pdf" if cairo not available
ggsave("amr_network_fisher.png", plot = p, width = 12, height = 8, dpi = 300)

cat("Saved: amr_network_fisher.pdf / .png\n")
