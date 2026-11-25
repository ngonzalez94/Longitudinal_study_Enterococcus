if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggmsa")
library(DECIPHER)
library(ggplot2)
library(patchwork)
library(readxl)
library(ggdendro)
install.packages(ggalt)
library(ggalt)
library(ggmsa)
library(ggseqlogo)
library(Biostrings)
setwd("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/genes comparisson/plot")
# ── Data setup ────────────────────────────────────────────────────────────────

meta_all <- read_excel("metadata.xlsx")
meta_all$ID_short <- sub("\\|.*", "", meta_all$ID)
meta_all$ST <- as.character(meta_all$ST)

fasta_files <- list(
  "erm(B)" = "combined_ermB_aligned_1.fas",
  "ant(6)" = "combined_ant6_aligned.fas",
  "aph(3)" = "combined_aph3_aligned.fas",
  "tet(M)" = "combined_tetM_aligned_1.fas",
  "tet(S)" = "combined_tetS_aligned.fas"
)

gene_names <- names(fasta_files)

# Shared colour scales
col_species <- c("E. faecalis" = "#08d1bd", "E. faecium" = "#ff9a63")
col_org     <- c("Plasmid" = "lightblue",     "Chromosome" = "darkred")


# ── Helper: load + prep one gene ─────────────────────────────────────────────

prep_gene <- function(gene_name, fasta_path, meta_all) {
  
  seqs <- readAAStringSet(fasta_path)
  names(seqs) <- sub("\\|.*", "", names(seqs))
  
  dist_matrix     <- DistanceMatrix(seqs, type = "matrix", includeTerminalGaps = FALSE)
  identity_matrix <- (1 - dist_matrix) * 100
  
  hc        <- hclust(as.dist(dist_matrix), method = "complete")
  seq_order <- hc$labels[hc$order]
  
  # Metadata for this gene — includes ST and Diet Stage
  meta_gene <- meta_all[!is.na(meta_all[[gene_name]]) & meta_all[[gene_name]] != "",
                        c("ID_short", "Species", "ST", "Diet Stage", gene_name)]
  colnames(meta_gene)[5] <- "Organisation"  # ← position 5 now, not 4
  meta_gene$Organisation <- ifelse(toupper(meta_gene$Organisation) == "P",
                                   "Plasmid", "Chromosome")
  meta_gene$ID_short <- factor(meta_gene$ID_short, levels = seq_order)
  
  list(seqs = seqs, dist = dist_matrix, identity = identity_matrix,
       hc = hc, seq_order = seq_order, meta = meta_gene)
}

for (g in gene_names) {
  meta_gene <- meta_all[!is.na(meta_all[[g]]) & meta_all[[g]] != "", "ID_short"]
  dups <- meta_gene$ID_short[duplicated(meta_gene$ID_short)]
  if (length(dups) > 0) cat(g, "duplicates:", dups, "\n")
}

for (g in gene_names) {
  seqs <- readAAStringSet(fasta_files[[g]])
  names(seqs) <- sub("\\|.*", "", names(seqs))
  dups <- names(seqs)[duplicated(names(seqs))]
  if (length(dups) > 0) cat(g, "duplicates:", dups, "\n")
}
# Precompute all genes
gene_data <- mapply(prep_gene,
                    gene_name = gene_names,
                    fasta_path = fasta_files,
                    MoreArgs = list(meta_all = meta_all),
                    SIMPLIFY = FALSE)

# ── Plot functions ────────────────────────────────────────────────────────────

# 1. Dendrogram + metadata strips
make_dendrogram <- function(gd, gene_name) {
  
  dendro <- dendro_data(as.dendrogram(gd$hc))
  meta   <- gd$meta
  
  p_dendro <- ggplot() +
    geom_segment(data = segment(dendro),
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.4, colour = "grey30") +
    scale_x_continuous(breaks = seq_along(gd$seq_order),
                       labels = gd$seq_order) +
    theme_minimal() +
    theme(axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.title.x       = element_blank(),
          axis.title.y       = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title         = element_text(size = 11, hjust = 0.5, face = "italic"),
          legend.position    = "none") +
    labs(title = gene_name, y = "Distance")
  
  p_species <- ggplot(meta, aes(x = ID_short, y = 1, fill = Species)) +
    geom_tile() +
    scale_fill_manual(values = col_species, name = "Species") +
    theme_void() + theme(legend.position = "none")
  
  p_org <- ggplot(meta, aes(x = ID_short, y = 1, fill = Organisation)) +
    geom_tile() +
    scale_fill_manual(values = col_org, name = "Organisation") +
    theme_void() + theme(legend.position = "none")
  
  p_st <- ggplot(meta, aes(x = ID_short, y = 1, fill = factor(ST))) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3", name = "ST") +
    theme_void() + theme(legend.position = "none")
  
  p_diet <- ggplot(meta, aes(x = ID_short, y = 1, fill = `Diet Stage`)) +
    geom_tile() +
    scale_fill_manual(values = col_diet, name = "Diet Stage") +
    theme_void() + theme(legend.position = "none")
  
  p_dendro / p_species / p_org / p_st / p_diet +
    plot_layout(heights = c(10, 0.4, 0.4, 0.4, 0.4))
} 

# Shared colour scales
col_species <- c("E. faecalis" = "#08d1bd", "E. faecium" = "#ff9a63")
col_org     <- c("Plasmid" = "lightblue",   "Chromosome" = "darkred")
col_diet    <- c("Pre-WM" = "#f7c59f",      "WM" = "#e05c5c")

dendro_plots <- mapply(make_dendrogram, gene_data, gene_names, SIMPLIFY = FALSE)

library(patchwork)
library(grid)

# ── Build individual legends (same as before, but slightly tighter) ─────────

dummy_species <- ggplot(data.frame(Species = c("E. faecalis", "E. faecium"),
                                   x = 1:2, y = 1:2),
                        aes(x, y, colour = Species)) +
  geom_point(size = 4) +
  scale_colour_manual(values = col_species, name = "Species") +
  theme_minimal() +
  theme(legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

dummy_org <- ggplot(data.frame(Organisation = c("Plasmid", "Chromosome"),
                               x = 1:2, y = 1:2),
                    aes(x, y, fill = Organisation)) +
  geom_tile() +
  scale_fill_manual(values = col_org, name = "Organisation") +
  theme_minimal() +
  theme(legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

dummy_diet <- ggplot(data.frame(Diet = c("Pre-WM", "WM"),
                                x = 1:2, y = 1:2),
                     aes(x, y, fill = Diet)) +
  geom_tile() +
  scale_fill_manual(values = col_diet, name = "Diet Stage") +
  theme_minimal() +
  theme(legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

st_levels <- sort(unique(meta_all$ST))

dummy_st <- ggplot(data.frame(ST = factor(st_levels),
                              x = seq_along(st_levels),
                              y = 1),
                   aes(x, y, fill = ST)) +
  geom_tile() +
  scale_fill_brewer(palette = "Set3", name = "ST") +
  theme_minimal() +
  theme(legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"))

library(scales)
st_levels <- sort(unique(meta_all$ST))
n_st      <- length(st_levels)
st_palette <- hue_pal()(n_st)          # ggplot's default hue wheel, evenly spaced
names(st_palette) <- st_levels


leg_species <- cowplot::get_legend(dummy_species)
leg_org     <- cowplot::get_legend(dummy_org)
leg_diet    <- cowplot::get_legend(dummy_diet)
leg_st      <- cowplot::get_legend(dummy_st)

# ── Stack with patchwork instead of cowplot::plot_grid ───────────────────────
# Heights proportional to (number of entries + 1 for the title row)

n_species <- 2 + 1
n_org     <- 2 + 1
n_diet    <- 2 + 1
n_st      <- length(st_levels) + 1

legend_gg <- wrap_elements(full = leg_species) /
  wrap_elements(full = leg_org)     /
  wrap_elements(full = leg_diet)    /
  wrap_elements(full = leg_st) +
  plot_layout(heights = c(n_species, n_org, n_diet, n_st))

legend_theme <- theme(
  legend.title       = element_text(size = 8, face = "bold"),
  legend.text        = element_text(size = 7),
  legend.key.size    = unit(0.3, "cm"),
  legend.spacing.y   = unit(0.05, "cm"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.margin  = margin(0, 0, 0, 0)
)

dummy_species <- dummy_species + legend_theme
dummy_org     <- dummy_org     + legend_theme
dummy_diet    <- dummy_diet    + legend_theme
dummy_st      <- dummy_st      + legend_theme


# ── Assemble final figure — legend gets its OWN full-height column ──────────

fig1_dendros <- (dendro_plots[["erm(B)"]] | dendro_plots[["ant(6)"]] | dendro_plots[["aph(3)"]]) /
  (dendro_plots[["tet(M)"]] | dendro_plots[["tet(S)"]] | plot_spacer())

fig1 <- fig1_dendros | legend_gg +
  plot_layout(widths = c(3, 3, 3, 0.6)) +
  plot_annotation(title = "Hierarchical clustering dendrograms",
                  tag_levels = "A")

ggsave("Fig1_dendrograms_complete.pdf", fig1, width = 20, height = 12)




# 2. Sequence logo
make_logo <- function(gd, gene_name) {
  
  # Convert to data frame for ggmsa
  tidy <- tidy_msa(gd$seqs)
  
  ggplot() +
    geom_msa(data = tidy,
             seq_name = FALSE,
             color = "Chemistry_AA") +
    labs(title = gene_name,
         x = "Position",
         y = NULL) +
    theme_minimal() +
    theme(plot.title  = element_text(size = 10, hjust = 0.5, face = "italic"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

# 3. Identity heatmap + metadata strips
make_heatmap <- function(gd, gene_name) {
  
  mat_df <- as.data.frame(as.table(gd$identity))
  colnames(mat_df) <- c("SeqA", "SeqB", "Identity")
  mat_df$SeqA <- factor(mat_df$SeqA, levels = gd$seq_order)
  mat_df$SeqB <- factor(mat_df$SeqB, levels = gd$seq_order)
  
  meta <- gd$meta
  
  p_main <- ggplot(mat_df, aes(x = SeqA, y = SeqB, fill = Identity)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkblue",
                        name = "% Identity", limits = c(0, 100)) +
    theme_minimal() +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5, face = "italic")) +
    labs(title = gene_name)
  
  p_species <- ggplot(meta, aes(x = ID_short, y = 1, fill = Species)) +
    geom_tile() +
    scale_fill_manual(values = col_species, name = "Species") +
    theme_void() + theme(legend.position = "right")
  
  p_org <- ggplot(meta, aes(x = ID_short, y = 1, fill = Organisation)) +
    geom_tile() +
    scale_fill_manual(values = col_org, name = "Organisation") +
    theme_void() + theme(legend.position = "right")
  
  p_st <- ggplot(meta, aes(x = ID_short, y = 1, fill = factor(ST))) +
    geom_tile() +
    scale_fill_brewer(palette = "Set3", name = "ST") +
    theme_void() + theme(legend.position = "right")
  
  p_main / p_species / p_org / p_st +
    plot_layout(heights = c(10, 0.4, 0.4, 0.4))
}

# 4. MSA (one per gene, saved individually)
make_msa <- function(gd, gene_name) {
  ggmsa(gd$seqs,
        color = "Chemistry_AA",
        font = "DroidSansMono",
        char_width = 0.5,
        seq_name = TRUE) +
    labs(title = gene_name) +
    theme(plot.title = element_text(size = 10, hjust = 0.5, face = "italic"))
}

# ── Build multipanel figures ──────────────────────────────────────────────────

# Helper to arrange 5 plots in 2x3 grid (last cell empty)
arrange_2x3 <- function(plot_list, gene_names) {
  p <- setNames(plot_list, gene_names)
  (p[[1]] | p[[2]] | p[[3]]) /
    (p[[4]] | p[[5]] | plot_spacer())
}

# Figure 1 — Dendrograms
dendro_plots <- mapply(make_dendrogram, gene_data, gene_names, SIMPLIFY = FALSE)

fig1 <- (dendro_plots[["erm(B)"]] | dendro_plots[["ant(6)"]] | dendro_plots[["aph(3)"]]) /
  (dendro_plots[["tet(M)"]] | dendro_plots[["tet(S)"]]  | legend_gg) +
  plot_annotation(title = "Hierarchical clustering dendrograms",
                  tag_levels = "A")

# Figure 2 — Sequence logos
logo_plots <- mapply(make_logo, gene_data, gene_names, SIMPLIFY = FALSE)
fig2 <- arrange_2x3(logo_plots, gene_names) +
  plot_annotation(title = "Sequence logos",
                  tag_levels = "A")

# Figure 3 — Heatmaps
heatmap_plots <- mapply(make_heatmap, gene_data, gene_names, SIMPLIFY = FALSE)
fig3 <- arrange_2x3(heatmap_plots, gene_names) +
  plot_annotation(title = "Pairwise identity heatmaps",
                  tag_levels = "A")







# ── Save figures ──────────────────────────────────────────────────────────────

ggsave("Fig1_dendrograms_legend_2.pdf", fig1, width = 18, height = 10)
ggsave("Fig2_logos.pdf",       fig2, width = 18, height = 10)
ggsave("Fig3_heatmaps_metadata_2.pdf",    fig3, width = 18, height = 10)

# Figures 4-8 — MSA one per gene
for (g in gene_names) {
  p_msa <- make_msa(gene_data[[g]], g)
  safe_name <- gsub("[()']", "", g)  # remove special chars for filename
  ggsave(paste0("Fig_MSA_", safe_name, ".pdf"),
         p_msa, width = 20, height = nrow(as.matrix(gene_data[[g]]$seqs)) * 0.4 + 2)
}

cat("All figures saved.\n")



##INDIVIDUAL DENDOGRAMS

# Load sequences

for (g in gene_names) {
  
  # Load sequences
  seqs <- readAAStringSet(fasta_files[[g]])
  names(seqs) <- sub("\\|.*", "", names(seqs))
  
  # Distance and clustering
  dist_matrix <- DistanceMatrix(seqs, type = "matrix", includeTerminalGaps = FALSE)
  hc          <- hclust(as.dist(dist_matrix), method = "complete")
  
  # Dendrogram data
  dendro    <- dendro_data(as.dendrogram(hc))
  labels_df <- dendro$labels
  
  # Add metadata — now including Diet Stage
  meta_gene <- meta_all[!is.na(meta_all[[g]]) & meta_all[[g]] != "",
                        c("ID_short", "Species", "Diet Stage", g)]
  colnames(meta_gene)[4] <- "Organisation"
  meta_gene$Organisation <- ifelse(toupper(meta_gene$Organisation) == "P",
                                   "Plasmid", "Chromosome")
  
  labels_df <- merge(labels_df, meta_gene, by.x = "label", by.y = "ID_short", all.x = TRUE)
  
  # Sequence order for strips
  seq_order <- hc$labels[hc$order]
  meta_gene$ID_short <- factor(meta_gene$ID_short, levels = seq_order)
  
  # Main dendrogram
  p_dendro <- ggplot() +
    geom_segment(data = segment(dendro),
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.4, colour = "grey30") +
    geom_point(data = labels_df,
               aes(x = x, y = 0, colour = Species, shape = Organisation),
               size = 3) +
    geom_text(data = labels_df,
              aes(x = x, y = -0.01, label = label, colour = Species),
              angle = 90, hjust = 1, size = 2.5) +
    scale_colour_manual(values = col_species, name = "Species") +
    scale_shape_manual(values = c("Plasmid" = 16, "Chromosome" = 16),
                       name = "Organisation") +
    scale_y_continuous(expand = expansion(mult = c(0.3, 0.05))) +
    theme_minimal() +
    theme(axis.text.x        = element_blank(),
          axis.ticks.x       = element_blank(),
          axis.title.x       = element_blank(),
          axis.title.y       = element_text(size = 9),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title         = element_text(size = 11, hjust = 0.5, face = "italic"),
          legend.position    = "right") +
    labs(title = g, y = "Distance")
  
  # Metadata strips below dendrogram
  p_species <- ggplot(meta_gene, aes(x = ID_short, y = 1, fill = Species)) +
    geom_tile() +
    scale_fill_manual(values = col_species, name = "Species") +
    theme_void() +
    theme(legend.position = "none")
  
  p_org <- ggplot(meta_gene, aes(x = ID_short, y = 1, fill = Organisation)) +
    geom_tile() +
    scale_fill_manual(values = col_org, name = "Organisation") +
    theme_void() +
    theme(legend.position = "none")
  
  p_diet <- ggplot(meta_gene, aes(x = ID_short, y = 1, fill = `Diet Stage`)) +
    geom_tile() +
    scale_fill_manual(values = col_diet, name = "Diet Stage") +
    theme_void() +
    theme(legend.position = "none")
  
  # Stack dendrogram + strips
  p_final <- p_dendro / p_species / p_org / p_diet +
    plot_layout(heights = c(10, 0.4, 0.4, 0.4))
  
  # Save
  safe_name <- gsub("[()']", "", g)
  n_seqs    <- length(seqs)
  ggsave(paste0("Dendrogram_", safe_name, ".pdf"),
         p_final,
         width  = max(8, n_seqs * 0.3),
         height = 8)
  
  cat("Saved dendrogram for", g, "\n")
}
