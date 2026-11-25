# ============================================================
#  AMR Heatmap per Farm  –  E. coli resistance gene analysis
#
#  Layout:
#    Rows    = isolates (y-axis), clustered within Diet Stage groups
#    Columns = genes (x-axis), split into AMR GENES | POINT MUTATIONS
#              and colour-annotated by antimicrobial Class below
#
#  Requirements:
#    install.packages(c("tidyverse","readxl","circlize","gridExtra"))
#    BiocManager::install("ComplexHeatmap")
# ============================================================

library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)

# ============================================================
# 0.  FILE PATHS  –  adjust these two lines
# ============================================================
META_FILE <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/Efaecium_metadata.xlsx"
AMR_FILE  <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/amrfinder_faecium.xlsx"

# ============================================================
# 1.  COLOUR MAPS
# ============================================================
diet_stage_order <- c("Pre-WM", "WM", "Post-WM")

diet_colors <- c(
  "Pre-WM"  = "#fee8c8",
  "WM"      = "#fdbb84",
  "Post-WM" = "#e34a33"
)

rp_label_map <- c(
  "A" = "CIP-TET-TYL-ERY",     "B" = "CIP-TET-TYL-AMP"
)

custom_rp_colors <- c(
  "CIP-TET-TYL-ERY" = "#e41a1c",
  "CIP-TET-TYL-AMP" = "#377eb8"
)

custom_st_colors <- c(
  "771"   = "#a50026", "870"   = "#ffffbf", "22"   = "#006837", "NM"   = "grey"
)

# 12-colour palette for AMR classes
class_palette <- c(
  "#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3",
  "#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd",
  "#ccebc5","#ffed6f"
)

# ============================================================
# 2.  LOAD & CLEAN DATA
# ============================================================
meta <- read_excel(META_FILE) %>%
  mutate(
    Isolate    = str_trim(Isolate),
    ST         = as.character(ST),
    Diet_Stage = factor(Diet_Stage, levels = diet_stage_order, ordered = TRUE),
    RP_label   = rp_label_map[as.character(Resistance_Profile)]
  )

amr_raw <- read_excel(AMR_FILE) %>%
  mutate(
    Sample            = str_trim(Sample),
    `Element subtype` = str_trim(`Element subtype`)
  )

message("Total rows in AMR file: ", nrow(amr_raw))

amr <- amr_raw %>%
  filter(`Element subtype` %in% c("AMR", "POINT"))

message("Rows after filtering to AMR/POINT: ", nrow(amr))

# Isolates in metadata with no AMR/POINT entries are legitimately all-absent
# (they may only have STRESS/BIOCIDE entries, or no entries at all)
no_amr_isos <- setdiff(meta$Isolate, unique(amr$Sample))
if (length(no_amr_isos) > 0) {
  message("ℹ ", length(no_amr_isos), " isolate(s) have no AMR/POINT genes — ",
          "will appear as all-absent (grey) rows in the heatmap.")
}

# ── Diagnostic: report any mismatches between the two files ──
in_meta_not_amr <- setdiff(meta$Isolate, unique(amr$Sample))
in_amr_not_meta <- setdiff(unique(amr$Sample), meta$Isolate)

if (length(in_meta_not_amr) > 0) {
  message("⚠ These isolates are in the METADATA but NOT in the AMR file")
  message("  (they will appear as all-absent rows, not NA):")
  message("  ", paste(in_meta_not_amr, collapse = ", "))
}
if (length(in_amr_not_meta) > 0) {
  message("⚠ These samples are in the AMR file but NOT in the metadata")
  message("  (they will be ignored):")
  message("  ", paste(in_amr_not_meta, collapse = ", "))
}

# ============================================================
# 3.  GENE ORDER  &  PRESENCE/ABSENCE MATRIX
# ============================================================
amr_dedup <- amr %>%
  distinct(Sample, `Gene symbol`, .keep_all = TRUE)

gene_info <- amr_dedup %>%
  distinct(`Gene symbol`, Class, `Element subtype`) %>%
  rename(gene = `Gene symbol`, element_subtype = `Element subtype`)

# Genes on x-axis: AMR block (sorted by Class) then POINT block (sorted by Class)
amr_genes   <- gene_info %>% filter(element_subtype == "AMR")   %>% arrange(Class)
point_genes <- gene_info %>% filter(element_subtype == "POINT") %>% arrange(Class)
ordered_genes <- c(amr_genes$gene, point_genes$gene)

n_amr_genes   <- nrow(amr_genes)
n_point_genes <- nrow(point_genes)

# Presence/absence matrix  (isolates = rows, genes = columns)
pres_mat <- amr_dedup %>%
  select(Sample, `Gene symbol`) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = `Gene symbol`, values_from = present,
              values_fill = 0L) %>%
  column_to_rownames("Sample")

# Add any genes not observed in any isolate
missing_genes <- setdiff(ordered_genes, colnames(pres_mat))
if (length(missing_genes) > 0) pres_mat[missing_genes] <- 0L

pres_mat <- pres_mat[, ordered_genes, drop = FALSE]

# Ensure every isolate in metadata has a row in pres_mat.
# If an isolate is missing here it means its name in the AMR file
# does not match the Isolate column in the metadata (see diagnostics above).
# Those isolates get all-absent rows so the plot never shows NA.
missing_isos <- setdiff(meta$Isolate, rownames(pres_mat))
if (length(missing_isos) > 0) {
  zero_rows <- matrix(0L,
                      nrow     = length(missing_isos),
                      ncol     = ncol(pres_mat),
                      dimnames = list(missing_isos, colnames(pres_mat)))
  pres_mat <- rbind(pres_mat, zero_rows)
  message("ℹ ", length(missing_isos),
          " isolate(s) had no matching rows in the AMR file and were treated as all-absent.\n",
          "  Check for name mismatches printed above.\n",
          "  Affected isolates: ", paste(missing_isos, collapse = ", "))
}

# Stable class → colour mapping
all_classes <- unique(c(amr_genes$Class, point_genes$Class))
class_color_map <- setNames(
  class_palette[(seq_along(all_classes) - 1) %% length(class_palette) + 1],
  all_classes
)

# Gene block factor (for column_split)
gene_block <- factor(
  c(rep("AMR GENES", n_amr_genes), rep("POINT MUTATIONS", n_point_genes)),
  levels = c("AMR GENES", "POINT MUTATIONS")
)

# Class vector for column annotation
class_vec_all <- c(amr_genes$Class, point_genes$Class)

# ============================================================
# 4.  MAIN PLOT FUNCTION
# ============================================================
plot_farm <- function(farm_id, save_pdf = TRUE) {
  
  # ── Subset & sort isolates ─────────────────────────────────
  fm <- meta %>%
    filter(Farm == farm_id) %>%
    arrange(Diet_Stage, Isolate)
  
  isolates  <- fm$Isolate
  diet_vec  <- as.character(fm$Diet_Stage)
  rp_vec    <- fm$RP_label
  st_vec    <- fm$ST
  
  # Matrix for this farm: isolates × genes
  fp_mat <- as.matrix(pres_mat[isolates, , drop = FALSE])
  
  # ── Drop genes absent in every isolate of this farm ────────
  gene_present <- colSums(fp_mat) > 0
  fp_mat       <- fp_mat[, gene_present, drop = FALSE]
  
  # Subset all gene-level vectors to match
  farm_genes      <- ordered_genes[gene_present]
  farm_gene_block <- gene_block[gene_present]
  farm_class_vec  <- class_vec_all[gene_present]
  
  # Rebuild column annotation for this farm's genes only
  farm_col_anno <- HeatmapAnnotation(
    Class = farm_class_vec,
    col   = list(Class = class_color_map),
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    which                = "column",
    height               = unit(5, "mm"),
    annotation_name_side = "left"
  )
  
  n_farm_genes <- ncol(fp_mat)
  
  # ── Clustering within each Diet Stage slice ────────────────
  isolate_block <- factor(diet_vec, levels = diet_stage_order)
  
  make_block_clust_fn <- function(isolate_to_block) {
    function(mat) {
      n    <- nrow(mat)
      big  <- 1e6
      dmat <- matrix(big, n, n, dimnames = list(rownames(mat), rownames(mat)))
      blocks <- isolate_to_block[rownames(mat)]
      lvls   <- unique(blocks[!is.na(blocks)])
      
      for (lv in lvls) {
        idx <- which(blocks == lv)
        # Always fill within-slice block, even for single isolates
        if (length(idx) == 1) {
          dmat[idx, idx] <- 0
          next
        }
        sub <- mat[idx, , drop = FALSE]
        # Use ALL gene columns (not just variable ones) for the distance —
        # shared absence is meaningful for grouping isolates with no genes
        d_sub <- tryCatch(
          as.matrix(dist(sub, method = "binary")),
          error = function(e) as.matrix(dist(sub, method = "euclidean"))
        )
        # binary dist returns NaN when a row is all-zero (0/0);
        # replace NaN with 0 so those isolates cluster together
        d_sub[is.nan(d_sub) | !is.finite(d_sub)] <- 0
        dmat[idx, idx] <- d_sub
      }
      diag(dmat) <- 0
      hclust(as.dist(dmat), method = "complete")
    }
  }
  
  isolate_to_block <- setNames(diet_vec, isolates)
  clust_fn <- make_block_clust_fn(isolate_to_block)
  
  # ── Row annotation (left): Diet Stage, RP, ST ──────────────
  farm_diet_cols <- diet_colors[intersect(diet_stage_order, unique(diet_vec))]
  farm_rp_cols   <- custom_rp_colors[intersect(names(custom_rp_colors), unique(na.omit(rp_vec)))]
  farm_st_cols   <- custom_st_colors[intersect(names(custom_st_colors), unique(st_vec))]
  
  row_anno <- rowAnnotation(
    `Diet Stage`         = diet_vec,
    `Resistance Profile` = rp_vec,
    `ST`                 = st_vec,
    col = list(
      `Diet Stage`         = farm_diet_cols,
      `Resistance Profile` = farm_rp_cols,
      `ST`                 = farm_st_cols
    ),
    annotation_name_rot  = 90,
    annotation_name_gp   = gpar(fontsize = 8, fontface = "bold"),
    annotation_name_side = "top",
    show_legend          = TRUE,
    width                = unit(rep(5, 3), "mm"),
    gap                  = unit(1, "mm")
  )
  
  # ── Heatmap ────────────────────────────────────────────────
  # - row_names_side = "left"  → isolate labels on left, next to dendrogram
  # - No fixed width/height   → heatmap scales to fill the PDF page
  ht <- Heatmap(
    fp_mat,
    
    col = c("0" = "#f0f0f0", "1" = "#2171b5"),
    
    # ── Rows = isolates ──────────────────────────────────────
    row_split            = isolate_block,
    cluster_row_slices   = FALSE,
    cluster_rows         = clust_fn,
    show_row_dend        = TRUE,
    row_dend_side        = "left",
    row_dend_width       = unit(15, "mm"),
    row_dend_reorder     = FALSE,
    row_title_gp         = gpar(fontsize = 0),  # hide Diet Stage slice labels
    row_title_rot        = 0,
    row_title_side       = "right",
    show_row_names       = TRUE,
    row_names_side       = "left",          # names on left, beside dendrogram
    row_names_gp         = gpar(fontsize = 7),
    row_gap              = unit(3, "mm"),
    
    # ── Columns = genes ───────────────────────────────────────
    column_split          = farm_gene_block,
    cluster_column_slices = FALSE,
    cluster_columns       = FALSE,
    column_title_gp       = gpar(fontsize = 10, fontface = "bold",
                                 col = c("#2171b5", "#d7301f")),
    column_title_side     = "top",
    show_column_names     = TRUE,
    column_names_gp       = gpar(fontsize = 8, fontface = "italic"),
    column_names_rot      = 45,
    column_names_side     = "bottom",
    column_gap            = unit(4, "mm"),
    
    # ── Annotations ───────────────────────────────────────────
    left_annotation       = row_anno,
    bottom_annotation     = farm_col_anno,
    
    # ── Appearance ────────────────────────────────────────────
    border   = TRUE,
    rect_gp  = gpar(col = "white", lwd = 0.6),
    
    # ── Legend ────────────────────────────────────────────────
    heatmap_legend_param = list(
      title     = "Gene",
      labels    = c("Absent", "Present"),
      at        = c(0, 1),
      legend_gp = gpar(fill = c("#f0f0f0", "#2171b5")),
      border    = "grey40"
    )
    # No fixed width/height → ComplexHeatmap fills available space
  )
  
  # ── Class colour legend (included in merge_legend) ─────────
  farm_classes <- unique(farm_class_vec)
  class_legend <- Legend(
    labels    = farm_classes,
    title     = "AMR Class",
    legend_gp = gpar(fill = class_color_map[farm_classes]),
    border    = "grey40",
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    nrow      = 2
  )
  
  plot_title <- paste0("Farm ", farm_id,
                       " \u2014 AMR Gene & Point Mutation Presence/Absence")
  
  draw_heatmap <- function() {
    # merge_legend = TRUE folds all annotation legends + class legend together
    draw(
      ht,
      annotation_legend_list = list(class_legend),
      merge_legend           = TRUE,
      heatmap_legend_side    = "bottom",
      annotation_legend_side = "bottom",
      padding                = unit(c(5, 15, 10, 5), "mm"),  # extra left padding for y title
      auto_adjust            = TRUE
    )
    # Add axis titles manually with grid
    # Y-axis title ("Isolates") — rotated, on the far left
    grid.text(
      "Isolates",
      x      = unit(4, "mm"),
      y      = unit(0.5, "npc"),
      rot    = 90,
      gp     = gpar(fontsize = 11, fontface = "bold")
    )
    # X-axis title ("Genes") — centred above the legend block
    grid.text(
      "Genes",
      x      = unit(0.5, "npc"),
      y      = unit(40, "mm"),
      rot    = 0,
      gp     = gpar(fontsize = 11, fontface = "bold")
    )
  }
  
  # Preview
  if (dev.cur() == 1) dev.new()
  draw_heatmap()
  message("  Previewing Farm ", farm_id, " \u2014 press Enter to continue...")
  readline()
  
  # Save as A4 landscape PDF (297 x 210 mm = 11.69 x 8.27 in)
  # For portrait A4 swap width/height: width=8.27, height=11.69
  if (save_pdf) {
    out_file <- paste0("heatmap_farm_faecium", farm_id, ".pdf")
    pdf(out_file, width = 16.54, height = 11.69, paper = "a4r")   # A4 landscape
    draw_heatmap()
    dev.off()
    message("  \u2713 Saved \u2192 ", out_file)
  }
}

# ============================================================
# 5.  RUN FOR ALL FARMS
#     Set save_pdf = FALSE if you only want to preview
# ============================================================
farms <- sort(unique(meta$Farm))
message("Farms found: ", paste(farms, collapse = ", "))

for (f in farms) {
  plot_farm(f, save_pdf = TRUE)
}

message("All done!")

