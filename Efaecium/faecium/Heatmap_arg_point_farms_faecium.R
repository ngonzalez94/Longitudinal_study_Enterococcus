# ==============================================================================
# Antimicrobial Resistance Gene Heatmap with Metadata Annotations
# Updated: Cell colors reflect detection source (Abricate/AMRFinder/Both)
#          Point mutations shown as a separate column block (AMRFinder, blue)
#          Resistance Profile codes A-G shown as descriptive labels in legend
# ==============================================================================

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(dplyr)
library(tidyr)
library(grid)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

main_data <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/ARG_faecium.xlsx", sheet = 1)

gene_classes <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/ARG_faecium.xlsx", sheet = 2)

amrfinder_genes_raw     <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/amrfinder_faecium.xlsx", sheet = 2)
amrfinder_mutations_raw <- read_excel("C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/amrfinder_faecium.xlsx", sheet = 1)

# =============================================================================
# CUSTOMIZE COLORS
# =============================================================================

# Resistance Profile: colors keyed on DESCRIPTIVE labels (used in legend)
custom_rp_colors <- c(
  "CIP-ERY-TET-TYL"     = "#e41a1c",
  "CIP-TET-TYL -AMP" = "#377eb8"
)

# Sequence Type colors
custom_st_colors <- c(
  "771"   = "#a50026", "870"  = "#fee08b", "22"  = "#66bd63", "-"   = "grey")

# Cell color scheme: 0=absent, 1=Abricate only, 2=AMRFinder only, 3=Both
source_colors <- colorRamp2(c(0, 1, 2, 3), c("white", "#FFD700", "#2166AC", "#33A02C"))

# SAVE FILES OR JUST PREVIEW?
save_files <- TRUE

# ==============================================================================
# 2. PREPARE ABRICATE (main_data) PRESENCE MATRIX
# ==============================================================================

metadata_cols <- c("Resistance_Profile", "Farm", "Diet_Stage", "ST")

# Mapping from Excel codes (A-G) to descriptive resistance profile labels
rp_label_map <- c(
  "A" = "CIP-ERY-TET-TYL",
  "B" = "CIP-TET-TYL-AMP"
)

# --- AUTO-DETECT ISOLATE COLUMN -----------------------------------------------
candidate_cols <- c("#FILE", "Isolate", "Sample", "Sample Name", "isolate", "sample",
                    "SampleName", "Strain", "strain", "ID", "id")
isolate_col <- NULL
for (cand in candidate_cols) {
  if (cand %in% colnames(main_data)) { isolate_col <- cand; break }
}
if (is.null(isolate_col)) {
  cat("ERROR: Could not auto-detect the isolate/sample column.\n")
  cat("Columns present in main_data:\n")
  print(colnames(main_data))
  stop("Set isolate_col manually to the correct column name shown above.")
}
cat("Auto-detected isolate column: '", isolate_col, "'\n\n", sep = "")
# ------------------------------------------------------------------------------

# Clean isolate names (remove _contig.fa suffix if present)
isolates_raw <- as.character(main_data[[isolate_col]])
isolates     <- gsub("_contig\\.fa$", "", isolates_raw)

# --- DIAGNOSTICS ---------------------------------------------------------------
cat("=== Data Diagnostics ===\n")
cat("Rows in main_data       :", nrow(main_data), "\n")
cat("Isolate names extracted :", length(isolates), "\n")
cat("NA isolate names        :", sum(is.na(isolates)), "\n")
cat("Empty isolate names     :", sum(nchar(trimws(isolates)) == 0, na.rm = TRUE), "\n")
cat("Duplicate isolate names :", sum(duplicated(isolates)), "\n")
if (any(duplicated(isolates))) {
  cat("  Duplicates found:\n")
  print(isolates[duplicated(isolates)])
}
cat("First 5 isolate names   :", paste(head(isolates, 5), collapse = ", "), "\n\n")
# ------------------------------------------------------------------------------

# Convert to plain data.frame (avoids tibble rowname issues)
metadata <- as.data.frame(main_data[, metadata_cols])
rownames(metadata) <- isolates

# Translate Resistance_Profile codes to descriptive labels
metadata$Resistance_Profile <- rp_label_map[as.character(metadata$Resistance_Profile)]
# Preserve any codes not in the map as-is
metadata$Resistance_Profile[is.na(metadata$Resistance_Profile)] <-
  as.character(main_data$Resistance_Profile[is.na(metadata$Resistance_Profile)])

exclude_cols <- c(isolate_col, "NUM_FOUND", metadata_cols)
gene_cols    <- setdiff(colnames(main_data), exclude_cols)
gene_data    <- as.data.frame(main_data[, gene_cols])

# Abricate presence/absence (1 = present, 0 = absent)
abricate_matrix <- apply(gene_data, 2, function(x) ifelse(x == "." | is.na(x), 0L, 1L))
rownames(abricate_matrix) <- isolates

# ==============================================================================
# 3. PREPARE AMRFINDER GENE PRESENCE MATRIX
# ==============================================================================

amrfinder_genes_raw$Sample <- gsub("_contig\\.fa$", "", as.character(amrfinder_genes_raw$Sample))

cat("=== AMRFinder Gene Data Diagnostics ===\n")
cat("Unique samples in AMRFinder genes sheet:", length(unique(amrfinder_genes_raw$Sample)), "\n")
cat("Samples in main_data but NOT in AMRFinder:",
    sum(!isolates %in% unique(amrfinder_genes_raw$Sample)), "\n")
cat("Samples in AMRFinder but NOT in main_data:",
    sum(!unique(amrfinder_genes_raw$Sample) %in% isolates), "\n\n")

amr_gene_long <- amrfinder_genes_raw %>%
  select(Sample, `Gene symbol`) %>%
  distinct() %>%
  mutate(present = 1L)

amr_gene_wide <- amr_gene_long %>%
  pivot_wider(names_from = `Gene symbol`, values_from = present, values_fill = 0L)

# Ensure every isolate row exists (left join from canonical list)
all_isolates_df <- data.frame(Sample = isolates, stringsAsFactors = FALSE)
amr_gene_wide   <- all_isolates_df %>% left_join(amr_gene_wide, by = "Sample")

amr_gene_matrix <- as.matrix(amr_gene_wide[, -1, drop = FALSE])
amr_gene_matrix[is.na(amr_gene_matrix)] <- 0L
rownames(amr_gene_matrix) <- isolates

# ==============================================================================
# 4. BUILD COMBINED SOURCE MATRIX FOR ARG GENES
#    0 = absent, 1 = Abricate only, 2 = AMRFinder only, 3 = Both
# ==============================================================================

all_arg_genes <- union(colnames(abricate_matrix), colnames(amr_gene_matrix))

expand_matrix <- function(mat, all_genes, row_names) {
  out    <- matrix(0L, nrow = length(row_names), ncol = length(all_genes),
                   dimnames = list(row_names, all_genes))
  shared <- intersect(colnames(mat), all_genes)
  out[, shared] <- mat[, shared]
  out
}

ab_full  <- expand_matrix(abricate_matrix, all_arg_genes, isolates)
amr_full <- expand_matrix(amr_gene_matrix, all_arg_genes, isolates)

# 0=absent | 1=ab only | 2=amr only | 3=both
source_matrix <- ab_full + 2L * amr_full

# ==============================================================================
# 5. PREPARE POINT MUTATIONS MATRIX (sheet 1, always AMRFinder = 2)
# ==============================================================================

amrfinder_mutations_raw$Sample <-
  gsub("_contig\\.fa$", "", as.character(amrfinder_mutations_raw$Sample))

mut_long <- amrfinder_mutations_raw %>%
  select(Sample, `Gene symbol`) %>%
  distinct() %>%
  mutate(present = 1L)

mut_wide <- mut_long %>%
  pivot_wider(names_from = `Gene symbol`, values_from = present, values_fill = 0L)

mut_wide        <- all_isolates_df %>% left_join(mut_wide, by = "Sample")
mut_matrix_raw  <- as.matrix(mut_wide[, -1, drop = FALSE])
mut_matrix_raw[is.na(mut_matrix_raw)] <- 0L
rownames(mut_matrix_raw) <- isolates

# Point mutations are AMRFinder-only: encode as 2 (blue)
mut_matrix <- mut_matrix_raw * 2L

# ==============================================================================
# 6. ORGANIZE GENES BY ANTIMICROBIAL CLASS
# ==============================================================================

# ==============================================================================
# MANUAL CLASS OVERRIDES
# If a gene is not found in Ecoli_mechanisms.xlsx it will be labelled "Unknown".
# Add entries below to manually assign classes to those genes, e.g.:
#   "gyrA"  = "FLUOROQUINOLONE",
#   "parC"  = "FLUOROQUINOLONE"
# ==============================================================================
manual_class_overrides <- c(
  "eat(A)_T450I" = "Pleuromutilin",
  "liaR_E75K" = "Daptomycin"
)

organize_by_class <- function(gene_names, gene_classes_df, overrides = c()) {
  info <- data.frame(Genes = gene_names, stringsAsFactors = FALSE)
  info <- left_join(info, gene_classes_df, by = "Genes")
  # Apply manual overrides first
  if (length(overrides) > 0) {
    override_idx <- info$Genes %in% names(overrides)
    info$Classification[override_idx] <- overrides[info$Genes[override_idx]]
  }
  info$Classification[is.na(info$Classification)] <- "Unknown"
  info[order(info$Classification), ]
}

arg_gene_info <- organize_by_class(all_arg_genes,          gene_classes, manual_class_overrides)
mut_gene_info <- organize_by_class(colnames(mut_matrix),   gene_classes, manual_class_overrides)

# Diagnostic: show which genes ended up as "Unknown"
arg_unknown <- arg_gene_info$Genes[arg_gene_info$Classification == "Unknown"]
mut_unknown <- mut_gene_info$Genes[mut_gene_info$Classification == "Unknown"]
if (length(arg_unknown) > 0) {
  cat("ARG genes with no class match (showing as Unknown):
")
  print(arg_unknown)
}
if (length(mut_unknown) > 0) {
  cat("Point mutation genes with no class match (showing as Unknown):
  Add these to manual_class_overrides or Ecoli_mechanisms.xlsx
")
  print(mut_unknown)
}

source_matrix <- source_matrix[, arg_gene_info$Genes, drop = FALSE]
mut_matrix    <- mut_matrix[,    mut_gene_info$Genes,  drop = FALSE]

# ==============================================================================
# 7. METADATA COLOR SCHEMES
# ==============================================================================

unique_profiles <- sort(unique(metadata$Resistance_Profile))
profile_colors  <- custom_rp_colors
missing_rp <- setdiff(unique_profiles, names(custom_rp_colors))
if (length(missing_rp) > 0) {
  cat("WARNING: Adding auto colors for Resistance Profiles:", paste(missing_rp, collapse = ", "), "\n")
  profile_colors <- c(profile_colors,
                      setNames(hcl.colors(length(missing_rp), "Set 2"), missing_rp))
}

unique_farms <- unique(metadata$Farm)
farm_colors  <- setNames(rainbow(length(unique_farms)), unique_farms)

# Diet Stage colors: Pre-WM = light blue, WM = blue, Post-WM = light green
diet_stage_levels <- c("Pre-WM", "WM", "Post-WM")
unique_diets <- diet_stage_levels
diet_colors  <- c(
  "Pre-WM"  = "#ADD8E6",   # light blue
  "WM"      = "#1F78B4",   # blue
  "Post-WM" = "#B2E0A0"    # light green
)

unique_st <- sort(unique(as.character(metadata$ST)))
st_colors  <- custom_st_colors
missing_st <- setdiff(unique_st, names(custom_st_colors))
if (length(missing_st) > 0) {
  cat("WARNING: Adding auto colors for STs:", paste(missing_st, collapse = ", "), "\n")
  st_colors <- c(st_colors,
                 setNames(hcl.colors(length(missing_st), "Pastel 1"), missing_st))
}

# ==============================================================================
# 8. CUSTOM LEGEND FOR DETECTION SOURCE
# ==============================================================================

source_legend <- Legend(
  title     = "Antimicrobial Detection Source",
  labels    = c("Absent", "Abricate", "AMRFinder", "Both"),
  legend_gp = gpar(fill = c("white", "#FFD700", "#2166AC", "#33A02C"), col = "gray60"),
  border    = "gray60"
)

# ==============================================================================
# 9. HELPER: BUILD ONE COMBINED HEATMAP OBJECT (ARGs + mutations side by side)
# ==============================================================================

make_heatmap <- function(sm, mm, s_gene_info, m_gene_info,
                         meta, row_split_var, row_title_str,
                         show_row_names = TRUE) {
  
  unique_classes <- unique(c(s_gene_info$Classification, m_gene_info$Classification))
  class_colors   <- setNames(rainbow(length(unique_classes)), unique_classes)
  
  row_ha <- rowAnnotation(
    `Resistance Profile` = meta$Resistance_Profile,
    `Diet Stage` = meta$Diet_Stage,
    ST = as.character(meta$ST),
    col = list(`Resistance Profile` = profile_colors, `Diet Stage` = diet_colors, ST = st_colors),
    annotation_name_side = "top",
    annotation_name_gp   = gpar(fontsize = 10),
    simple_anno_size     = unit(0.5, "cm"),
    gap                  = unit(0.5, "mm")
  )
  
  # Colour bar on top of ARGs block — one bar, legend suppressed (built manually below)
  arg_col_ha <- HeatmapAnnotation(
    `Antimicrobial Class` = s_gene_info$Classification,
    col                        = list(`Antimicrobial Class` = class_colors),
    show_annotation_name       = FALSE,   # name hidden; legend shown separately
    simple_anno_size           = unit(0.4, "cm"),
    show_legend                = FALSE
  )
  
  mut_col_ha <- HeatmapAnnotation(
    `Antimicrobial Class` = m_gene_info$Classification,
    col                        = list(`Antimicrobial Class` = class_colors),
    show_annotation_name       = FALSE,
    simple_anno_size           = unit(0.4, "cm"),
    show_legend                = FALSE
  )
  
  ht_arg <- Heatmap(
    sm,
    name                = "Detection Source",
    col                 = source_colors,
    show_heatmap_legend = FALSE,
    row_names_side      = "left",
    row_names_gp        = gpar(fontsize = if (show_row_names) 8 else 1),
    show_row_names      = show_row_names,
    row_title           = row_title_str,
    row_title_gp        = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows        = TRUE,
    show_row_dend       = TRUE,
    row_dend_width      = unit(2, "cm"),
    row_split           = factor(as.character(row_split_var),
                                 levels = c("Pre-WM", "WM", "Post-WM"),
                                 ordered = TRUE),
    row_gap             = unit(0.5, "mm"),
    column_names_side   = "bottom",
    column_names_gp     = gpar(fontsize = 7),
    column_names_rot    = 45,
    column_title        = "ARGs",
    column_title_gp     = gpar(fontsize = 11, fontface = "bold"),
    cluster_columns     = FALSE,
    show_column_dend    = FALSE,
    column_split        = s_gene_info$Classification,
    column_gap          = unit(0.5, "mm"),
    left_annotation     = row_ha,
    top_annotation      = arg_col_ha,
    border              = TRUE,
    rect_gp             = gpar(col = "gray90", lwd = 0.5)
  )
  
  ht_mut <- Heatmap(
    mm,
    name                = "Detection Source (mut)",
    col                 = source_colors,
    show_heatmap_legend = FALSE,
    show_row_names      = FALSE,
    cluster_rows        = FALSE,
    column_names_side   = "bottom",
    column_names_gp     = gpar(fontsize = 7),
    column_names_rot    = 45,
    column_title        = "Point Mutations",
    column_title_gp     = gpar(fontsize = 11, fontface = "bold"),
    cluster_columns     = FALSE,
    show_column_dend    = FALSE,
    column_split        = m_gene_info$Classification,
    column_gap          = unit(0.5, "mm"),
    top_annotation      = mut_col_ha,
    border              = TRUE,
    rect_gp             = gpar(col = "gray90", lwd = 0.5)
  )
  
  # Build the Antimicrobial Class legend manually (covers both ARG + mutation classes)
  # padding pushes it below the other legends on the right
  class_legend <- Legend(
    title     = "Antimicrobial Class",
    labels    = unique_classes,
    legend_gp = gpar(fill = class_colors[unique_classes], col = "gray60"),
    border    = "gray60",
    gap       = unit(8, "mm")   # vertical gap above this legend block
  )
  
  list(heatmap = ht_arg + ht_mut, class_legend = class_legend)
}

# ==============================================================================
# 10. PER-FARM HEATMAPS
# ==============================================================================

metadata$Diet_Stage <- factor(as.character(metadata$Diet_Stage),
                              levels = c("Pre-WM", "WM", "Post-WM"),
                              ordered = TRUE)
farms <- sort(unique(metadata$Farm))

cat("\n=== Creating per-farm heatmaps ===\n")

for (farm in farms) {
  
  cat("\nProcessing Farm", farm, "...\n")
  
  idx          <- which(metadata$Farm == farm)
  farm_meta    <- metadata[idx, , drop = FALSE]
  farm_sm      <- source_matrix[idx, , drop = FALSE]
  farm_mm      <- mut_matrix[idx,    , drop = FALSE]
  
  if (nrow(farm_sm) == 0) { cat("  No isolates — skipping.\n"); next }
  
  sm_keep <- colSums(farm_sm > 0) > 0
  mm_keep <- colSums(farm_mm > 0) > 0
  
  farm_sm       <- farm_sm[, sm_keep, drop = FALSE]
  farm_mm       <- farm_mm[, mm_keep, drop = FALSE]
  farm_arg_info <- arg_gene_info[sm_keep, ]
  farm_mut_info <- mut_gene_info[mm_keep, ]
  
  cat("  ARGs:", sum(sm_keep), " | Mutations:", sum(mm_keep), "\n")
  
  # Re-apply factor with all three levels in correct order.
  # Using factor() with all levels (even absent ones) forces ComplexHeatmap
  # to respect Pre-WM -> WM -> Post-WM regardless of what is present in this farm.
  farm_meta$Diet_Stage <- factor(as.character(farm_meta$Diet_Stage),
                                 levels = c("Pre-WM", "WM", "Post-WM"),
                                 ordered = TRUE)
  
  ht_result <- make_heatmap(
    sm            = farm_sm,
    mm            = farm_mm,
    s_gene_info   = farm_arg_info,
    m_gene_info   = farm_mut_info,
    meta          = farm_meta,
    row_split_var = farm_meta$Diet_Stage,
    row_title_str = paste("Farm", farm, "- Isolates")
  )
  ht           <- ht_result$heatmap
  class_legend <- ht_result$class_legend
  
  draw(ht,
       annotation_legend_side = "right",
       merge_legend           = TRUE,
       annotation_legend_list = list(source_legend, class_legend))
  
  if (save_files) {
    pdf_out <- paste0("ARG_heatmap_Farm_", farm, ".pdf")
    png_out <- paste0("ARG_heatmap_Farm_", farm, ".png")
    pdf(pdf_out, width = 20, height = 12)
    draw(ht, annotation_legend_side = "right", merge_legend = TRUE,
         annotation_legend_list = list(source_legend, class_legend))
    dev.off()
    png(png_out, width = 2000, height = 1200, res = 100)
    draw(ht, annotation_legend_side = "right", merge_legend = TRUE,
         annotation_legend_list = list(source_legend, class_legend))
    dev.off()
    cat("  Saved:", pdf_out, "&", png_out, "\n")
  }
}

# ==============================================================================
# 11. COMBINED ALL-FARMS HEATMAP
# ==============================================================================

cat("\n=== Creating combined all-farms heatmap ===\n")

sm_keep_all <- colSums(source_matrix > 0) > 0
mm_keep_all <- colSums(mut_matrix    > 0) > 0

comb_sm       <- source_matrix[, sm_keep_all, drop = FALSE]
comb_mm       <- mut_matrix[,    mm_keep_all, drop = FALSE]
comb_arg_info <- arg_gene_info[sm_keep_all, ]
comb_mut_info <- mut_gene_info[mm_keep_all, ]

metadata$Farm <- factor(metadata$Farm, levels = sort(unique(metadata$Farm)))

unique_classes_all <- unique(c(comb_arg_info$Classification, comb_mut_info$Classification))
class_colors_all   <- setNames(rainbow(length(unique_classes_all)), unique_classes_all)

row_ha_all <- rowAnnotation(
  Farm = metadata$Farm,
  `Resistance Profile` = metadata$Resistance_Profile,
  `Diet Stage` = metadata$Diet_Stage,
  ST   = as.character(metadata$ST),
  col  = list(Farm = farm_colors, `Resistance Profile` = profile_colors, `Diet Stage` = diet_colors, ST = st_colors),
  annotation_name_side = "top",
  annotation_name_gp   = gpar(fontsize = 10),
  simple_anno_size     = unit(0.5, "cm"),
  gap                  = unit(2, "mm")
)

# Build combined class legend for all farms plot
all_classes_combined <- unique(c(comb_arg_info$Classification, comb_mut_info$Classification))
class_legend_all <- Legend(
  title     = "Antimicrobial Class",
  labels    = all_classes_combined,
  legend_gp = gpar(fill = class_colors_all[all_classes_combined], col = "gray60"),
  border    = "gray60",
  gap       = unit(8, "mm")
)

# Colour bars for combined heatmap (legend suppressed — built manually)
arg_col_ha_all <- HeatmapAnnotation(
  `Antimicrobial Class` = comb_arg_info$Classification,
  col                        = list(`Antimicrobial Class` = class_colors_all),
  show_annotation_name       = FALSE,
  simple_anno_size           = unit(0.4, "cm"),
  show_legend                = FALSE
)

mut_col_ha_all <- HeatmapAnnotation(
  `Antimicrobial Class` = comb_mut_info$Classification,
  col                        = list(`Antimicrobial Class` = class_colors_all),
  show_annotation_name       = FALSE,
  simple_anno_size           = unit(0.4, "cm"),
  show_legend                = FALSE
)

ht_arg_all <- Heatmap(
  comb_sm,
  name                   = "Detection Source",
  col                    = source_colors,
  show_heatmap_legend    = FALSE,
  row_names_side         = "left",
  row_names_gp           = gpar(fontsize = 7),
  row_title              = "Isolates (All Farms)",
  row_title_gp           = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows           = TRUE,
  cluster_row_slices     = FALSE,
  show_row_dend          = TRUE,
  row_dend_width         = unit(2, "cm"),
  row_split              = metadata$Farm,
  row_gap                = unit(3, "mm"),
  column_names_side      = "bottom",
  column_names_gp        = gpar(fontsize = 7),
  column_names_rot       = 45,
  column_title           = "ARGs",
  column_title_gp        = gpar(fontsize = 11, fontface = "bold"),
  cluster_columns        = FALSE,
  show_column_dend       = FALSE,
  column_split           = comb_arg_info$Classification,
  column_gap             = unit(2, "mm"),
  left_annotation        = row_ha_all,
  top_annotation         = arg_col_ha_all,
  border                 = TRUE,
  rect_gp                = gpar(col = "gray90", lwd = 0.5)
)

ht_mut_all <- Heatmap(
  comb_mm,
  name                = "Detection Source (mut)",
  col                 = source_colors,
  show_heatmap_legend = FALSE,
  show_row_names      = FALSE,
  cluster_rows        = FALSE,
  column_names_side   = "bottom",
  column_names_gp     = gpar(fontsize = 7),
  column_names_rot    = 45,
  column_title        = "Point Mutations",
  column_title_gp     = gpar(fontsize = 11, fontface = "bold"),
  cluster_columns     = FALSE,
  show_column_dend    = FALSE,
  column_split        = comb_mut_info$Classification,
  column_gap          = unit(2, "mm"),
  top_annotation      = mut_col_ha_all,
  border              = TRUE,
  rect_gp             = gpar(col = "gray90", lwd = 0.5)
)

ht_combined_all <- ht_arg_all + ht_mut_all

draw(ht_combined_all,
     annotation_legend_side = "right",
     merge_legend           = TRUE,
     annotation_legend_list = list(source_legend, class_legend_all))

if (save_files) {
  pdf("ARG_heatmap_ALL_FARMS.pdf", width = 24, height = 16)
  draw(ht_combined_all, annotation_legend_side = "right", merge_legend = TRUE,
       annotation_legend_list = list(source_legend, class_legend_all))
  dev.off()
  png("ARG_heatmap_ALL_FARMS.png", width = 2400, height = 1600, res = 100)
  draw(ht_combined_all, annotation_legend_side = "right", merge_legend = TRUE,
       annotation_legend_list = list(source_legend, class_legend_all))
  dev.off()
  cat("Saved: ARG_heatmap_ALL_FARMS.pdf & ARG_heatmap_ALL_FARMS.png\n")
}

cat("\n=== ALL DONE! ===\n")