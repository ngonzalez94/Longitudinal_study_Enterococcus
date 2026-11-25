# =============================================================================
# Plasmid / Contig Gene Map Plotter using gggenes
# Input: Bakta TSV annotation files
# =============================================================================

library(gggenes)
library(ggplot2)
library(dplyr)
library(purrr)

# --- 1. Folder containing your Bakta*.tsv files ---
input_folder <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/E.coli isolates/Plasmids/"

# --- 2. Custom display names ---
# Key = filename WITHOUT extension, Value = label on y-axis
contig_labels <- c(
  "Bakta_tet_amino_phen"       = "Tetracycline-Aminoglycosides-Phenicol",
  "Bakta_sulfo_amino_trimet"   = "Sulfonamides-Aminoglycosides-Trimethoprim",
  "Bakta_beta_amino_sulfo"     = "Beta-lactams-Aminoglycosides-Sulfonamides",
  "Bakta_sulfo_amino_tet_sul1" = "Sulfonamides-Aminoglycosides-Tetracycline (A)",
  "Bakta_sulfo_amino_tet_sul2" = "Sulfonamides-Aminoglycosides-Tetracycline (B)"
)

# --- 3. Top-to-bottom plotting order (use display names above) ---
contig_order <- c(
  "Sulfonamides-Aminoglycosides-Tetracycline (B)",
  "Sulfonamides-Aminoglycosides-Tetracycline (A)",
  "Beta-lactams-Aminoglycosides-Sulfonamides",
  "Sulfonamides-Aminoglycosides-Trimethoprim",
  "Tetracycline-Aminoglycosides-Phenicol"
)

# --- 4. Gene lists per class ---
# Exact gene names as they appear in your Gene column (case-insensitive matching)

resistance_genes <- c(
  "sul1", "sul2", "sul3",
  "aadA1", "aadA2", "aadB",
  "tetR(A)", "tet(A)", "tet(B)", "tetC", "tetD",
  "blaTEM-1", "catA1",
  "dfrA1", "dfrA5", "dfrA57",
  "armA", "aac(6')", "aph(3'')-Ib", "aph(6)-Id"
)

transposon_genes <- c(
  "Transposase", "tnp", "istB",
  "Tn3-like element",
  "tnsA", "tnsB", "tnsC", "tnsD", "tnsE",
  "Transposon resolvase"
)

integron_genes <- c(
  "intI1", "intI2", "intI3",
  "Integrase catalytic region"
)

transcription_genes <- c(
  "lysR", "arsR"
)

# --- 5. Colours ---
# Keys MUST match the internal class names exactly:
# "resistance", "transposon", "integron", "transcription", "other"
class_colours <- c(
  "resistance"    = "darkred",
  "transposon"    = "darkorange",
  "integron"      = "darkblue",
  "transcription" = "darkgreen",
  "other"         = "#BBBBBB"
)

# --- 6. What text to show on gene arrows ---
# "gene_name" | "product" | "locus_tag"
label_column <- "gene_name"

# --- 7. Plot aesthetics ---
plot_title               <- "Bakta Plasmid / Contig Comparison"
arrow_height             <- 4     # arrow height in mm (try 3–6)
font_size_genes          <- 2.5   # size of gene name labels
font_size_axis           <- 10
output_width             <- 16
output_height_per_contig <- 3.0   # increased slightly to give labels more room

# How far labels sit from the arrow (in vjust units):
#   label_above — positive number, labels above the arrow (e.g. 1.2, 2.0)
#   label_below — negative number, labels below the arrow (e.g. -1.2, -2.0)
# Odd-numbered genes go above, even-numbered go below, so long names don't overlap.
label_above <-  1.5
label_below <- -1.5

# =============================================================================
# MAIN SCRIPT — do not edit below unless you know what you're changing
# =============================================================================

# Pre-lowercase gene lists once (used inside classify_genes)
res_l  <- tolower(trimws(resistance_genes))
tnp_l  <- tolower(trimws(transposon_genes))
intg_l <- tolower(trimws(integron_genes))
trn_l  <- tolower(trimws(transcription_genes))

# Internal class names — these MUST match the keys in class_colours above
class_levels <- c("resistance", "transposon", "integron", "transcription", "other")

classify_genes <- function(gene_vec) {
  g   <- tolower(trimws(as.character(gene_vec)))
  out <- rep("other", length(g))
  out[g %in% res_l]                   <- "resistance"
  out[g %in% tnp_l  & out == "other"] <- "transposon"
  out[g %in% intg_l & out == "other"] <- "integron"
  out[g %in% trn_l  & out == "other"] <- "transcription"
  out[is.na(g) | g == ""]             <- "other"
  factor(out, levels = class_levels)
}

# Legend display names — keys match class_levels
legend_labels <- c(
  resistance    = "Resistance genes",
  transposon    = "Transposon / IS element",
  integron      = "Integron / integrase region",
  transcription = "Transcriptional regulator",
  other         = "Other / hypothetical"
)

# Smart reader: skips Bakta comment lines before the header
read_bakta_tsv <- function(filepath) {
  raw_lines  <- readLines(filepath, warn = FALSE)
  header_idx <- which(!grepl("^#", raw_lines))[1]
  if (is.na(header_idx))
    stop("No non-comment lines found in: ", basename(filepath))
  message("  Skipped ", header_idx - 1, " comment line(s)")
  data_lines <- raw_lines[header_idx:length(raw_lines)]
  read.table(text = paste(data_lines, collapse = "\n"),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", fill = TRUE, check.names = FALSE)
}

# Discover files
txt_files <- list.files(path = input_folder,
                        pattern = "^Bakta.*\\.(txt|tsv)$",
                        full.names = TRUE)
if (length(txt_files) == 0)
  stop("No Bakta*.txt or Bakta*.tsv files found in: ", input_folder)

message("Found ", length(txt_files), " file(s):")
message(paste(" -", basename(txt_files), collapse = "\n"))

# Process one file
read_contig_file <- function(filepath) {
  df <- read_bakta_tsv(filepath)
  message("File: ", basename(filepath))
  message("  Raw columns: ", paste(names(df), collapse = " | "))
  
  # Rename columns robustly
  cn <- tolower(trimws(gsub("#", "", names(df))))
  rn <- function(patterns, new_name) {
    i <- which(cn %in% patterns)
    if (length(i)) names(df)[i[1]] <<- new_name
  }
  rn(c("sequence id","sequence_id","seqid","contig","molecule","chr","replicon"), "molecule")
  rn(c("start"),                                                                   "start")
  rn(c("stop","end"),                                                              "end")
  rn(c("strand"),                                                                  "strand")
  rn(c("gene"),                                                                    "gene_name")
  rn(c("product"),                                                                 "product")
  rn(c("locus tag","locus_tag","locustag"),                                        "locus_tag")
  rn(c("type","feat_type","feature"),                                              "feat_type")
  
  message("  Renamed: ", paste(names(df), collapse = " | "))
  
  if (!"molecule" %in% names(df))
    stop("Could not find Sequence Id column in: ", basename(filepath),
         "\nColumns: ", paste(names(df), collapse = ", "))
  
  # CDS only
  if ("feat_type" %in% names(df))
    df <- df[tolower(trimws(df$feat_type)) == "cds", ]
  message("  CDS rows: ", nrow(df))
  if (nrow(df) == 0) { warning("No CDS rows in: ", basename(filepath)); return(NULL) }
  
  # Display label
  base_name <- tools::file_path_sans_ext(basename(filepath))
  df$molecule <- if (base_name %in% names(contig_labels)) {
    contig_labels[[base_name]]
  } else {
    message("  NOTE: '", base_name, "' not in contig_labels — using filename")
    base_name
  }
  
  # Gene names
  if (!"gene_name" %in% names(df)) df$gene_name <- ""
  df$gene_name[is.na(df$gene_name)] <- ""
  
  # Classify
  df$class <- classify_genes(df$gene_name)
  
  # Labels — only non-"other" genes get a label
  raw_label <- switch(label_column,
                      gene_name = df$gene_name,
                      product   = if ("product"   %in% names(df)) df$product   else df$gene_name,
                      locus_tag = if ("locus_tag" %in% names(df)) df$locus_tag else df$gene_name,
                      df$gene_name
  )
  raw_label[is.na(raw_label)] <- ""
  df$label_text <- ifelse(as.character(df$class) == "other", "", raw_label)
  
  # Strand
  df$strand_fwd <- trimws(as.character(df$strand)) %in% c("+", "1", "forward", "Forward")
  
  # Positions
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  
  df[, c("molecule","start","end","strand_fwd","gene_name","class","label_text")]
}

# Load all
all_list <- lapply(txt_files, read_contig_file)
all_list <- Filter(Negate(is.null), all_list)
all_data <- do.call(rbind, all_list)

# Re-enforce factor levels after rbind
all_data$class <- factor(as.character(all_data$class), levels = class_levels)

message("\nTotal rows: ", nrow(all_data))
message("Unique contigs: ", paste(unique(all_data$molecule), collapse = ", "))
message("Class breakdown:"); print(table(all_data$class, useNA = "always"))
if (nrow(all_data) == 0) stop("No data loaded.")

# Plotting order
present       <- unique(as.character(all_data$molecule))
extras        <- setdiff(present, contig_order)
if (length(extras)) {
  warning("Contigs not in contig_order, appended at bottom: ", paste(extras, collapse = ", "))
  contig_order <- c(contig_order, extras)
}
ordered       <- contig_order[contig_order %in% present]
all_data$molecule <- factor(all_data$molecule, levels = rev(ordered))

# Plot
p <- ggplot(all_data,
            aes(xmin    = start,
                xmax    = end,
                y       = molecule,
                fill    = class,
                label   = label_text,
                forward = strand_fwd)) +
  geom_gene_arrow(
    arrowhead_height  = unit(arrow_height,       "mm"),
    arrowhead_width   = unit(arrow_height * 0.6, "mm"),
    arrow_body_height = unit(arrow_height * 0.7, "mm")
  ) +
  # Labels alternate above/below to avoid overlapping long gene names.
  # Odd-indexed genes (within each contig) go above; even go below.
  geom_text(
    data = {
      lab <- subset(all_data, label_text != "")
      # Number genes within each contig in order of position
      lab$gene_idx <- ave(lab$start, lab$molecule,
                          FUN = function(x) rank(x, ties.method = "first"))
      lab$vjust_val <- ifelse(lab$gene_idx %% 2 == 1, label_above, label_below)
      lab
    },
    aes(x     = (start + end) / 2,
        y     = molecule,
        label = label_text,
        vjust = vjust_val),
    size         = font_size_genes,
    fontface     = "italic",
    inherit.aes  = FALSE
  ) +
  facet_wrap(~ molecule, ncol = 1, scales = "free_x") +
  scale_fill_manual(
    values   = class_colours,       # keys: resistance/transposon/integron/transcription/other
    labels   = legend_labels,       # human-readable names for the legend
    name     = "Gene class",
    na.value = "#BBBBBB",           # stray NAs → grey
    drop     = TRUE
  ) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = 0.02)) +
  labs(title = plot_title, x = "Position (bp)", y = NULL) +
  theme_genes() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0),
    axis.text.y     = element_text(size = font_size_axis, face = "italic"),
    legend.position = "bottom",
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 9),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    panel.spacing.y  = unit(8, "mm")   # vertical gap between contig panels
  )

print(p)

# Save
n_contigs   <- length(levels(all_data$molecule))
output_file <- file.path(input_folder, "Bakta_contig_map.pdf")
ggsave(filename = output_file, plot = p,
       width  = output_width,
       height = output_height_per_contig * n_contigs + 1.5,
       units  = "in", device = "pdf")
message("Plot saved to: ", output_file)