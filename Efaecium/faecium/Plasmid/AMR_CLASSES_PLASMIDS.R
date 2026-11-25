# =============================================================================
# Plasmid / Contig Gene Map Plotter using gggenes
# Reads all txt files starting with "Bkata" from a specified folder
# Features: custom contig names, custom class colours, custom order
# =============================================================================

library(gggenes)
library(ggplot2)
library(dplyr)
library(purrr)

# =============================================================================
# !! USER SETTINGS — Edit everything in this section !!
# =============================================================================

# --- 1. Folder containing your Bkata*.txt files ---
input_folder <- "C:/Users/ngonzalez/OneDrive - Massey University/Thesis/Chapter 2 - Longitudinal Study/Enterococcus/faecium/Plasmid/"   # <-- change this

# --- 2. Custom display names for each contig ---
# Key   = filename WITHOUT extension (e.g. "Bkata_contig1")
# Value = label you want shown on the plot
contig_labels <- c(
  "Bakta_amino_linco_streptogramin_1" = "Aminoglycosides- Lincosamides- Streptogramin (A)",
  "Bakta_amino_linco_streptogramin_2" = "Aminoglycosides- Lincosamides- Streptogramin (B)",
  "Bakta_linco_streptogramin_aminoglyco" = "Lincosamides- Streptogramin- Aminoglyscosides",
  "Bakta_tetra_amino_mls" = "Tetracycline- Aminoglycosides- MLS" 
  
  
  # Add / remove lines to match your files
)

# --- 3. Desired plotting order (use the display names from contig_labels) ---
# The first name in the vector will appear at the TOP of the plot.

contig_order <- c(
  "Aminoglycosides- Lincosamides- Streptogramin (A)",
  "Aminoglycosides- Lincosamides- Streptogramin (B)",
  "Lincosamides- Streptogramin- Aminoglyscosides",
  "Tetracycline- Aminoglycosides- MLS" 
)
#--- 4. Gene class definitions ---
# List every Gene name (from the "Gene" column) that belongs to each class.
# Matching is case-insensitive. Must be the exact gene name as in your file.
# Genes not listed in any class → assigned to "other" (grey).
resistance_genes <- c(
  "sul1", "sul2", "sul3",
  "aadA1", "aadA2", "aadB",
  "tetR(A)", "tet(A)", "tet(B)", "tetC", "tetD","aph(3')-IIIa",
  "blaTEM-1", "catA1","erm(B)","tet(S)","tet(M)","tet(L)","ant(6)-Ia", "lsa(E)","lnu(B)",
  "dfrA1", "dfrA5", "dfrA57","sat4","spw",
  "armA", "aac(6')", "aph(3'')-Ib", "aph(6)-Id"
  # Add more resistance gene names as needed
)

transposon_genes <- c(
  "Transposase", "tnp", "istB",
  "Tn3-like element", "Putative transposase",
  "tnsA", "tnsB", "tnsC", "tnsD", "tnsE", "Transposon resolvase", "Tn3 resolvase"
  # Add more transposon-related gene names as needed
)

integron_genes <- c(
  "intI1", "intI2", "intI3", "Integrase catalytic region"
  # Add more integron gene names as needed
)

transcriptional_regulator <- c(
  "lysR", "arsR"
  # Add other transcriptional element / plasmid backbone gene names as needed
)

origin_replication <- c(
  "Origin of replication"
  # Add other transcriptional element / plasmid backbone gene names as needed
)


# --- 5. Colours for each class ---
class_colours <- c(
  "Resistance genes"  = "darkred",   # red
  "Transposon genes"  = "darkorange",   # orange
  "Integrons/integrase region"    = "darkblue",   # teal
  "Transcriptional regulators"      = "darkgreen",   # steel blue
  "Origin of replication" = "purple",
  "Other"       = "#BBBBBB"    # grey
)

# --- 6. What to label genes with on the arrows ---
# Options: "gene_name"  — uses the Gene column (blank for hypotheticals)
#          "product"    — uses the Product column (always has text)
#          "locus_tag"  — uses the Locus Tag column
label_column <- "gene_name"   # <-- change to "product" or "locus_tag" if preferred

# --- 7. Plot aesthetics ---
plot_title               <- "Bakta Plasmid / Contig Comparison"
arrow_height             <- 8   # height of gene arrows (0–1)
font_size_genes          <- 4   # size of gene labels inside arrows
font_size_axis           <- 13    # size of y-axis contig labels
output_width             <- 16    # PDF width in inches
output_height_per_contig <- 2.5   # PDF height per contig (inches)


label_above <-  3
label_below <- -2.5
# =============================================================================
# MAIN SCRIPT
# =============================================================================

# Pre-compute lowercased gene lists once so the lookup is always fast & safe
res_l  <- tolower(trimws(resistance_genes))
tnp_l  <- tolower(trimws(transposon_genes))
intg_l <- tolower(trimws(integron_genes))
trn_l  <- tolower(trimws(transcriptional_regulator))
ori_l  <- tolower(trimws(origin_replication))

classify_genes <- function(gene_vec) {
  # Vectorised classification — no sapply, no scoping issues
  g <- tolower(trimws(as.character(gene_vec)))
  out <- rep("Other", length(g))
  out[g %in% res_l]                    <- "Resistance genes"
  out[g %in% tnp_l & out == "Other"]   <- "Transposon genes"
  out[g %in% intg_l & out == "Other"]  <- "Integrons/integrase region"
  out[g %in% trn_l & out == "Other"]   <- "Transcriptional regulators"
  out[g %in% ori_l & out == "Other"]   <- "Origin of replication"
  out[is.na(g) | g == ""]              <- "Other"
  # Return as factor with fixed levels — survives bind_rows/map_dfr
  factor(out, levels = c("Resistance genes", "Transposon genes", "Integrons/integrase region", "Transcriptional regulators", "Origin of replication", "Other"))
}


# -- Smart file reader: skips leading Bakta comment lines --
read_bakta_tsv <- function(filepath) {
  raw_lines  <- readLines(filepath, warn = FALSE)
  header_idx <- which(!grepl("^#", raw_lines))[1]
  if (is.na(header_idx))
    stop("No non-comment lines found in: ", basename(filepath))
  message("  Skipped ", header_idx - 1, " comment line(s) before header")
  data_lines <- raw_lines[header_idx:length(raw_lines)]
  read.table(text = paste(data_lines, collapse = "\n"),
             header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             quote = "", fill = TRUE, check.names = FALSE)
}

# -- Discover files --
txt_files <- list.files(path = input_folder,
                        pattern = "^Bakta.*\\.(txt|tsv)$",
                        full.names = TRUE)

if (length(txt_files) == 0)
  stop("No files matching 'Bkata*.txt' or 'Bkata*.tsv' found in: ", input_folder)

message("Found ", length(txt_files), " file(s):")
message(paste(" -", basename(txt_files), collapse = "\n"))

# -- Read and process each file --
read_contig_file <- function(filepath) {
  
  df <- read_bakta_tsv(filepath)
  
  message("File: ", basename(filepath))
  message("  Raw columns: ", paste(names(df), collapse = " | "))
  
  # Robust column renaming (strip #, lowercase, trim)
  clean_name <- function(x) tolower(trimws(gsub("#", "", x)))
  cn <- clean_name(names(df))
  
  rename_col <- function(patterns, new_name) {
    idx <- which(cn %in% patterns)
    if (length(idx) > 0) names(df)[idx[1]] <<- new_name
  }
  rename_col(c("sequence id","sequence_id","seqid","contig","molecule","chr","replicon"), "molecule")
  rename_col(c("start"),                                                                   "start")
  rename_col(c("stop","end"),                                                              "end")
  rename_col(c("strand"),                                                                  "strand")
  rename_col(c("gene"),                                                                    "gene_name")
  rename_col(c("product"),                                                                 "product")
  rename_col(c("locus tag","locus_tag","locustag"),                                        "locus_tag")
  rename_col(c("type","feat_type","feature"),                                              "feat_type")
  
  message("  Renamed columns: ", paste(names(df), collapse = " | "))
  
  if (!"molecule" %in% names(df))
    stop("Could not find a molecule/sequence-ID column in: ", basename(filepath),
         "\nColumns found: ", paste(names(df), collapse = ", "))
  
  # CDS only
  if ("feat_type" %in% names(df)) {
    df <- df[tolower(trimws(df$feat_type)) == "cds", ]
    message("  CDS rows: ", nrow(df))
  }
  if (nrow(df) == 0) { warning("No CDS rows in: ", basename(filepath)); return(NULL) }
  
  # Display label
  base_name <- tools::file_path_sans_ext(basename(filepath))
  display_label <- if (base_name %in% names(contig_labels)) {
    contig_labels[[base_name]]
  } else {
    message("  NOTE: '", base_name, "' not in contig_labels — using filename")
    base_name
  }
  df$molecule <- display_label
  
  # Gene name column (fill blanks with "")
  if (!"gene_name" %in% names(df)) df$gene_name <- ""
  df$gene_name[is.na(df$gene_name)] <- ""
  
  # Class — vectorised, returns a proper factor
  df$class <- classify_genes(df$gene_name)
  
  # Label: show gene name only for classified (non-other) genes
  raw_label <- switch(label_column,
                      gene_name  = df$gene_name,
                      product    = if ("product"   %in% names(df)) df$product   else df$gene_name,
                      locus_tag  = if ("locus_tag" %in% names(df)) df$locus_tag else df$gene_name,
                      df$gene_name
  )
  raw_label[is.na(raw_label)] <- ""
  df$label_text <- ifelse(as.character(df$class) == "Other", "", raw_label)
  
  # Strand → logical (TRUE = forward)
  df$strand_fwd <- trimws(as.character(df$strand)) %in% c("+", "1", "forward", "Forward")
  
  # Numeric positions
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  
  # Return only the columns we need (keeps bind_rows clean)
  df[, c("molecule","start","end","strand_fwd","gene_name","class","label_text")]
}

# -- Load all files --
all_list <- lapply(txt_files, read_contig_file)
all_list <- Filter(Negate(is.null), all_list)  # drop any empty files
all_data <- do.call(rbind, all_list)            # base-R bind — preserves factors

# Re-apply factor (do.call/rbind can sometimes drop levels too)
all_data$class <- factor(as.character(all_data$class),
                         levels = c("Resistance genes","Transposon genes","Integrons/integrase region","Transcriptional regulators","Origin of replication", "Other"))

message("\nTotal rows: ", nrow(all_data))
message("Unique contigs: ", paste(unique(all_data$molecule), collapse = ", "))
message("Class breakdown:"); print(table(all_data$class, useNA = "always"))

if (nrow(all_data) == 0) stop("No data loaded.")

# -- Plotting order --
present_labels   <- unique(as.character(all_data$molecule))
missing_in_order <- setdiff(present_labels, contig_order)
if (length(missing_in_order) > 0) {
  warning("Contigs not in contig_order, appended at bottom: ",
          paste(missing_in_order, collapse = ", "))
  contig_order <- c(contig_order, missing_in_order)
}
contig_order_present <- contig_order[contig_order %in% present_labels]
# Reverse so first entry appears at top of plot
all_data$molecule <- factor(all_data$molecule, levels = rev(contig_order_present))

# -- Legend labels (only for classes present in data) --
legend_labels <- c(
  resistance = "Resistance genes",
  transposon = "Transposon genes",
  integron   = "Integrons/integrase region",
  transcription     = "Transcriptional regulators",
  origin = "Origin of replication",
  other      = "Other"
)


# -- Plot --
p <- ggplot(all_data,
            aes(xmin    = start,
                xmax    = end,
                y       = molecule,
                fill    = class,
                label   = label_text,
                forward = strand_fwd)) +
  geom_gene_arrow(
    arrowhead_height  = unit(arrow_height,       "mm"),
    arrowhead_width   = unit(1.75,               "mm"),
    arrow_body_height = unit(arrow_height * 0.8, "mm")
  ) +
  geom_text(
    data = {
      lab <- subset(all_data, label_text != "")
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
  scale_fill_manual(
    values   = class_colours,
    labels   = legend_labels,
    name     = "Gene class",
    na.value = "#BBBBBB",
    drop     = TRUE
  ) +
  scale_x_continuous(labels = scales::comma,
                     expand = expansion(mult = 0.02)) +
  #labs(title = plot_title, x = "Position (bp)", y = NULL) +
  labs(x = "Position (bp)", y = NULL) +
  theme_genes() +
  theme(
    # plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    axis.text.y      = element_text(size = font_size_axis, face = "italic"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    legend.text      = element_text(size = 9),
    strip.background = element_blank(),
    strip.text       = element_blank()
  )

print(p)

# -- Save PDF --
n_contigs   <- length(levels(all_data$molecule))
output_file <- file.path(input_folder, "Bakta_contig_map.pdf")
ggsave(filename = output_file, plot = p,
       width  = output_width,
       height = output_height_per_contig * n_contigs + 1.5,
       units  = "in", device = "pdf")
message("Plot saved to: ", output_file)