#!/usr/bin/env Rscript
# Whole-genome synteny dot plot: hybrid assembly (Unicycler) vs. public
# reference genome, per strain. Segments are MUMmer/NucDiff alignment blocks
# (one row per aligned region, from <strain>.coords), reference position on x,
# hybrid-assembly position on y (contigs stacked with a cumulative offset so
# each query contig gets its own y-range, standard practice for multi-contig
# whole-genome alignment plots). Colored by alignment orientation -- large
# off-diagonal or anti-diagonal segments reveal relocations/inversions.
#
# Usage:
#   Rscript plot_assembly_dotplot.R <strain> <display_name> [comparison_dir] [ori_ter_json] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(RColorBrewer)
  library(jsonlite)
  library(sysfonts)
  library(showtext)
})

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 2) {
  stop("Usage: plot_assembly_dotplot.R <strain> <display_name> [comparison_dir] [ori_ter_json] [out_dir]")
}
strain         <- argv[1]
display_name   <- argv[2]
comparison_dir <- if (length(argv) >= 3) argv[3] else file.path("analysis/assembly_comparison", strain)
ori_ter_json   <- if (length(argv) >= 4) argv[4] else file.path("mapping/reference", strain, "ori_ter.json")
out_dir        <- if (length(argv) >= 5) argv[5] else file.path(comparison_dir, "figures")
individual_dir <- file.path(out_dir, "individual")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

# ---- fonts / colors (matches analysis/replication_profile conventions) -----
font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

dark2 <- brewer.pal(3, "Dark2")
palette_strand <- c(Forward = dark2[1], Reverse = dark2[2])

base_theme <- theme_bw(base_size = 13, base_family = FONT) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12, family = FONT),
    plot.subtitle = element_text(size = 10, color = "grey30", family = FONT),
    plot.caption = element_text(size = 8.5, family = FONT, hjust = 0, color = "grey20"),
    plot.caption.position = "plot",
    legend.position = "bottom",
    text = element_text(family = FONT)
  )

# ---- inputs -----------------------------------------------------------------
coords_path <- file.path(comparison_dir, "nucdiff", paste0(strain, ".coords"))
ani_tsv     <- file.path("analysis/assembly_comparison", "summary_ani.tsv")

ori_ter <- fromJSON(ori_ter_json)

ani_summary <- read_tsv(ani_tsv, show_col_types = FALSE)
ani_row <- ani_summary %>% filter(PairID == strain)
ani_value <- if (nrow(ani_row) == 1) as.numeric(ani_row$ANI[1]) else NA_real_

# MUMmer/NucDiff show-coords output: space-delimited with literal "|" column
# separators and no header; last two fields (ref_id, query_seq_num) are
# tab-separated. Strip "|" tokens and parse positionally.
read_coords <- function(path) {
  lines <- readLines(path)
  lines <- lines[trimws(lines) != ""]
  parsed <- lapply(lines, function(l) {
    toks <- strsplit(trimws(l), "\\s+")[[1]]
    toks[toks != "|"]
  })
  n <- length(parsed[[1]])
  df <- as.data.frame(do.call(rbind, parsed), stringsAsFactors = FALSE)
  colnames(df) <- c("ref_start", "ref_end", "query_start", "query_end",
                     "ref_len", "query_len", "identity",
                     "ref_total_len", "query_total_len",
                     "ref_cov", "query_cov", "ref_frame", "query_frame",
                     "ref_id", "query_seq_num")[seq_len(n)]
  num_cols <- c("ref_start", "ref_end", "query_start", "query_end", "ref_len",
                "query_len", "identity", "ref_total_len", "query_total_len",
                "ref_cov", "query_cov", "ref_frame", "query_frame")
  df %>% mutate(across(all_of(num_cols), as.numeric))
}

coords <- read_coords(coords_path)

# ---- stack query contigs along y with a cumulative offset -------------------
contig_lengths <- coords %>%
  group_by(query_seq_num) %>%
  summarise(len = max(query_total_len), .groups = "drop") %>%
  arrange(desc(len)) %>%
  mutate(offset = cumsum(lag(len, default = 0)) + (row_number() - 1) * (max(len) * 0.02))

coords <- coords %>%
  left_join(contig_lengths %>% select(query_seq_num, offset), by = "query_seq_num") %>%
  mutate(
    strand = ifelse(query_frame >= 0, "Forward", "Reverse"),
    y_start = (query_start + offset) / 1e6,
    y_end   = (query_end + offset) / 1e6
  )

boundaries <- contig_lengths %>%
  mutate(y_top = (offset + len) / 1e6, y_label = offset / 1e6 + (len / 1e6) / 2)

n_contigs <- nrow(contig_lengths)
main_contig_len <- contig_lengths$len[1]
total_len <- sum(contig_lengths$len)
largest_contig_frac <- main_contig_len / total_len
is_fragmented <- largest_contig_frac < 0.9

# ---- caption / title ---------------------------------------------------------
n_blocks <- nrow(coords)
mean_identity <- weighted.mean(coords$identity, coords$ref_len)
ani_label <- if (is.na(ani_value)) "NA" else sprintf("%.4f%%", ani_value)

title_txt <- paste0("Whole-genome synteny: ", display_name, " hybrid assembly vs. reference (FastANI = ", ani_label, ")")
caption_txt <- paste0(
  "Strain ", display_name, ". Each segment is one MUMmer/NucDiff alignment block between the ",
  "reference genome (x, Mb) and the Unicycler hybrid assembly (y, Mb); ", n_contigs,
  " assembly contig(s) are stacked with a cumulative offset (grey dashed boundaries, largest first). ",
  n_blocks, " alignment blocks shown, mean %identity (block-length-weighted) = ", sprintf("%.2f%%", mean_identity),
  ". Colored by alignment orientation: off-diagonal blocks indicate relocated/reshuffled sequence; ",
  "anti-diagonal (Reverse) blocks indicate inversions. Dashed vertical lines mark oriC (position 0) and ter.",
  if (is_fragmented) paste0(
    " CAUTION: this hybrid assembly is fragmented (largest contig = ", sprintf("%.1f%%", largest_contig_frac * 100),
    " of total length across ", n_contigs, " contigs) -- most structural differences from this comparison are ",
    "likely assembly artifacts (contig-boundary effects) rather than confirmed biological rearrangement."
  ) else ""
)

# ---- plot ---------------------------------------------------------------------
p <- ggplot(coords) +
  geom_vline(xintercept = ori_ter$oriC / 1e6, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = ori_ter$ter / 1e6, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_hline(data = boundaries, aes(yintercept = y_top), linetype = "dashed", color = "grey75", linewidth = 0.3) +
  geom_segment(aes(x = ref_start / 1e6, xend = ref_end / 1e6, y = y_start, yend = y_end, color = strand),
               alpha = 0.85, linewidth = 0.9, lineend = "round") +
  annotate("text", x = ori_ter$oriC / 1e6, y = Inf, label = "oriC", vjust = 1.3, hjust = -0.1,
           size = 3.0, color = "grey35", family = FONT) +
  annotate("text", x = ori_ter$ter / 1e6, y = Inf, label = "ter", vjust = 1.3, hjust = -0.1,
           size = 3.0, color = "grey45", family = FONT) +
  scale_color_manual(values = palette_strand, name = "Alignment orientation") +
  scale_x_continuous(labels = label_number(suffix = " Mb")) +
  scale_y_continuous(labels = label_number(suffix = " Mb")) +
  labs(x = "Reference genome position", y = "Hybrid assembly position", title = title_txt) +
  base_theme

if (is_fragmented) {
  warning_label <- sprintf("FRAGMENTED ASSEMBLY (largest contig %.0f%% of total, %d contigs)\nStructural calls not confirmed biological rearrangement",
                            largest_contig_frac * 100, n_contigs)
  p <- p + annotate("label", x = -Inf, y = -Inf, label = warning_label, hjust = -0.02, vjust = -0.3,
                     size = 3.1, family = FONT, fontface = "bold", color = "#B22222",
                     fill = alpha("#FFF3F3", 0.9))
}

p_final <- p + labs(caption = paste(strwrap(caption_txt, width = 130), collapse = "\n"))

out_name <- paste0(strain, "_dotplot")
ggsave(file.path(individual_dir, paste0(out_name, ".png")), p_final, width = 8.5, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(individual_dir, paste0(out_name, ".pdf")), p_final, width = 8.5, height = 7.5, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(out_name, ".png")), p_final, width = 8.5, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(out_name, ".pdf")), p_final, width = 8.5, height = 7.5, device = cairo_pdf)

message("Saved: ", out_name, ".png / .pdf to ", out_dir, " (and individual/)")
