#!/usr/bin/env Rscript
# Zoomed breakpoint diagrams for AB30's 2 inversion calls -- the only
# inversions found across all 5 strains (see analysis/assembly_comparison/summary_table.csv).
# Each panel is a small synteny ribbon: reference segment (top, forward) vs.
# hybrid-assembly contig segment (bottom, reverse orientation), with the
# crossed connecting ribbon being the standard visual signature of an
# inversion. Built as a case-study "hero figure", but captioned with the
# evidence for caution: both events sit on a 107kb contig (not the 1.5Mb main
# chromosome scaffold) in a mid-sized fragmented assembly (see
# plot_assembly_dotplot.R / plot_assembly_circos.R caveats for AB30), local
# alignment identity here (~97.4-97.9%) is below the genome-wide mean
# (99.76%), and the two inversion blocks map to opposite ends of that same
# contig -- a pattern consistent with a contig-end/assembly-junction artifact,
# not necessarily a confirmed biological inversion. Recommend long-read
# breakpoint validation before treating as confirmed.
#
# Usage: Rscript plot_assembly_locus.R [strain] [display_name] [comparison_dir] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(sysfonts)
  library(showtext)
})

argv <- commandArgs(trailingOnly = TRUE)
strain         <- if (length(argv) >= 1) argv[1] else "AB30"
display_name   <- if (length(argv) >= 2) argv[2] else "AB030"
comparison_dir <- if (length(argv) >= 3) argv[3] else file.path("analysis/assembly_comparison", strain)
out_dir        <- if (length(argv) >= 4) argv[4] else "analysis/assembly_comparison/figures_locus"
individual_dir <- file.path(out_dir, "individual")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

dark2 <- brewer.pal(8, "Dark2")
inv_color <- dark2[4]

base_theme <- theme_bw(base_size = 13, base_family = FONT) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", size = 11, family = FONT),
    plot.subtitle = element_text(size = 9, color = "grey30", family = FONT),
    plot.caption = element_text(size = 8.5, family = FONT, hjust = 0, color = "grey20"),
    plot.caption.position = "plot",
    text = element_text(family = FONT)
  )

# ---- the 2 inversion events, read directly from NucDiff's struct GFF --------
struct_gff <- file.path(comparison_dir, "nucdiff", "results", paste0(strain, "_ref_struct.gff"))
lines <- readLines(struct_gff)
lines <- lines[!startsWith(lines, "#") & grepl("Name=inversion", lines)]
parts <- strsplit(lines, "\t")

parse_attr <- function(attrs, key) sub(paste0(".*", key, "=([^;]+).*"), "\\1", attrs)

inv <- do.call(rbind, lapply(parts, function(p) {
  attrs <- p[9]
  data.frame(
    ref_start = as.numeric(p[4]), ref_end = as.numeric(p[5]),
    query_coord = parse_attr(attrs, "query_coord"),
    query_seq = parse_attr(attrs, "query_sequence"),
    id = parse_attr(attrs, "ID")
  )
}))
inv$query_start <- as.numeric(sub("-.*", "", inv$query_coord))
inv$query_end   <- as.numeric(sub(".*-", "", inv$query_coord))

# ---- ribbon diagram builder --------------------------------------------------
ribbon_panel <- function(row, panel_label) {
  ref_span <- row$ref_end - row$ref_start
  q_lo <- min(row$query_start, row$query_end); q_hi <- max(row$query_start, row$query_end)

  ref_df <- data.frame(x = c(row$ref_start, row$ref_end), y = c(1, 1))
  # query is reverse-oriented (inversion): draw with start/end swapped left-right
  qry_df <- data.frame(x = c(row$ref_start, row$ref_end), y = c(0, 0))  # positioned under ref span for the ribbon
  ribbon_df <- data.frame(
    x = c(row$ref_start, row$ref_end, row$ref_start, row$ref_end),
    y = c(1, 1, 0, 0),
    grp = 1
  )
  # crossed connector: top-left -> bottom-right, top-right -> bottom-left
  cross_df <- data.frame(
    x = c(row$ref_start, row$ref_end, row$ref_end, row$ref_start),
    y = c(1, 1, 0, 0)
  )

  ggplot() +
    geom_polygon(data = cross_df, aes(x = x, y = y), fill = inv_color, alpha = 0.25) +
    geom_segment(aes(x = row$ref_start, xend = row$ref_end, y = 1, yend = 1),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.3, color = "grey20") +
    geom_segment(aes(x = row$ref_end, xend = row$ref_start, y = 0, yend = 0),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.3, color = inv_color) +
    annotate("text", x = mean(c(row$ref_start, row$ref_end)), y = 1.35,
             label = sprintf("Reference: %s bp (forward)", format(ref_span, big.mark = ",")),
             size = 3.0, family = FONT, color = "grey20") +
    annotate("text", x = mean(c(row$ref_start, row$ref_end)), y = -0.35,
             label = sprintf("Contig %s: %s-%s bp (reverse)", row$query_seq,
                              format(q_lo, big.mark = ","), format(q_hi, big.mark = ",")),
             size = 3.0, family = FONT, color = inv_color) +
    scale_x_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = 0.15)) +
    scale_y_continuous(limits = c(-0.6, 1.6)) +
    labs(x = paste0("Reference position (bp), ", row$id), y = NULL, title = panel_label) +
    base_theme
}

p1 <- ribbon_panel(inv[1, ], sprintf("Inversion 1: ref %s-%s bp", format(inv$ref_start[1], big.mark=","), format(inv$ref_end[1], big.mark=",")))
p2 <- ribbon_panel(inv[2, ], sprintf("Inversion 2: ref %s-%s bp", format(inv$ref_start[2], big.mark=","), format(inv$ref_end[2], big.mark=",")))

caption_txt <- paste0(
  "Strain ", display_name, ". Both NucDiff-called inversions sit on the same 106,967bp assembly contig ",
  "(contig 5 of 175), not the 1,515,311bp main chromosome scaffold -- this assembly's largest contig is only ",
  "29.0% of total assembly length (see Panel C completeness figure), i.e. it is fragmented. Local alignment identity ",
  "at these breakpoints (97.4% and 97.9%) is below the genome-wide mean for this comparison (99.76%). The two ",
  "inversion blocks map to opposite ends of the same 107kb contig (positions ~2 and ~106,900) onto two reference ",
  "windows only ~460bp apart -- a pattern consistent with a contig-end/assembly-junction artifact rather than a ",
  "confirmed internal biological inversion. CAUTION: treat as a candidate finding pending confirmation via raw ",
  "long-read breakpoint support (read alignment spanning the junction), not as a confirmed rearrangement."
)

combined <- (p1 | p2) +
  plot_annotation(
    title = paste0("Candidate Inversion Breakpoints: ", display_name, " Hybrid Assembly vs. Reference"),
    caption = paste(strwrap(caption_txt, width = 165), collapse = "\n"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, family = FONT, hjust = 0),
      plot.caption = element_text(size = 9, family = FONT, hjust = 0, color = "grey20"),
      plot.caption.position = "plot"
    )
  )

save_individual <- function(p, out_name) {
  ggsave(file.path(individual_dir, paste0(out_name, ".png")), p, width = 5.5, height = 4.5, dpi = 300, bg = "white")
  ggsave(file.path(individual_dir, paste0(out_name, ".pdf")), p, width = 5.5, height = 4.5, device = cairo_pdf)
}
save_individual(p1, paste0(strain, "_inversion1"))
save_individual(p2, paste0(strain, "_inversion2"))

out_name <- paste0(strain, "_inversion_loci")
ggsave(file.path(out_dir, paste0(out_name, ".png")), combined, width = 11, height = 6, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(out_name, ".pdf")), combined, width = 11, height = 6, device = cairo_pdf)

message("Saved: ", out_name, ".png / .pdf to ", out_dir, " (+ 2 individual panels)")
