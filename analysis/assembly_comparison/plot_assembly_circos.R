#!/usr/bin/env Rscript
# Circular ideogram per strain: reference genome as a single circular sector
# (matching the project's oriC/ter-normalized circular-chromosome framing),
# with SNP density, indel density, and structural-rearrangement tracks from
# NucDiff output. NucDiff structural calls are split into "genuine"
# rearrangement types (inversion, duplication, relocation, translocation,
# reshuffling) vs. assembly-boundary artifact types (unaligned_beginning/end,
# collapsed_repeat, circular_genome_start) -- only genuine types are drawn on
# the rearrangement track, since boundary artifacts are a property of contig
# fragmentation, not biology (see plot_assembly_dotplot.R's fragmentation
# caveat for the same distinction).
#
# Usage:
#   Rscript plot_assembly_circos.R <strain> <display_name> [comparison_dir] [ori_ter_json] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(RColorBrewer)
  library(circlize)
  library(sysfonts)
  library(showtext)
})

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 2) {
  stop("Usage: plot_assembly_circos.R <strain> <display_name> [comparison_dir] [ori_ter_json] [out_dir]")
}
strain         <- argv[1]
display_name   <- argv[2]
comparison_dir <- if (length(argv) >= 3) argv[3] else file.path("analysis/assembly_comparison", strain)
ori_ter_json   <- if (length(argv) >= 4) argv[4] else file.path("mapping/reference", strain, "ori_ter.json")
out_dir        <- if (length(argv) >= 5) argv[5] else file.path(comparison_dir, "figures")
individual_dir <- file.path(out_dir, "individual")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

dark2 <- brewer.pal(8, "Dark2")
col_snp    <- dark2[1]
col_indel  <- dark2[2]
rearr_types <- c("inversion", "duplication", "relocation-insertion", "relocation-overlap", "translocation")
rearr_colors <- setNames(colorRampPalette(dark2[3:8])(length(rearr_types)), rearr_types)
# "reshuffling-part_*_gr_0" calls are grouped here with the assembly-artifact
# types, not treated as genuine rearrangement: for a circular chromosome
# compared against another circular sequence with a different arbitrary FASTA
# start point, NucDiff (which has no notion of circularity) reports the
# resulting whole-genome offset as exactly two complementary "reshuffling"
# blocks spanning the full genome -- confirmed in this data (mff, AB42, AB30
# each show exactly reshuffling-part_1_gr_0 + reshuffling-part_2_gr_0, whose
# coordinates partition the entire chromosome). That is a rotation artifact of
# circular-sequence representation, not a real structural rearrangement.
artifact_types <- c("unaligned_beginning", "unaligned_end", "collapsed_repeat",
                     "collapsed_tandem_repeat", "circular_genome_start",
                     "reshuffling-part_1_gr_0", "reshuffling-part_2_gr_0")

# ---- inputs -------------------------------------------------------------
ori_ter <- fromJSON(ori_ter_json)
genome_length <- ori_ter$genome_length

snps_gff   <- file.path(comparison_dir, "nucdiff", "results", paste0(strain, "_ref_snps.gff"))
struct_gff <- file.path(comparison_dir, "nucdiff", "results", paste0(strain, "_ref_struct.gff"))
ani_tsv    <- file.path("analysis/assembly_comparison", "summary_ani.tsv")

read_gff <- function(path) {
  lines <- readLines(path)
  lines <- lines[!startsWith(lines, "#") & trimws(lines) != ""]
  if (length(lines) == 0) return(data.frame(start = numeric(0), end = numeric(0), type = character(0)))
  parts <- strsplit(lines, "\t")
  start <- as.numeric(sapply(parts, `[`, 4))
  end   <- as.numeric(sapply(parts, `[`, 5))
  attrs <- sapply(parts, `[`, 9)
  type  <- sub(".*Name=([^;]+).*", "\\1", attrs)
  data.frame(start = start, end = end, type = type)
}

snps <- read_gff(snps_gff)
struct <- read_gff(struct_gff)

ani_summary <- read_tsv(ani_tsv, show_col_types = FALSE)
ani_row <- ani_summary %>% filter(PairID == strain)
ani_value <- if (nrow(ani_row) == 1) as.numeric(ani_row$ANI[1]) else NA_real_
ani_label <- if (is.na(ani_value)) "NA" else sprintf("%.4f%%", ani_value)

# ---- fragmentation check (same threshold/logic as plot_assembly_dotplot.R) --
coords_path <- file.path(comparison_dir, "nucdiff", paste0(strain, ".coords"))
contig_lens <- {
  lines <- readLines(coords_path)
  lines <- lines[trimws(lines) != ""]
  toks_list <- lapply(lines, function(l) { t <- strsplit(trimws(l), "\\s+")[[1]]; t[t != "|"] })
  qtot <- as.numeric(sapply(toks_list, `[`, 9))
  qseq <- sapply(toks_list, function(t) t[length(t)])
  agg <- tapply(qtot, qseq, max)
  sort(agg, decreasing = TRUE)
}
n_contigs <- length(contig_lens)
largest_contig_frac <- contig_lens[1] / sum(contig_lens)
is_fragmented <- largest_contig_frac < 0.9

# ---- bin SNPs / indels into 10kb windows ---------------------------------
window <- 10000
n_windows <- ceiling(genome_length / window)
bin_of <- function(pos) pmin(floor(pos / window) + 1, n_windows)

snp_pos <- snps$start[snps$type == "substitution"]
indel_pos <- snps$start[snps$type %in% c("insertion", "deletion")]

snp_counts <- tabulate(bin_of(snp_pos), n_windows)
indel_counts <- tabulate(bin_of(indel_pos), n_windows)
bin_mid <- (seq_len(n_windows) - 0.5) * window

rearr <- struct %>% filter(type %in% rearr_types)
n_artifact <- sum(struct$type %in% artifact_types)

# ---- plot -----------------------------------------------------------------
title_txt <- paste0("Genome-wide variant landscape: ", display_name, " (FastANI = ", ani_label, ")")
rearr_summary <- if (nrow(rearr) > 0) paste(names(table(rearr$type)), table(rearr$type), sep = "=", collapse = ", ") else "none"
caption_txt <- paste0(
  "Strain ", display_name, ". Rings (outer to inner): genome position (Mb, oriC/ter marked), ",
  "substitution density per 10kb window (max ", max(snp_counts), "), insertion+deletion density per 10kb window (max ", max(indel_counts), "), and ",
  "structural rearrangements (", nrow(rearr), " shown: ", rearr_summary,
  "). ", n_artifact, " NucDiff structural record(s) of assembly-artifact type ",
  "(unaligned_beginning/end, collapsed_repeat(s), circular_genome_start, and whole-genome ",
  "reshuffling-part_1/2_gr_0 pairs -- the latter reflect the assembly's circular start point ",
  "differing from the reference's, not real rearrangement) are excluded from the rearrangement ring.",
  if (is_fragmented) paste0(
    " CAUTION: this hybrid assembly is fragmented (largest contig = ", sprintf("%.1f%%", largest_contig_frac * 100),
    " of total length across ", n_contigs, " contigs) -- SNP/indel density here is likely inflated by ",
    "assembly/sequencing error in low-confidence contigs, not confirmed biological variation."
  ) else ""
)

draw_circos <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 0, track.margin = c(0.005, 0.005), cell.padding = c(0, 0, 0, 0))
  circos.initialize(sectors = "genome", xlim = c(0, genome_length))

  # Track 1: position axis + oriC/ter markers (ylim has headroom above 1 so
  # the oriC/ter text labels don't fall outside the plotting region)
  circos.track(ylim = c(0, 1.6), track.height = 0.09, bg.border = NA, panel.fun = function(x, y) {
    circos.axis(h = 0, major.at = seq(0, genome_length, by = 5e5),
                labels = paste0(seq(0, genome_length, by = 5e5) / 1e6, " Mb"),
                labels.cex = 0.55, minor.ticks = 0)
  })
  circos.lines(rep(ori_ter$oriC, 2), c(0, 1), sector.index = "genome", track.index = 1, col = "grey30", lwd = 1.5)
  circos.lines(rep(ori_ter$ter, 2), c(0, 1), sector.index = "genome", track.index = 1, col = "grey55", lwd = 1.5)
  circos.text(ori_ter$oriC, 1.35, "oriC", sector.index = "genome", track.index = 1, cex = 0.75, family = FONT, col = "grey30",
              facing = "bending.outside", niceFacing = TRUE)
  circos.text(ori_ter$ter, 1.35, "ter", sector.index = "genome", track.index = 1, cex = 0.75, family = FONT, col = "grey45",
              facing = "bending.outside", niceFacing = TRUE)

  # Track 2: SNP density (no numeric axis -- max noted in caption to avoid clutter)
  max_snp <- max(1, max(snp_counts))
  circos.track(ylim = c(0, max_snp), track.height = 0.16, panel.fun = function(x, y) {
    circos.lines(bin_mid, snp_counts, type = "h", col = col_snp, lwd = 0.6, area = TRUE, border = col_snp)
  })

  # Track 3: indel density
  max_indel <- max(1, max(indel_counts))
  circos.track(ylim = c(0, max_indel), track.height = 0.16, panel.fun = function(x, y) {
    circos.lines(bin_mid, indel_counts, type = "h", col = col_indel, lwd = 0.6, area = TRUE, border = col_indel)
  })

  # Track 4: structural rearrangements (genuine types only)
  circos.track(ylim = c(0, 1), track.height = 0.1, bg.border = "grey85", panel.fun = function(x, y) {})
  if (nrow(rearr) > 0) {
    min_width <- genome_length * 0.0025
    for (i in seq_len(nrow(rearr))) {
      s <- rearr$start[i]; e <- max(rearr$end[i], rearr$start[i] + min_width)
      circos.rect(s, 0, e, 1, sector.index = "genome", track.index = 4,
                  col = rearr_colors[[rearr$type[i]]], border = NA)
    }
  }
}

save_plot <- function(path_fn) {
  path_fn()
  # Reserve a wide bottom margin (outside the circular plotting region) for
  # the wrapped caption, so it can never overlap the circle or the legends.
  cap_lines <- strwrap(caption_txt, width = 108)
  par(mar = c(length(cap_lines) + 1.5, 1, 3, 1), family = FONT)
  draw_circos()
  title(main = title_txt, cex.main = 1.05, family = FONT)
  # Legends sit in the empty corners above the circle and below the title --
  # confirmed clear of both the plot and the caption margin.
  legend_types <- unique(rearr$type)
  if (length(legend_types) > 0) {
    legend("topleft", legend = legend_types, fill = rearr_colors[legend_types],
           cex = 0.65, bty = "n", title = "Rearrangement type", text.font = 1, inset = c(0, 0.06))
  }
  legend("topright", legend = c("Substitution density", "Indel density"),
         fill = c(col_snp, col_indel), cex = 0.65, bty = "n", title = "Density tracks", inset = c(0, 0.06))
  for (i in seq_along(cap_lines)) {
    mtext(cap_lines[i], side = 1, line = i, adj = 0, cex = 0.55, family = FONT, col = "grey20")
  }
  dev.off()
}

out_name <- paste0(strain, "_circos")
save_plot(function() png(file.path(individual_dir, paste0(out_name, ".png")), width = 8, height = 9, units = "in", res = 300, bg = "white"))
save_plot(function() cairo_pdf(file.path(individual_dir, paste0(out_name, ".pdf")), width = 8, height = 9))
save_plot(function() png(file.path(out_dir, paste0(out_name, ".png")), width = 8, height = 9, units = "in", res = 300, bg = "white"))
save_plot(function() cairo_pdf(file.path(out_dir, paste0(out_name, ".pdf")), width = 8, height = 9))

message("Saved: ", out_name, ".png / .pdf to ", out_dir, " (and individual/)")
