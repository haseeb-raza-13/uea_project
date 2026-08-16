#!/usr/bin/env Rscript
# Cross-strain summary: 3-panel composite from analysis/assembly_comparison/summary_table.csv.
# Panel A: variant-type counts per strain (log10 scale -- substitution counts
# span ~60 to ~11,000 across strains, a linear scale would make the small
# strains invisible). Panel B: FastANI %, lollipop style (a truncated-baseline
# bar chart would misrepresent magnitude for values this tightly clustered
# near 100%). Panel C: assembly completeness (largest contig as % of total
# assembly length, contig count annotated per bar) -- the evidence behind the
# fragmentation caveat on the AB30/Lac4 dot plots and circos ideograms.
#
# Usage: Rscript plot_assembly_summary.R [summary_csv] [out_dir]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(RColorBrewer)
  library(patchwork)
  library(sysfonts)
  library(showtext)
})

argv <- commandArgs(trailingOnly = TRUE)
summary_csv <- if (length(argv) >= 1) argv[1] else "analysis/assembly_comparison/summary_table.csv"
out_dir     <- if (length(argv) >= 2) argv[2] else "analysis/assembly_comparison/figures_summary"
individual_dir <- file.path(out_dir, "individual")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

dark2 <- brewer.pal(8, "Dark2")

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

display_names <- c(mff = "ATCC 17978-mff", VU = "ATCC 17978-VU", AB30 = "AB030", AB42 = "AB042", Lac4 = "LAC-4")
strain_order <- c("mff", "VU", "AB42", "AB30", "Lac4")  # ascending fragmentation, for readability

d <- read_csv(summary_csv, show_col_types = FALSE) %>%
  mutate(display = factor(display_names[strain], levels = display_names[strain_order]))

# ---- Panel A: variant-type counts (log10 scale) -----------------------------
variant_cols <- c(substitutions = "Substitutions", insertions = "Insertions", deletions = "Deletions",
                   inversions = "Inversions", relocations = "Relocations", reshufflings = "Reshufflings")

d_long <- d %>%
  select(display, all_of(names(variant_cols))) %>%
  pivot_longer(-display, names_to = "category", values_to = "count") %>%
  mutate(category = factor(variant_cols[category], levels = unname(variant_cols)))

panel_a <- ggplot(d_long, aes(x = display, y = count + 1, fill = category)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.75), vjust = -0.4,
            size = 2.3, family = FONT, angle = 90, hjust = -0.05) +
  scale_y_log10(labels = label_number(), expand = expansion(mult = c(0, 0.28))) +
  scale_fill_brewer(palette = "Dark2", name = "Variant type") +
  labs(x = NULL, y = "Count + 1 (log10 scale)",
       title = "A. NucDiff variant-type counts by strain") +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ---- Panel B: FastANI %, lollipop (tight range near 100%) -------------------
ani_min <- floor(min(d$ani_percent) * 100) / 100 - 0.02
panel_b <- ggplot(d, aes(x = display, y = ani_percent)) +
  geom_segment(aes(xend = display, y = ani_min, yend = ani_percent), color = "grey60", linewidth = 0.6) +
  geom_point(size = 3.5, color = dark2[1]) +
  geom_text(aes(label = sprintf("%.4f%%", ani_percent)), vjust = -1.1, size = 3.0, family = FONT) +
  coord_cartesian(ylim = c(ani_min, 100)) +
  labs(x = NULL, y = "FastANI (%)", title = "B. Average Nucleotide Identity by strain") +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ---- Panel C: assembly completeness ------------------------------------------
d <- d %>% mutate(frac_pct = largest_contig_frac * 100,
                   fragmented = frac_pct < 90)
panel_c <- ggplot(d, aes(x = display, y = frac_pct, fill = fragmented)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text(aes(label = paste0("n=", n_contigs, " contigs")), vjust = -0.6, size = 2.9, family = FONT) +
  scale_fill_manual(values = c(`FALSE` = dark2[3], `TRUE` = "#B22222"), guide = "none") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = "Largest contig (% of total assembly length)",
       title = "C. Hybrid assembly completeness by strain",
       subtitle = "Dashed line: 90% threshold used to flag fragmented assemblies in other figures") +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_individual <- function(p, out_name, w = 6.5, h = 5.5) {
  ggsave(file.path(individual_dir, paste0(out_name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(individual_dir, paste0(out_name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
}
save_individual(panel_a, "summary_panelA_variant_counts", w = 7.5)
save_individual(panel_b, "summary_panelB_ani")
save_individual(panel_c, "summary_panelC_completeness")

caption_txt <- paste0(
  "Cross-strain summary of hybrid assembly (Unicycler) vs. public reference genome comparisons (FastANI + NucDiff), ",
  "all 5 strains. Panel A: NucDiff variant-type counts (translocations omitted, 0 in all strains; reshuffling counts ",
  "reflect whole-genome start-point offsets between circular sequences, not confirmed internal rearrangement -- see ",
  "per-strain dot plots/ideograms). Panel B: FastANI average nucleotide identity. Panel C: assembly completeness ",
  "(largest contig as % of total assembly length); AB30 and Lac4 fall well below the 90% completeness threshold, ",
  "meaning their elevated variant counts in Panel A are substantially assembly-fragmentation artifacts rather than ",
  "confirmed biological variation."
)

combined <- (panel_a / (panel_b | panel_c)) +
  plot_annotation(
    title = "Hybrid Assembly vs. Reference Genome: Cross-Strain Comparison Summary",
    caption = paste(strwrap(caption_txt, width = 170), collapse = "\n"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, family = FONT, hjust = 0),
      plot.caption = element_text(size = 9.5, family = FONT, hjust = 0, color = "grey20"),
      plot.caption.position = "plot"
    )
  )

ggsave(file.path(out_dir, "assembly_comparison_summary.png"), combined, width = 11, height = 11, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "assembly_comparison_summary.pdf"), combined, width = 11, height = 11, device = cairo_pdf)

message("Saved: assembly_comparison_summary.png / .pdf to ", out_dir, " (+ 3 individual panels)")
