#!/usr/bin/env Rscript
# Publication-ready marker-frequency / replication-profile plots, following
# Skovgaard et al. 2011 (Genome Research) methodology: x = genomic position
# normalized per replichore (oriC = 0, ter = +-1), y = log2(marker frequency)
# relative to the genome-wide median window count. Exponential-phase cultures
# are expected to show a "tent" peak at oriC tapering to ter; stationary-phase
# cultures (non-replicating) should be flat.
#
# Points are plotted independently at both 1kb and 10kb resolution; no line
# connects them, since adjacent genomic windows are not literally linked
# observations. Each panel is annotated with an Ori:Ter ratio (median marker
# frequency near oriC divided by median marker frequency near ter, from the
# 10kb series), the same style of summary statistic reported in the reference
# paper. The y-axis uses the same fixed range on every panel, across every
# strain, so figures are directly comparable.
#
# Usage:
#   Rscript plot_replication_profile.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]
#
# strain        short code used in file/column names, e.g. mff
# display_name  full strain name for titles/captions, e.g. "ATCC 17978-mff"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(RColorBrewer)
  library(patchwork)
  library(sysfonts)
  library(showtext)
})

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 2) {
  stop("Usage: plot_replication_profile.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]")
}
strain        <- argv[1]
display_name  <- argv[2]
coverage_dir  <- if (length(argv) >= 3) argv[3] else file.path("mapping/coverage", strain)
metadata_csv  <- if (length(argv) >= 4) argv[4] else "mapping/coverage/sample_metadata.csv"
out_dir       <- if (length(argv) >= 5) argv[5] else file.path("analysis/replication_profile", strain, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- fonts --------------------------------------------------------------
font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

# ---- colors (colorblind-safe, RColorBrewer Dark2) ------------------------
dark2 <- brewer.pal(3, "Dark2")
palette_replicate <- c(A = dark2[1], B = dark2[2])
palette_treatment <- c(drug = dark2[1], ND = dark2[2])

# ---- fixed y-axis, shared across every strain and phase for comparability --
Y_LIMITS <- c(-5, 3)
Y_BREAKS <- seq(-5, 3, 1)
LABEL_Y  <- 2.8

# ---- data ------------------------------------------------------------------
metadata <- read_csv(metadata_csv, show_col_types = FALSE) %>%
  filter(strain == !!strain)

read_sample_tsv <- function(pid_sample, window_label) {
  path <- file.path(coverage_dir, paste0(pid_sample, "_", window_label, ".tsv"))
  df <- read_tsv(path, show_col_types = FALSE)
  df$pid_sample <- pid_sample
  df
}

load_resolution <- function(window_label) {
  do.call(rbind, lapply(metadata$pid_sample, read_sample_tsv, window_label = window_label)) %>%
    left_join(metadata, by = "pid_sample") %>%
    mutate(abs_pos = abs(norm_position))
}

data_10kb <- load_resolution("10kb")
data_1kb  <- load_resolution("1kb")

average_replicates <- function(df, trt, phs) {
  sub <- df %>% filter(treatment == trt, phase == phs)
  sub %>%
    group_by(window_mid, norm_position, abs_pos) %>%
    summarise(log2_ratio = mean(log2_ratio), .groups = "drop") %>%
    mutate(treatment = trt, phase = phs)
}

# ---- Ori:Ter ratio ----------------------------------------------------------
# Median log2(marker frequency) in windows near oriC (|position| <= 0.1) versus
# near ter (|position| >= 0.9), converted back to a linear-scale fold ratio.
# Median (not mean) so a handful of extreme windows near either boundary don't
# dominate the estimate.
ORI_ZONE <- 0.1
ratio_records <- list()

compute_ratio <- function(df) {
  ori_val <- median(df$log2_ratio[abs(df$norm_position) <= ORI_ZONE])
  ter_val <- median(df$log2_ratio[abs(df$norm_position) >= (1 - ORI_ZONE)])
  2^(ori_val - ter_val)
}

format_ratio <- function(label, ratio) {
  sprintf("%s  Ori:Ter = %.2f", label, ratio)
}

record_ratio <- function(phs, panel, series, ratio) {
  ratio_records[[length(ratio_records) + 1]] <<- data.frame(
    phase = phs, panel = panel, series = series, ori_ter_ratio = ratio
  )
}

# ---- shared plot elements --------------------------------------------------
base_theme <- theme_bw(base_size = 13, base_family = FONT) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, family = FONT),
    plot.subtitle = element_text(size = 10, color = "grey30", family = FONT),
    legend.position = "bottom",
    text = element_text(family = FONT)
  )

ori_ter_layers <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4),
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.4),
  annotate("text", x = 0.03, y = LABEL_Y, label = "oriC", hjust = 0,
           size = 3.2, color = "grey35", family = FONT),
  annotate("text", x = c(-0.97, 0.97), y = LABEL_Y, label = c("ter", "ter"),
           hjust = c(0, 1), size = 3.2, color = "grey55", family = FONT)
)

axis_layers <- list(
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)),
  scale_y_continuous(limits = Y_LIMITS, breaks = Y_BREAKS),
  labs(x = "Normalized position relative to oriC (oriC = 0, ter = ±1)",
       y = expression(log[2]*"(marker frequency)"))
)

ratio_annotation <- function(label) {
  annotate("label", x = Inf, y = Inf, label = label, hjust = 1.02, vjust = 1.3,
           size = 3.1, family = FONT, fill = alpha("white", 0.75))
}

# ---- panel builders ---------------------------------------------------------
panel_AB <- function(trt, phs, title) {
  sub10 <- data_10kb %>% filter(treatment == trt, phase == phs)
  sub1  <- data_1kb  %>% filter(treatment == trt, phase == phs)

  ratio_a <- compute_ratio(sub10 %>% filter(replicate == "A"))
  ratio_b <- compute_ratio(sub10 %>% filter(replicate == "B"))
  record_ratio(phs, title, paste0(trt, "_A"), ratio_a)
  record_ratio(phs, title, paste0(trt, "_B"), ratio_b)
  ratio_label <- paste(format_ratio("A:", ratio_a), format_ratio("B:", ratio_b), sep = "\n")

  ggplot() +
    geom_point(data = sub1, aes(x = norm_position, y = log2_ratio, color = replicate),
               alpha = 0.12, size = 0.5, show.legend = FALSE) +
    geom_point(data = sub10, aes(x = norm_position, y = log2_ratio, color = replicate),
               alpha = 0.8, size = 1.3) +
    ori_ter_layers + axis_layers +
    scale_color_manual(values = palette_replicate, name = "Biological replicate") +
    ratio_annotation(ratio_label) +
    labs(title = title) +
    base_theme
}

panel_comparison <- function(phs, title) {
  sub10 <- rbind(average_replicates(data_10kb, "drug", phs), average_replicates(data_10kb, "ND", phs))
  sub1  <- rbind(average_replicates(data_1kb, "drug", phs), average_replicates(data_1kb, "ND", phs))

  ratio_drug <- compute_ratio(sub10 %>% filter(treatment == "drug"))
  ratio_nd   <- compute_ratio(sub10 %>% filter(treatment == "ND"))
  record_ratio(phs, title, "drug_avg", ratio_drug)
  record_ratio(phs, title, "ND_avg", ratio_nd)
  ratio_label <- paste(format_ratio("Drug:", ratio_drug), format_ratio("Untreated:", ratio_nd), sep = "\n")

  ggplot() +
    geom_point(data = sub1, aes(x = norm_position, y = log2_ratio, color = treatment),
               alpha = 0.12, size = 0.5, show.legend = FALSE) +
    geom_point(data = sub10, aes(x = norm_position, y = log2_ratio, color = treatment),
               alpha = 0.8, size = 1.3) +
    ori_ter_layers + axis_layers +
    scale_color_manual(values = palette_treatment, name = "Treatment",
                        labels = c(drug = "Drug-treated", ND = "Untreated (no drug)")) +
    ratio_annotation(ratio_label) +
    labs(title = title) +
    base_theme
}

# ---- build one composite figure per growth phase ---------------------------
build_phase_figure <- function(phs, phase_label, out_name) {
  pA <- panel_AB("drug", phs, "Drug-treated: replicate A vs B")
  pB <- panel_AB("ND",   phs, "Untreated (no drug): replicate A vs B")
  pC <- panel_comparison(phs, "Drug-treated vs untreated (replicate-averaged)")

  combined <- (pA | pB | pC) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  caption <- paste0(
    "Strain ", display_name, ", ", phase_label, " phase. oriC = 0, ter = ±1. ",
    "Ori:Ter is the median marker-frequency ratio between windows near oriC and near ter."
  )

  final <- combined +
    plot_annotation(
      title = paste0("Replication Profile: Acinetobacter baumannii Strain ", display_name, ", ", phase_label, " Phase"),
      caption = caption,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, family = FONT, hjust = 0),
        plot.caption = element_text(size = 10, family = FONT, hjust = 0, color = "grey20"),
        plot.caption.position = "plot"
      )
    )

  ggsave(file.path(out_dir, paste0(out_name, ".png")), final, width = 15, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(out_name, ".pdf")), final, width = 15, height = 5.5, device = cairo_pdf)
  message("Saved: ", out_name, ".png / .pdf")
}

build_phase_figure("stat", "Stationary", "stationary_phase")
build_phase_figure("exp",  "Exponential", "exponential_phase")

ratio_table <- do.call(rbind, ratio_records)
write_csv(ratio_table, file.path(out_dir, "..", "ori_ter_ratios.csv"))

message("Both composite figures written to ", out_dir)
message("Ori:Ter ratio table written to ", file.path(out_dir, "..", "ori_ter_ratios.csv"))
