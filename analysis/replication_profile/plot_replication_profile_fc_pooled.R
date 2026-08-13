#!/usr/bin/env Rscript
# Fold-change variant of Panel C only (Panels A and B are identical to the
# log2 version, plot_replication_profile.R). Instead of plotting Meropenem and
# untreated read density side by side, Panel C here plots, per 1kb genomic
# window, log2(Meropenem density / untreated density) -- "dividing Meropenem
# by control, then log2 of that ratio" -- against normalized genomic position.
#
# Statistical design ("Design 1: pooled/overall test"): with only 2 biological
# replicates per treatment, a per-window t-test/Mann-Whitney U is essentially
# meaningless (Mann-Whitney U can never reach p < 0.05 at n=2 vs n=2). Instead,
# ALL windows for a strain/phase are pooled as the sample (n in the hundreds to
# thousands), and ONE overall t-test and ONE overall Mann-Whitney U test compare
# the full set of window-level mean densities, Meropenem vs untreated, answering
# "does Meropenem shift read density overall, genome-wide?" That single p-value
# pair is annotated on the panel, alongside vertical fold-change threshold lines
# at log2FC = -1 and +1. See plot_replication_profile_fc_volcano.R for the
# classic per-window volcano variant (Design 2).
#
# Usage:
#   Rscript plot_replication_profile_fc_pooled.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]

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
  stop("Usage: plot_replication_profile_fc_pooled.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]")
}
strain        <- argv[1]
display_name  <- argv[2]
coverage_dir  <- if (length(argv) >= 3) argv[3] else file.path("mapping/coverage", strain)
metadata_csv  <- if (length(argv) >= 4) argv[4] else "mapping/coverage/sample_metadata.csv"
out_dir       <- if (length(argv) >= 5) argv[5] else file.path("analysis/replication_profile", strain, "figures_fc_pooled")
individual_dir <- file.path(out_dir, "individual")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

# ---- fonts --------------------------------------------------------------
font_add_google("Lato", "Lato")
showtext_auto()
showtext_opts(dpi = 300)
FONT <- "Lato"

# ---- colors (colorblind-safe, RColorBrewer Dark2) ------------------------
dark2 <- brewer.pal(3, "Dark2")
palette_replicate <- c(A = dark2[1], B = dark2[2])
palette_fc <- c("Enriched in Meropenem" = dark2[1],
                "Depleted in Meropenem" = dark2[2],
                "No substantial change" = "grey70")

# ---- data (identical to plot_replication_profile.R, for Panels A/B) --------
metadata <- read_csv(metadata_csv, show_col_types = FALSE) %>%
  filter(strain == !!strain)

read_sample_tsv <- function(pid_sample, window_label) {
  path <- file.path(coverage_dir, paste0(pid_sample, "_", window_label, ".tsv"))
  df <- read_tsv(path, show_col_types = FALSE)
  df$pid_sample <- pid_sample
  df
}

load_resolution <- function(window_label, window_kb) {
  do.call(rbind, lapply(metadata$pid_sample, read_sample_tsv, window_label = window_label)) %>%
    left_join(metadata, by = "pid_sample") %>%
    mutate(abs_pos = abs(norm_position),
           log2_depth = log2(raw_count / window_kb))
}

data_10kb <- load_resolution("10kb", 10)
data_1kb  <- load_resolution("1kb", 1)

# ---- fixed y-axis for Panels A/B, shared across both phases for this strain --
all_depth <- c(data_1kb$log2_depth, data_10kb$log2_depth)
y_min <- floor(min(all_depth, na.rm = TRUE))
y_max <- ceiling(max(all_depth, na.rm = TRUE)) + 1
Y_LIMITS <- c(y_min, y_max)
Y_BREAKS <- pretty(Y_LIMITS, n = 6)
LABEL_Y  <- y_max - 0.2

# ---- Ori:Ter ratio (Panels A/B only) ----------------------------------------
ORI_ZONE <- 0.1

compute_ratio <- function(df) {
  ori_val <- median(df$log2_depth[abs(df$norm_position) <= ORI_ZONE])
  ter_val <- median(df$log2_depth[abs(df$norm_position) >= (1 - ORI_ZONE)])
  2^(ori_val - ter_val)
}

format_ratio <- function(label, ratio) {
  sprintf("%s  Ori:Ter = %.2f", label, ratio)
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

# ---- shared plot elements --------------------------------------------------
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
       y = expression(log[2]*"(read-start density, reads per kb)"))
)

ratio_annotation <- function(label) {
  annotate("label", x = Inf, y = Inf, label = label, hjust = 1.02, vjust = 1.3,
           size = 3.1, family = FONT, fill = alpha("white", 0.75))
}

panel_caption <- function(phase_label) {
  paste0(
    "Strain ", display_name, ", ", phase_label, " phase. oriC = 0, ter = ±1. ",
    "Ori:Ter is the median read-start density ratio between windows near oriC and near ter."
  )
}

# ---- Panel A/B builder (unchanged from plot_replication_profile.R) ---------
panel_AB <- function(trt, phs, title) {
  sub10 <- data_10kb %>% filter(treatment == trt, phase == phs)
  sub1  <- data_1kb  %>% filter(treatment == trt, phase == phs)

  ratio_a <- compute_ratio(sub10 %>% filter(replicate == "A"))
  ratio_b <- compute_ratio(sub10 %>% filter(replicate == "B"))
  ratio_label <- paste(format_ratio("A:", ratio_a), format_ratio("B:", ratio_b), sep = "\n")

  ggplot() +
    geom_point(data = sub1, aes(x = norm_position, y = log2_depth, color = replicate),
               alpha = 0.12, size = 0.5, show.legend = FALSE) +
    geom_point(data = sub10, aes(x = norm_position, y = log2_depth, color = replicate),
               alpha = 0.8, size = 1.3) +
    ori_ter_layers + axis_layers +
    scale_color_manual(values = palette_replicate, name = "Biological replicate") +
    ratio_annotation(ratio_label) +
    labs(title = title) +
    base_theme
}

# ---- fold-change data (Panel C) --------------------------------------------
build_fc_data <- function(phs) {
  sub <- data_1kb %>% filter(phase == phs)
  pick <- function(trt, rep, out_col) {
    sub %>% filter(treatment == trt, replicate == rep) %>%
      select(window_mid, norm_position, raw_count) %>%
      rename(!!out_col := raw_count)
  }
  drug_A <- pick("drug", "A", "drug_A")
  drug_B <- pick("drug", "B", "drug_B")
  nd_A   <- pick("ND",   "A", "nd_A")
  nd_B   <- pick("ND",   "B", "nd_B")

  drug_A %>%
    inner_join(drug_B, by = c("window_mid", "norm_position")) %>%
    inner_join(nd_A,   by = c("window_mid", "norm_position")) %>%
    inner_join(nd_B,   by = c("window_mid", "norm_position")) %>%
    mutate(
      density_drug = (drug_A + drug_B) / 2,
      density_nd   = (nd_A + nd_B) / 2,
      log2FC = log2(density_drug / density_nd)
    )
}

fc_data <- list(stat = build_fc_data("stat"), exp = build_fc_data("exp"))

all_fc <- unlist(lapply(fc_data, function(d) d$log2FC))
fc_max_abs <- ceiling(max(abs(all_fc[is.finite(all_fc)]), na.rm = TRUE) * 10) / 10 + 0.5
X_LIMITS <- c(-fc_max_abs, fc_max_abs)

fc_pooled_summary_records <- list()
fc_pooled_window_records  <- list()

record_fc_pooled_summary <- function(phs, n, n_excluded, tt, mw, n_up, n_down, n_ns) {
  fc_pooled_summary_records[[length(fc_pooled_summary_records) + 1]] <<- data.frame(
    strain = strain, phase = phs, n_windows = n, n_excluded = n_excluded,
    t_test_p = tt, mwu_p = mw, n_up = n_up, n_down = n_down, n_ns = n_ns
  )
}

record_fc_pooled_windows <- function(phs, fc) {
  fc_pooled_window_records[[length(fc_pooled_window_records) + 1]] <<- fc %>%
    mutate(strain = strain, phase = phs) %>%
    select(strain, phase, window_mid, norm_position, density_drug, density_nd, log2FC, category)
}

# ---- Panel C builder: log2FC per window vs genomic position ----------------
panel_fc_pooled <- function(phs, phase_label, title) {
  fc <- fc_data[[phs]]
  n_total <- length(unique((data_1kb %>% filter(phase == phs))$window_mid))
  n <- nrow(fc)
  n_excluded <- n_total - n

  fc <- fc %>% mutate(category = case_when(
    log2FC >  1 ~ "Enriched in Meropenem",
    log2FC < -1 ~ "Depleted in Meropenem",
    TRUE        ~ "No substantial change"
  ))
  n_up   <- sum(fc$category == "Enriched in Meropenem")
  n_down <- sum(fc$category == "Depleted in Meropenem")
  n_ns   <- sum(fc$category == "No substantial change")

  tt <- t.test(fc$density_drug, fc$density_nd)$p.value
  mw <- suppressWarnings(wilcox.test(fc$density_drug, fc$density_nd)$p.value)

  record_fc_pooled_summary(phs, n, n_excluded, tt, mw, n_up, n_down, n_ns)
  record_fc_pooled_windows(phs, fc)

  ann_label <- sprintf(
    "t-test p = %s\nMann-Whitney U p = %s\nn = %d windows (up: %d, down: %d, ns: %d)",
    format_p(tt), format_p(mw), n, n_up, n_down, n_ns
  )

  caption_c <- paste0(
    "Strain ", display_name, ", ", phase_label, " phase. Panel C: log2 fold-change ",
    "(Meropenem/untreated) per 1kb genomic window (mean of replicates A and B), plotted ",
    "against normalized genomic position (y; oriC = 0, ter = ± 1). Vertical lines mark ",
    "2-fold-change thresholds (log2FC = -1, 0, +1). ", n, " windows used (", n_excluded,
    " excluded for zero coverage in at least one condition): ", n_up, " above +1, ", n_down,
    " below -1, ", n_ns, " within ± 1. Overall t-test p = ", format_p(tt),
    "; overall Mann-Whitney U p = ", format_p(mw), " (both compare per-window mean read ",
    "density, Meropenem vs untreated, pooled across all windows genome-wide as the sample)."
  )

  p <- ggplot(fc, aes(x = log2FC, y = norm_position, color = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = -1, linetype = "dashed", color = unname(palette_fc["Depleted in Meropenem"]), linewidth = 0.5) +
    geom_vline(xintercept = 1,  linetype = "dashed", color = unname(palette_fc["Enriched in Meropenem"]), linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.1) +
    annotate("text", x = X_LIMITS[1] + 0.1, y = 0.05, label = "oriC", hjust = 0,
             size = 3.2, color = "grey35", family = FONT) +
    annotate("text", x = X_LIMITS[1] + 0.1, y = c(-0.95, 0.95), label = c("ter", "ter"), hjust = 0,
             size = 3.2, color = "grey55", family = FONT) +
    scale_color_manual(values = palette_fc, name = "log2FC category") +
    scale_x_continuous(limits = X_LIMITS, breaks = pretty(X_LIMITS, n = 6)) +
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = expression(log[2]*"(Meropenem / untreated)"),
         y = "Normalized position relative to oriC (oriC = 0, ter = ±1)",
         title = title) +
    annotate("label", x = Inf, y = Inf, label = ann_label, hjust = 1.02, vjust = 1.3,
             size = 3.0, family = FONT, fill = alpha("white", 0.75)) +
    base_theme

  list(plot = p, caption = caption_c)
}

# ---- save one standalone panel ----------------------------------------------
save_individual <- function(p, out_name, caption_text) {
  p_standalone <- p + labs(caption = paste(strwrap(caption_text, width = 95), collapse = "\n"))
  ggsave(file.path(individual_dir, paste0(out_name, ".png")), p_standalone,
         width = 6.2, height = 5.8, dpi = 300, bg = "white")
  ggsave(file.path(individual_dir, paste0(out_name, ".pdf")), p_standalone,
         width = 6.2, height = 5.8, device = cairo_pdf)
}

# ---- build one composite figure per growth phase ---------------------------
build_phase_figure <- function(phs, phase_label, out_name) {
  pA <- panel_AB("drug", phs, "Meropenem-treated: replicate A vs B")
  pB <- panel_AB("ND",   phs, "Untreated (no drug): replicate A vs B")
  fc_result <- panel_fc_pooled(phs, phase_label, "Meropenem vs untreated: log2 fold-change")
  pC <- fc_result$plot

  caption_ab <- panel_caption(phase_label)
  cap_wrapped <- paste(strwrap(fc_result$caption, width = 175), collapse = "\n")

  save_individual(pA, paste0(out_name, "_panelA_meropenem_repA_vs_repB"), caption_ab)
  save_individual(pB, paste0(out_name, "_panelB_untreated_repA_vs_repB"), caption_ab)
  save_individual(pC, paste0(out_name, "_panelC_fc_pooled"), fc_result$caption)

  combined <- (pA | pB | pC) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  final <- combined +
    plot_annotation(
      title = paste0("Replication Profile (Fold-Change, Pooled Test): Acinetobacter baumannii Strain ", display_name, ", ", phase_label, " Phase"),
      caption = cap_wrapped,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, family = FONT, hjust = 0),
        plot.caption = element_text(size = 10, family = FONT, hjust = 0, color = "grey20"),
        plot.caption.position = "plot"
      )
    )

  ggsave(file.path(out_dir, paste0(out_name, ".png")), final, width = 15, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(out_name, ".pdf")), final, width = 15, height = 5.5, device = cairo_pdf)
  message("Saved: ", out_name, ".png / .pdf (composite) + 3 individual panels")
}

build_phase_figure("stat", "Stationary", "stationary_phase")
build_phase_figure("exp",  "Exponential", "exponential_phase")

write_csv(do.call(rbind, fc_pooled_summary_records), file.path(out_dir, "..", "fc_pooled_summary.csv"))
write_csv(do.call(rbind, fc_pooled_window_records), file.path(out_dir, "..", "fc_pooled_windows.csv"))

message("Composite + individual figures written to ", out_dir, " (individual/ subfolder)")
message("Fold-change summary written to ", file.path(out_dir, "..", "fc_pooled_summary.csv"))
message("Per-window fold-change values written to ", file.path(out_dir, "..", "fc_pooled_windows.csv"))
