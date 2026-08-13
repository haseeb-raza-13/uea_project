#!/usr/bin/env Rscript
# Fold-change variant of Panel C only (Panels A and B are identical to the
# log2 version, plot_replication_profile.R). Panel C here is a classic
# differential-expression-style volcano plot: per 1kb genomic window,
# x = log2(Meropenem density / untreated density), y = -log10(p-value) from a
# per-window Welch t-test on the 2-vs-2 replicate values.
#
# Statistical design ("Design 2: classic per-window volcano"). Caveat, made
# explicit in the plot caption: with only 2 biological replicates per
# treatment, the per-window t-test has just 2 degrees of freedom, and a
# per-window Mann-Whitney U test is severely under-powered at n=2 vs n=2 (the
# exact test's minimum attainable two-sided p-value is 1/3 when the four
# values are untied; read counts are discrete and frequently tied in practice,
# which makes R fall back to a normal approximation that can occasionally read
# below that nominal floor). Both p-values are nonetheless computed and
# written to the results CSV for every window (t-test raw and BH-adjusted,
# Mann-Whitney U raw); only the t-test p-value drives the plotted y-axis,
# since Mann-Whitney U cannot produce a meaningful continuous axis at this
# replicate count. See plot_replication_profile_fc_pooled.R for the
# statistically well-powered pooled/overall-test variant (Design 1).
#
# Usage:
#   Rscript plot_replication_profile_fc_volcano.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]

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
  stop("Usage: plot_replication_profile_fc_volcano.R <strain> <display_name> [coverage_dir] [metadata_csv] [out_dir]")
}
strain        <- argv[1]
display_name  <- argv[2]
coverage_dir  <- if (length(argv) >= 3) argv[3] else file.path("mapping/coverage", strain)
metadata_csv  <- if (length(argv) >= 4) argv[4] else "mapping/coverage/sample_metadata.csv"
out_dir       <- if (length(argv) >= 5) argv[5] else file.path("analysis/replication_profile", strain, "figures_fc_volcano")
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
palette_volcano <- c("Up (Meropenem-enriched)"   = dark2[1],
                      "Down (Meropenem-depleted)" = dark2[2],
                      "Not significant"           = "grey70")

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

# ---- fold-change + per-window significance data (Panel C) ------------------
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

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

  joined <- drug_A %>%
    inner_join(drug_B, by = c("window_mid", "norm_position")) %>%
    inner_join(nd_A,   by = c("window_mid", "norm_position")) %>%
    inner_join(nd_B,   by = c("window_mid", "norm_position")) %>%
    mutate(
      density_drug = (drug_A + drug_B) / 2,
      density_nd   = (nd_A + nd_B) / 2,
      log2FC = log2(density_drug / density_nd)
    )

  pvals <- mapply(function(a1, a2, b1, b2) {
    x <- c(a1, a2); y <- c(b1, b2)
    tt <- tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
    mw <- tryCatch(suppressWarnings(wilcox.test(x, y)$p.value), error = function(e) NA_real_)
    c(tt = tt, mw = mw)
  }, joined$drug_A, joined$drug_B, joined$nd_A, joined$nd_B)

  joined$t_test_p        <- pvals["tt", ]
  joined$mwu_p            <- pvals["mw", ]
  joined$t_test_p_adj_BH  <- p.adjust(joined$t_test_p, method = "BH")
  joined$neg_log10_p      <- -log10(joined$t_test_p)
  joined$category <- case_when(
    joined$log2FC >  1 & joined$t_test_p < 0.05 ~ "Up (Meropenem-enriched)",
    joined$log2FC < -1 & joined$t_test_p < 0.05 ~ "Down (Meropenem-depleted)",
    TRUE ~ "Not significant"
  )
  joined
}

fc_data <- list(stat = build_fc_data("stat"), exp = build_fc_data("exp"))

all_fc <- unlist(lapply(fc_data, function(d) d$log2FC))
fc_max_abs <- ceiling(max(abs(all_fc[is.finite(all_fc)]), na.rm = TRUE) * 10) / 10 + 0.5
X_LIMITS <- c(-fc_max_abs, fc_max_abs)

all_neg_log10_p <- unlist(lapply(fc_data, function(d) d$neg_log10_p))
Y_MAX_V <- ceiling(max(all_neg_log10_p[is.finite(all_neg_log10_p)], na.rm = TRUE) * 10) / 10 + 0.3
Y_LIMITS_V <- c(0, Y_MAX_V)

fc_volcano_window_records <- list()

record_fc_volcano_windows <- function(phs, fc) {
  fc_volcano_window_records[[length(fc_volcano_window_records) + 1]] <<- fc %>%
    mutate(strain = strain, phase = phs) %>%
    select(strain, phase, window_mid, norm_position, density_drug, density_nd, log2FC,
           t_test_p, t_test_p_adj_BH, mwu_p, category)
}

# ---- Panel C builder: volcano plot ------------------------------------------
panel_fc_volcano <- function(phs, phase_label, title) {
  fc <- fc_data[[phs]]
  n_total <- length(unique((data_1kb %>% filter(phase == phs))$window_mid))
  n <- nrow(fc)
  n_excluded <- n_total - n

  n_up      <- sum(fc$category == "Up (Meropenem-enriched)")
  n_down    <- sum(fc$category == "Down (Meropenem-depleted)")
  n_mwu_sig <- sum(fc$mwu_p < 0.05, na.rm = TRUE)

  record_fc_volcano_windows(phs, fc)

  caption_c <- paste0(
    "Strain ", display_name, ", ", phase_label, " phase. Panel C: volcano plot, ",
    "log2(Meropenem/untreated) per 1kb window (x) vs -log10(t-test p-value) (y); n = ", n,
    " windows (", n_excluded, " excluded for zero coverage in at least one condition); ",
    n_up, " significantly up, ", n_down, " significantly down (|log2FC| > 1 and t-test p < 0.05). ",
    "Vertical lines mark log2FC = -1/+1; horizontal line marks p = 0.05. Caution: p-values are ",
    "per-window from only 2 replicates per condition (t-test df = 2); BH-adjusted p-values are in ",
    "the results CSV. Mann-Whitney U p-values (also in the CSV) are of limited value at this ",
    "replicate count: the exact test's minimum attainable two-sided p-value is 1/3 with no ties among ",
    "the four values, but read counts are discrete and frequently tied, which makes R fall back to a ",
    "normal approximation that can occasionally read below that nominal floor. Either way it is ",
    "under-powered here and is not used to drive the y-axis (", n_mwu_sig, " of ", n,
    " windows show raw Mann-Whitney U p < 0.05, for reference)."
  )

  p <- ggplot(fc, aes(x = log2FC, y = neg_log10_p, color = category)) +
    geom_vline(xintercept = -1, linetype = "dashed", color = unname(palette_volcano["Down (Meropenem-depleted)"]), linewidth = 0.5) +
    geom_vline(xintercept = 1,  linetype = "dashed", color = unname(palette_volcano["Up (Meropenem-enriched)"]), linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.1) +
    scale_color_manual(values = palette_volcano, name = "Significance (t-test)") +
    scale_x_continuous(limits = X_LIMITS, breaks = pretty(X_LIMITS, n = 6)) +
    scale_y_continuous(limits = Y_LIMITS_V) +
    labs(x = expression(log[2]*"(Meropenem / untreated)"),
         y = expression(-log[10]*"(t-test p-value)"),
         title = title) +
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
  fc_result <- panel_fc_volcano(phs, phase_label, "Meropenem vs untreated: volcano plot")
  pC <- fc_result$plot

  caption_ab <- panel_caption(phase_label)
  cap_wrapped <- paste(strwrap(fc_result$caption, width = 175), collapse = "\n")

  save_individual(pA, paste0(out_name, "_panelA_meropenem_repA_vs_repB"), caption_ab)
  save_individual(pB, paste0(out_name, "_panelB_untreated_repA_vs_repB"), caption_ab)
  save_individual(pC, paste0(out_name, "_panelC_fc_volcano"), fc_result$caption)

  combined <- (pA | pB | pC) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  final <- combined +
    plot_annotation(
      title = paste0("Replication Profile (Fold-Change, Volcano Plot): Acinetobacter baumannii Strain ", display_name, ", ", phase_label, " Phase"),
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

write_csv(do.call(rbind, fc_volcano_window_records), file.path(out_dir, "..", "fc_volcano_windows.csv"))

message("Composite + individual figures written to ", out_dir, " (individual/ subfolder)")
message("Per-window volcano results written to ", file.path(out_dir, "..", "fc_volcano_windows.csv"))
