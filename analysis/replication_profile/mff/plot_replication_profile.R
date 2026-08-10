#!/usr/bin/env Rscript
# Publication-ready marker-frequency / replication-profile plots, following
# Skovgaard et al. 2011 (Genome Research) methodology: x = genomic position
# normalized per replichore (oriC = 0, ter = +-1), y = log2(marker frequency)
# relative to the genome-wide median window count. Exponential-phase cultures
# are expected to show a "tent" peak at oriC tapering to ter; stationary-phase
# cultures (non-replicating) should be flat.
#
# Usage:
#   Rscript plot_replication_profile.R [coverage_dir] [metadata_csv] [strain] [out_dir]
# Defaults are set up for the mff strain as run from the project root.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

argv <- commandArgs(trailingOnly = TRUE)
coverage_dir <- if (length(argv) >= 1) argv[1] else "mapping/coverage/mff"
metadata_csv <- if (length(argv) >= 2) argv[2] else "mapping/coverage/sample_metadata.csv"
strain        <- if (length(argv) >= 3) argv[3] else "mff"
out_dir       <- if (length(argv) >= 4) argv[4] else file.path("analysis/replication_profile", strain, "figures")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    left_join(metadata, by = "pid_sample")
}

data_10kb <- load_resolution("10kb")
data_1kb  <- load_resolution("1kb")

# average the two biological replicates (A & B) per treatment/phase, joining on window_mid
average_replicates <- function(df, trt, phs) {
  sub <- df %>% filter(treatment == trt, phase == phs)
  sub %>%
    group_by(window_mid, norm_position) %>%
    summarise(log2_ratio = mean(log2_ratio), .groups = "drop") %>%
    mutate(treatment = trt, phase = phs)
}

palette_replicate <- c(A = "#1b9e77", B = "#d95f02")
palette_treatment <- c(drug = "#7570b3", ND = "#e7298a")

base_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ori_ter_layers <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30"),
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60"),
  annotate("text", x = 0, y = Inf, label = "oriC", vjust = 1.5, hjust = -0.1, size = 3.5, color = "grey30"),
  annotate("text", x = c(-1, 1), y = Inf, label = c("ter", "ter"), vjust = 1.5, hjust = c(-0.1, 1.1), size = 3.5, color = "grey50")
)

axis_layers <- list(
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)),
  labs(x = "Normalized position relative to oriC (oriC = 0, ter = ±1)",
       y = expression(log[2]*"(marker frequency)"))
)

save_plot <- function(p, name) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = 7.5, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = 7.5, height = 5)
  message("Saved: ", name, ".png / .pdf")
}

plot_AB <- function(trt, phs, title, name) {
  sub10 <- data_10kb %>% filter(treatment == trt, phase == phs)
  sub1  <- data_1kb  %>% filter(treatment == trt, phase == phs)

  p <- ggplot() +
    geom_point(data = sub1, aes(x = norm_position, y = log2_ratio, color = replicate),
               alpha = 0.15, size = 0.6, show.legend = FALSE) +
    geom_line(data = sub10, aes(x = norm_position, y = log2_ratio, color = replicate), linewidth = 0.9) +
    geom_point(data = sub10, aes(x = norm_position, y = log2_ratio, color = replicate), size = 1.2) +
    ori_ter_layers + axis_layers +
    scale_color_manual(values = palette_replicate, name = "Replicate") +
    labs(title = title, subtitle = paste0(strain, " — ", trt, ", ", phs)) +
    base_theme
  save_plot(p, name)
}

plot_comparison <- function(phs, title, name) {
  sub10_drug <- average_replicates(data_10kb, "drug", phs)
  sub10_nd   <- average_replicates(data_10kb, "ND", phs)
  sub10 <- rbind(sub10_drug, sub10_nd)
  sub1_drug <- average_replicates(data_1kb, "drug", phs)
  sub1_nd   <- average_replicates(data_1kb, "ND", phs)
  sub1 <- rbind(sub1_drug, sub1_nd)

  p <- ggplot() +
    geom_point(data = sub1, aes(x = norm_position, y = log2_ratio, color = treatment),
               alpha = 0.15, size = 0.6, show.legend = FALSE) +
    geom_line(data = sub10, aes(x = norm_position, y = log2_ratio, color = treatment), linewidth = 0.9) +
    geom_point(data = sub10, aes(x = norm_position, y = log2_ratio, color = treatment), size = 1.2) +
    ori_ter_layers + axis_layers +
    scale_color_manual(values = palette_treatment, name = "Treatment", labels = c(drug = "Drug", ND = "Untreated (ND)")) +
    labs(title = title, subtitle = paste0(strain, " — ", phs, ", replicate-averaged (A & B)")) +
    base_theme
  save_plot(p, name)
}

plot_AB("drug", "stat", "Replication profile: drug-treated, stationary phase (A vs B)", "01_drug_stationary_A_vs_B")
plot_AB("ND",   "stat", "Replication profile: untreated, stationary phase (A vs B)",    "02_ND_stationary_A_vs_B")
plot_comparison("stat", "Replication profile: drug vs untreated, stationary phase",     "03_drug_vs_ND_stationary")
plot_AB("drug", "exp",  "Replication profile: drug-treated, exponential phase (A vs B)", "04_drug_exponential_A_vs_B")
plot_AB("ND",   "exp",  "Replication profile: untreated, exponential phase (A vs B)",    "05_ND_exponential_A_vs_B")
plot_comparison("exp",  "Replication profile: drug vs untreated, exponential phase",     "06_drug_vs_ND_exponential")

message("All 6 plots written to ", out_dir)
