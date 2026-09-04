#!/usr/bin/env Rscript
# =============================================================================
# plot_final_vs_tifseq_comparison.R
# =============================================================================
# Compares this study's FINAL merged UTR annotation (final_utr.tsv) against an
# INDEPENDENT orthogonal benchmark: the Pelechano, Wei & Steinmetz (2013)
# TIF-seq major transcript isoforms (Pelechano_UTR.csv, derived by
# derive_pelechano_utr.py). Restricted to genes where both sources report a
# non-zero UTR length.
#
# This mirrors plot_suppl_fig3.R (the DRS-vs-Nagalakshmi comparison) but uses
# SIGNED differences in both directions, because TIF-seq is an independent
# benchmark and is NOT part of the merge that built the final annotation — the
# final call can be either longer OR shorter than the TIF-seq estimate.
#
# Panels:
#   A  Histogram of (Final - TIF-seq) differences, 5' UTR
#   B  Histogram of (Final - TIF-seq) differences, 3' UTR
#   C  Grouped bar chart: gene counts in each category (5' and 3')
#   D  Scatter: TIF-seq vs Final UTR length (log scale, both ends) + Spearman r_s
#
# Categories (threshold = 20 nt):
#   "TIF-seq longer (>20 nt)"  final call shorter than TIF-seq
#   "Within +/-20 nt"          approximately equal (noise level)
#   "Final longer (>20 nt)"    final call longer than TIF-seq
#
# Colour scheme (same palette as plot_suppl_fig3.R):
#   TIF-seq longer:  #AA3377  (purple)
#   Within +/-20:    #AAAAAA  (grey)
#   Final longer:    #228833  (green)
#
# Output: pdf, png, svg -> --outdir
#
# Usage:
#   Rscript plot_final_vs_tifseq_comparison.R \
#     --final_utr /path/final_utr.tsv \
#     --ref_utr   /path/Pelechano_UTR.csv \
#     --outdir    /path/Figures \
#     --outname   final_vs_tifseq_comparison \
#     --threshold 20
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# -- Argument parsing ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("No value supplied for", flag))
  args[i + 1]
}

fin_file <- get_arg("--final_utr")
ref_file <- get_arg("--ref_utr")
out_dir  <- get_arg("--outdir")
if (is.null(fin_file) || is.null(ref_file) || is.null(out_dir))
  stop("Required: --final_utr, --ref_utr, --outdir")
out_name  <- get_arg("--outname",  "final_vs_tifseq_comparison")
threshold <- as.integer(get_arg("--threshold", "20"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", out_dir, "\n")
cat("Threshold (nt):  ", threshold, "\n\n")

# -- Load data ----------------------------------------------------------------
fin <- read.table(fin_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(fin) <- c("Gene", "fin_five", "fin_three")

ref <- read.csv(ref_file, header = TRUE, stringsAsFactors = FALSE)
colnames(ref) <- c("Gene", "ref_five", "ref_three")
ref$ref_five  <- suppressWarnings(as.numeric(ref$ref_five))
ref$ref_three <- suppressWarnings(as.numeric(ref$ref_three))

# -- Merge and build per-end tables -------------------------------------------
df <- inner_join(fin, ref, by = "Gene")

cat_labels <- c(
  paste0("TIF-seq longer (>", threshold, " nt)"),
  paste0("Within \u00b1",     threshold, " nt"),
  paste0("Final longer (>",   threshold, " nt)")
)

classify <- function(diff) {
  case_when(
    diff <= -threshold ~ cat_labels[1],   # final shorter -> TIF-seq longer
    diff >=  threshold ~ cat_labels[3],   # final longer
    TRUE               ~ cat_labels[2]
  )
}

df5 <- df %>%
  filter(!is.na(ref_five), ref_five > 0, fin_five > 0) %>%
  mutate(diff = fin_five - ref_five, category = classify(diff), UTR = "5' UTR")

df3 <- df %>%
  filter(!is.na(ref_three), ref_three > 0, fin_three > 0) %>%
  mutate(diff = fin_three - ref_three, category = classify(diff), UTR = "3' UTR")

cat("5' UTR - genes with both annotations non-zero:", nrow(df5), "\n")
cat("3' UTR - genes with both annotations non-zero:", nrow(df3), "\n\n")

# Spearman correlation (rank-based, robust to outliers) per end
rs5 <- suppressWarnings(cor(df5$ref_five,  df5$fin_five,  method = "spearman"))
rs3 <- suppressWarnings(cor(df3$ref_three, df3$fin_three, method = "spearman"))
cat(sprintf("Spearman r_s  5' = %.3f   3' = %.3f\n\n", rs5, rs3))
cat("Category breakdown:\n")
cat("5' UTR:\n"); print(table(factor(df5$category, levels = cat_labels)))
cat("3' UTR:\n"); print(table(factor(df3$category, levels = cat_labels)))

# -- Colours and factor order -------------------------------------------------
cat_colours <- c("#AA3377", "#AAAAAA", "#228833")
names(cat_colours) <- cat_labels
df5$category <- factor(df5$category, levels = cat_labels)
df3$category <- factor(df3$category, levels = cat_labels)

# -- Shared theme -------------------------------------------------------------
base_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey60"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 10),
    plot.title       = element_text(size = 11, face = "bold"),
    plot.tag         = element_text(size = 13, face = "bold")
  )

# -- Helper: histogram panel --------------------------------------------------
make_hist <- function(dat, utr_label, xlim_cap = 600) {
  ggplot(dat, aes(x = diff, fill = category)) +
    geom_histogram(binwidth = 10, colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = 0,          colour = "black",  lwd = 0.8) +
    geom_vline(xintercept =  threshold, colour = "grey40", lwd = 0.6, linetype = "dashed") +
    geom_vline(xintercept = -threshold, colour = "grey40", lwd = 0.6, linetype = "dashed") +
    scale_fill_manual(values = cat_colours) +
    coord_cartesian(xlim = c(-xlim_cap, xlim_cap)) +
    labs(
      title = utr_label,
      x     = "Final UTR - TIF-seq UTR (nt)",
      y     = "Number of genes",
      tag   = if (utr_label == "5' UTR") "A" else "B"
    ) +
    annotate("text",
             x = c(-xlim_cap * 0.72, xlim_cap * 0.72),
             y = Inf, vjust = 1.3, size = 3.0,
             colour = c("#AA3377", "#228833"), fontface = "bold",
             label = c(
               paste0("TIF-seq longer\nn = ", sum(dat$category == cat_labels[1])),
               paste0("Final longer\nn = ",   sum(dat$category == cat_labels[3]))
             )) +
    annotate("text",
             x = 0, y = Inf, vjust = 1.3, size = 3.0,
             colour = "#666666", fontface = "bold",
             label = paste0("Within \u00b1", threshold, " nt\nn = ",
                            sum(dat$category == cat_labels[2]))) +
    annotate("text",
             x = xlim_cap, y = 0, hjust = 1, vjust = 0, size = 2.8,
             colour = "grey50",
             label = paste0("x-axis capped at \u00b1", xlim_cap, " nt")) +
    base_theme
}

pA <- make_hist(df5, "5' UTR")
pB <- make_hist(df3, "3' UTR")

# -- Panel C: grouped bar chart -----------------------------------------------
df_bar <- bind_rows(df5, df3) %>%
  count(UTR, category) %>%
  mutate(UTR = factor(UTR, levels = c("5' UTR", "3' UTR")))

pC <- ggplot(df_bar, aes(x = UTR, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cat_colours, name = NULL,
                    guide = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Gene counts by category", x = NULL, y = "Number of genes", tag = "C") +
  base_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

# -- Panel D: scatter TIF-seq vs Final (log scale, faceted) + Spearman --------
df_scatter <- bind_rows(df5, df3) %>%
  mutate(
    ref_val = if_else(UTR == "5' UTR", ref_five, ref_three),
    fin_val = if_else(UTR == "5' UTR", fin_five, fin_three),
    UTR = factor(UTR, levels = c("5' UTR", "3' UTR"))
  ) %>%
  filter(ref_val > 0, fin_val > 0)

rs_lab <- data.frame(
  UTR = factor(c("5' UTR", "3' UTR"), levels = c("5' UTR", "3' UTR")),
  label = c(sprintf("r[s] == %.2f", rs5), sprintf("r[s] == %.2f", rs3))
)

pD <- ggplot(df_scatter, aes(x = ref_val, y = fin_val, colour = category)) +
  geom_abline(slope = 1, intercept = 0, colour = "black",
              linetype = "dashed", lwd = 0.7) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = cat_colours) +
  facet_wrap(~UTR) +
  geom_text(data = rs_lab, aes(x = 1.5, y = Inf, label = label),
            parse = TRUE, inherit.aes = FALSE, hjust = 0, vjust = 1.6,
            size = 3.4, colour = "grey20") +
  labs(title = "TIF-seq vs Final UTR lengths",
       x = "TIF-seq UTR length (nt, log scale)",
       y = "Final UTR length (nt, log scale)", tag = "D") +
  base_theme +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

# -- Assemble -----------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    caption = paste0(
      "Comparison of the final merged annotation against the independent ",
      "Pelechano et al. (2013) TIF-seq benchmark, restricted to genes where both ",
      "report a non-zero UTR length (5' UTR: n = ", nrow(df5),
      "; 3' UTR: n = ", nrow(df3), "). ",
      "Differences are signed (Final - TIF-seq); a ", threshold,
      " nt threshold excludes differences plausibly attributable to noise. ",
      "Histogram x-axes capped at \u00b1600 nt. Dashed line in panel D = y = x. ",
      "Spearman r_s 5' = ", sprintf("%.2f", rs5), ", 3' = ", sprintf("%.2f", rs3), ". ",
      "Purple: TIF-seq longer; green: final annotation longer."
    ),
    theme = theme(plot.caption = element_text(size = 7.5, colour = "grey40",
                                              lineheight = 1.3))
  )

# -- Save ---------------------------------------------------------------------
w <- 11; h <- 9
out_base <- file.path(out_dir, out_name)
ggsave(paste0(out_base, ".pdf"), fig, width = w, height = h, device = "pdf")
cat("Saved:", paste0(out_base, ".pdf"), "\n")
ggsave(paste0(out_base, ".png"), fig, width = w, height = h, dpi = 300, device = "png")
cat("Saved:", paste0(out_base, ".png"), "\n")
ggsave(paste0(out_base, ".svg"), fig, width = w, height = h, device = "svg")
cat("Saved:", paste0(out_base, ".svg"), "\n")
cat("\nDone.\n")
