#!/usr/bin/env Rscript
# =============================================================================
# plot_drs_vs_nagalakshmi_comparison.R
# =============================================================================
# Produces a 4-panel supplementary figure comparing the pre-merge DRS-derived
# UTR annotation (DRS_UTR_corrected.tsv) against the Nagalakshmi et al. (2008)
# annotation, restricted to genes where both sources report a non-zero length.
#
# Panels:
#   A  Histogram of (DRS − Naga) differences, 5′ UTR
#   B  Histogram of (DRS − Naga) differences, 3′ UTR
#   C  Grouped bar chart: gene counts in each category (5′ and 3′)
#   D  Scatter: Nagalakshmi vs DRS UTR length (log scale, both UTR ends)
#
# Categories (threshold = 20 nt):
#   "DRS shorter (>20 nt)"   DRS call more conservative than Naga
#   "Within ±20 nt"          Approximately equal (noise level)
#   "DRS longer (>20 nt)"    DRS extends beyond Naga
#
# Colour scheme:
#   DRS shorter:  #AA3377  (purple — Nagalakshmi retained by merge)
#   Within ±20:   #AAAAAA  (grey)
#   DRS longer:   #228833  (green — DRS call retained by merge)
#
# Output: pdf, png, svg — all written to --outdir
#
# Usage:
#   Rscript plot_drs_vs_nagalakshmi_comparison.R \
#     --drs_utr   /path/DRS_UTR_corrected.tsv \
#     --ref_utr   /path/reference_UTR.csv \
#     --outdir    /path/Figures \
#     --outname   drs_vs_nag_comparison \
#     --threshold 20
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ── Argument parsing ----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) stop(paste("No value supplied for", flag))
  args[i + 1]
}

drs_file  <- get_arg("--drs_utr")
nag_file  <- get_arg("--ref_utr")
out_dir   <- get_arg("--outdir")

if (is.null(drs_file) || is.null(nag_file) || is.null(out_dir))
  stop("Required: --drs_utr, --ref_utr, --outdir")
out_name  <- get_arg("--outname",  "drs_vs_nag_comparison")
threshold <- as.integer(get_arg("--threshold", "20"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory:", out_dir, "\n")
cat("Threshold (nt):  ", threshold, "\n\n")

# ── Load data -----------------------------------------------------------------
drs <- read.table(drs_file, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE)
colnames(drs) <- c("Gene", "drs_five", "drs_three")

nag <- read.csv(nag_file, header = TRUE, stringsAsFactors = FALSE)
colnames(nag) <- c("Gene", "nag_five", "nag_three")
nag$nag_five  <- suppressWarnings(as.numeric(nag$nag_five))
nag$nag_three <- suppressWarnings(as.numeric(nag$nag_three))

# ── Merge and filter ----------------------------------------------------------
# Keep only genes with a real (non-NA, >0) annotation on both sides per end
df <- inner_join(drs, nag, by = "Gene")

df5 <- df %>%
  filter(!is.na(nag_five), nag_five > 0, drs_five > 0) %>%
  mutate(
    diff = drs_five - nag_five,
    category = case_when(
      diff <= -threshold ~ paste0("DRS shorter (>", threshold, " nt)"),
      diff >=  threshold ~ paste0("DRS longer (>",  threshold, " nt)"),
      TRUE               ~ paste0("Within \u00b1", threshold, " nt")
    ),
    UTR = "5' UTR"
  )

df3 <- df %>%
  filter(!is.na(nag_three), nag_three > 0, drs_three > 0) %>%
  mutate(
    diff = drs_three - nag_three,
    category = case_when(
      diff <= -threshold ~ paste0("DRS shorter (>", threshold, " nt)"),
      diff >=  threshold ~ paste0("DRS longer (>",  threshold, " nt)"),
      TRUE               ~ paste0("Within \u00b1", threshold, " nt")
    ),
    UTR = "3' UTR"
  )

cat("5' UTR — genes with both annotations non-zero:", nrow(df5), "\n")
cat("3' UTR — genes with both annotations non-zero:", nrow(df3), "\n\n")
cat("Category breakdown:\n")
cat("5' UTR:\n");  print(table(df5$category))
cat("3' UTR:\n");  print(table(df3$category))

# ── Colours and category levels -----------------------------------------------
cat_labels <- c(
  paste0("DRS shorter (>", threshold, " nt)"),
  paste0("Within \u00b1",  threshold, " nt"),
  paste0("DRS longer (>",  threshold, " nt)")
)
cat_colours <- c("#AA3377", "#AAAAAA", "#228833")
names(cat_colours) <- cat_labels

# Ensure consistent factor order throughout
df5$category <- factor(df5$category, levels = cat_labels)
df3$category <- factor(df3$category, levels = cat_labels)

# ── Shared theme --------------------------------------------------------------
base_theme <- theme_bw(base_size = 11) +
  theme(
    strip.background  = element_rect(fill = "grey92", colour = "grey60"),
    strip.text        = element_text(face = "bold"),
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    axis.title        = element_text(size = 10),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.tag          = element_text(size = 13, face = "bold")
  )

# ── Helper: build histogram panel for one UTR end ----------------------------
make_hist <- function(dat, utr_label, xlim_cap = 600) {

  # count annotations: compute outside ggplot for placement
  counts <- dat %>%
    count(category) %>%
    mutate(
      x_pos = c(-xlim_cap * 0.72, 0, xlim_cap * 0.72),
      label = paste0(category, "\nn = ", n)
    )

  ggplot(dat, aes(x = diff, fill = category)) +
    geom_histogram(binwidth = 10, colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = 0,          colour = "black",   lwd = 0.8) +
    geom_vline(xintercept =  threshold, colour = "grey40",  lwd = 0.6,
               linetype = "dashed") +
    geom_vline(xintercept = -threshold, colour = "grey40",  lwd = 0.6,
               linetype = "dashed") +
    scale_fill_manual(values = cat_colours) +
    coord_cartesian(xlim = c(-xlim_cap, xlim_cap)) +
    labs(
      title = utr_label,
      x     = "DRS UTR - Reference UTR (nt)",
      y     = "Number of genes",
      tag   = if (utr_label == "5' UTR") "A" else "B"
    ) +
    annotate("text",
             x = c(-xlim_cap * 0.72, xlim_cap * 0.72),
             y = Inf, vjust = 1.3, size = 3.0,
             colour = c("#AA3377", "#228833"),
             fontface = "bold",
             label = c(
               paste0("DRS shorter\nn = ",
                      sum(dat$category == cat_labels[1])),
               paste0("DRS longer\nn = ",
                      sum(dat$category == cat_labels[3]))
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

# ── Panel A: 5' UTR histogram -------------------------------------------------
pA <- make_hist(df5, "5' UTR")

# ── Panel B: 3' UTR histogram -------------------------------------------------
pB <- make_hist(df3, "3' UTR")

# ── Panel C: Grouped bar chart ------------------------------------------------
df_bar <- bind_rows(df5, df3) %>%
  count(UTR, category) %>%
  mutate(UTR = factor(UTR, levels = c("5' UTR", "3' UTR")))

pC <- ggplot(df_bar, aes(x = UTR, y = n, fill = category)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           colour = "white", linewidth = 0.3) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = cat_colours,
                    name   = NULL,
                    guide  = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Gene counts by category",
    x     = NULL,
    y     = "Number of genes",
    tag   = "C"
  ) +
  base_theme +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.45, "cm"))

# ── Panel D: Scatter Naga vs DRS (log scale, faceted) ------------------------
df_scatter <- bind_rows(df5, df3) %>%
  # rename so we have nag and drs columns
  mutate(
    nag = if_else(UTR == "5' UTR",
                  nag_five,
                  nag_three),
    drs_val = if_else(UTR == "5' UTR",
                      drs_five,
                      drs_three),
    UTR = factor(UTR, levels = c("5' UTR", "3' UTR"))
  ) %>%
  filter(nag > 0, drs_val > 0)

# For the y=x annotation: range across both facets
xy_range <- range(c(df_scatter$nag, df_scatter$drs_val), na.rm = TRUE)

pD <- ggplot(df_scatter,
             aes(x = nag, y = drs_val, colour = category)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "black", linetype = "dashed", lwd = 0.7) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = cat_colours) +
  facet_wrap(~UTR) +
  labs(
    title = "Reference vs DRS UTR lengths",
    x     = "Reference UTR length (nt, log scale)",
    y     = "DRS UTR length (nt, log scale)",
    tag   = "D"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold", size = 10)
  ) +
  annotate("text",
           x = max(df_scatter$nag, na.rm = TRUE),
           y = max(df_scatter$nag, na.rm = TRUE),
           hjust = 1.05, vjust = -0.3, size = 2.8,
           colour = "grey40", label = "y = x")

# ── Assemble with patchwork ---------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
  plot_annotation(
    caption = paste0(
      "Comparison restricted to genes where both DRS and Nagalakshmi report a non-zero UTR length ",
      "(5' UTR: n = ", nrow(df5), "; 3' UTR: n = ", nrow(df3), "). ",
      "Categories defined by a ", threshold, " nt threshold to exclude differences ",
      "plausibly attributable to annotation noise. ",
      "Histogram x-axes capped at \u00b1600 nt; extreme values fall outside the displayed range. ",
      "Dashed black line in panel D = y = x (equal UTR lengths). ",
      "Purple: Nagalakshmi annotation longer (retained by merge). ",
      "Green: DRS annotation longer (DRS call retained by merge)."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, colour = "grey40",
                                  lineheight = 1.3)
    )
  )

# ── Save ---------------------------------------------------------------------
w <- 11   # inches
h <- 9    # inches

out_base <- file.path(out_dir, out_name)

ggsave(paste0(out_base, ".pdf"), fig, width = w, height = h, device = "pdf")
cat("Saved:", paste0(out_base, ".pdf"), "\n")

ggsave(paste0(out_base, ".png"), fig, width = w, height = h,
       dpi = 300, device = "png")
cat("Saved:", paste0(out_base, ".png"), "\n")

ggsave(paste0(out_base, ".svg"), fig, width = w, height = h, device = "svg")
cat("Saved:", paste0(out_base, ".svg"), "\n")

cat("\nDone.\n")
