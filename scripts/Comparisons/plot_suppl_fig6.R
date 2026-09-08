#!/usr/bin/env Rscript
# =============================================================================
# plot_suppl_fig6.R
# =============================================================================
# Supplementary Figure S6. Over-extension and under-extension against
# TIF-seq, reported explicitly in both directions, for three annotations
# at three tolerances -- the report Reviewer 2 asked for in place of a
# single under-extension tail folded into "otherwise concordant or longer".
#
# Panels
#   A  Both tails at the 20 nt tolerance, for the released annotation,
#      the prior reference, and the raw DRS calls. The prior reference
#      is the control the manuscript's first submission never showed.
#   B  Both rates as a continuous function of the tolerance, so the
#      reader can see the choice of tolerance is not load-bearing.
#   C  Dispersion of the benchmark around its own point estimate: sets
#      the floor below which a tolerance measures TIF-seq isoform
#      heterogeneity rather than annotation error.
#   D  The envelope test: whether each released boundary falls inside
#      the range of ends TIF-seq actually observed for that gene.
#
# Inputs are the outputs of tifseq_rates.py. Nothing is recomputed here.
#
# Usage:
#   Rscript plot_suppl_fig6.R --indir Data/TIFseq --outdir Figures_out
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  args[i + 1]
}
IN  <- get_arg("--indir")
OUT <- get_arg("--outdir")
if (is.null(IN) || is.null(OUT)) stop("Required: --indir, --outdir")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COL_UNDER <- "#AA3377"; COL_CONC <- "#AAAAAA"; COL_OVER <- "#228833"
ANNOT_LAB <- c(final = "Final\nmerged", nagalakshmi = "Prior\nreference", drs = "DRS calls\nonly")
ANNOT_COL <- c(final = "#0077BB", nagalakshmi = "#EE7733", drs = "#009988")
END_LAB   <- c(`5prime` = "5' UTR", `3prime` = "3' UTR")

base_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10), plot.title = element_text(size = 11, face = "bold"),
        plot.tag = element_text(size = 13, face = "bold"))
fend <- function(x) factor(END_LAB[x], levels = END_LAB)

# -- Panel A. Both tails at 20 nt, three annotations --------------------------
rates <- read.delim(file.path(IN, "tifseq_rates.tsv"))
dfA <- rates %>%
  filter(tolerance_nt == 20) %>%
  select(annotation, end, n, Under = under_extension_pct, Concordant = concordant_pct,
         Over = over_extension_pct) %>%
  pivot_longer(c(Under, Concordant, Over), names_to = "cat", values_to = "pct") %>%
  mutate(cat = factor(cat, levels = c("Under", "Concordant", "Over")),
         annotation = factor(ANNOT_LAB[annotation], levels = ANNOT_LAB), end = fend(end))

cat_cols <- c(Under = COL_UNDER, Concordant = COL_CONC, Over = COL_OVER)
cat_full <- c(Under = "Under-extension: TIF-seq longer by 20 nt or more",
              Concordant = "Concordant: within 20 nt",
              Over = "Over-extension: TIF-seq shorter by 20 nt or more")

pA <- ggplot(dfA, aes(x = annotation, y = pct, fill = cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.8),
            vjust = -0.35, size = 2.7) +
  facet_wrap(~end) +
  scale_fill_manual(values = cat_cols, labels = cat_full, name = NULL, guide = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)), limits = c(0, 80)) +
  labs(title = "Both tails at a 20 nt tolerance", x = NULL, y = "Genes (%)", tag = "A") +
  base_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 8.5),
        legend.key.size = unit(0.4, "cm"), axis.text.x = element_text(size = 8.5))

# -- Panel B. Rates as a function of the tolerance -----------------------------
curve <- read.delim(file.path(IN, "tifseq_rate_curve.tsv"))
dfB <- curve %>%
  filter(tolerance_nt > 0) %>%
  select(annotation, end, tolerance_nt, Over = over_pct, Under = under_pct) %>%
  pivot_longer(c(Over, Under), names_to = "dir", values_to = "pct") %>%
  mutate(end = fend(end),
         dir = factor(dir, levels = c("Over", "Under"), labels = c("Over-extension", "Under-extension")),
         annotation = factor(ANNOT_LAB[annotation], levels = ANNOT_LAB))

pB <- ggplot(dfB, aes(tolerance_nt, pct, colour = annotation, linetype = dir)) +
  geom_vline(xintercept = 20, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 50, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~end) +
  scale_colour_manual(values = setNames(ANNOT_COL, ANNOT_LAB), name = NULL) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(breaks = c(0, 20, 50, 100, 150, 200)) +
  labs(title = "The tolerance is not load bearing", x = "Tolerance (nt)", y = "Genes (%)", tag = "B") +
  base_theme +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size = 8.5), legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.02, "cm"))

# -- Panel C. Dispersion of the benchmark around its own estimate -------------
rescur <- read.delim(file.path(IN, "tifseq_resolution_curve.tsv")) %>%
  rename(`5prime` = median_frac_reads_within_5prime, `3prime` = median_frac_reads_within_3prime) %>%
  pivot_longer(c(`5prime`, `3prime`), names_to = "end", values_to = "frac") %>%
  mutate(end = fend(end), pct = 100 * frac)
marks <- rescur %>% filter(tolerance_nt %in% c(20, 50, 100))

pC <- ggplot(rescur, aes(tolerance_nt, pct, colour = end)) +
  geom_vline(xintercept = 20, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 50, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  geom_line(linewidth = 0.8) +
  geom_point(data = marks, size = 1.6) +
  geom_text(data = marks, aes(label = sprintf("%.0f%%", pct)), hjust = -0.25, vjust = 2.1,
            size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c("#0077BB", "#CC3311"), name = NULL) +
  scale_x_continuous(breaks = c(0, 20, 50, 100, 150, 200)) +
  scale_y_continuous(limits = c(0, 102), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "TIF-seq dispersion around its own major isoform",
       x = "Distance from the gene's own major isoform (nt)",
       y = "Share of that gene's spanning reads (%)", tag = "C") +
  base_theme +
  theme(legend.position = c(0.82, 0.28),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"))

# -- Panel D. The envelope test ------------------------------------------------
env <- read.delim(file.path(IN, "tifseq_envelope.tsv"))
ENV_LAB <- c("95th percentile of reads" = "Within the 95th read percentile",
             "99th percentile of reads" = "Within the 99th read percentile",
             "most extreme observed" = "Within the most extreme isoform")
ENV_COL <- c("Within the 95th read percentile" = "#1B7837",
             "Within the 99th read percentile" = "#7FBC91",
             "Within the most extreme isoform" = "#D9F0D3")
dfD <- env %>% mutate(end = fend(end), envelope = factor(ENV_LAB[envelope], levels = rev(ENV_LAB)))

pD <- ggplot(dfD, aes(x = envelope, y = inside_pct, fill = envelope)) +
  geom_col(width = 0.68, colour = "grey35", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", inside_pct)), hjust = -0.15, size = 3.0) +
  facet_wrap(~end) + coord_flip() +
  scale_fill_manual(values = ENV_COL, guide = "none") +
  scale_y_continuous(limits = c(0, 118), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0))) +
  labs(title = "Released boundaries inside the observed TIF-seq isoform range",
       x = NULL, y = "Genes (%)", tag = "D") +
  base_theme + theme(axis.text.y = element_text(size = 8.5))

# -- Assemble ------------------------------------------------------------------
n5 <- rates$n[rates$annotation == "final" & rates$end == "5prime"][1]
n3 <- rates$n[rates$annotation == "final" & rates$end == "3prime"][1]
b5 <- round(100 - env$inside_pct[env$end == "5prime" & env$envelope == "most extreme observed"], 1)
b3 <- round(100 - env$inside_pct[env$end == "3prime" & env$envelope == "most extreme observed"], 1)
wrap_caption <- function(s, width = 165) paste(strwrap(gsub("[[:space:]]+", " ", s), width = width), collapse = "\n")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(caption = wrap_caption(paste0(
    "Over-extension and under-extension of transcript boundaries relative to the ",
    "independent TIF-seq benchmark (Pelechano et al., 2013), reported explicitly ",
    "in both directions. Differences are signed, annotation minus TIF-seq, over ",
    "genes where both report a non-zero UTR (final annotation: n = ", n5,
    " at the 5' end, n = ", n3, " at the 3' end). ",
    "(A) The prior reference over-extends relative to TIF-seq at almost the same ",
    "rate as the merged annotation, while the raw DRS calls agree with TIF-seq ",
    "most closely of the three. (B) Both rates vary smoothly with the tolerance ",
    "and the ordering of the three annotations does not change. ",
    "(C) The benchmark's own dispersion: a 20 nt tolerance is narrower than the ",
    "interquartile spread of TIF-seq ends within a single gene. ",
    "(D) Most released boundaries fall inside the range of ends TIF-seq observed ",
    "for that gene; ", b5, "% at the 5' end and ", b3,
    "% at the 3' end exceed every isoform TIF-seq recorded.")),
    theme = theme(plot.caption = element_text(size = 7.5, colour = "grey40", lineheight = 1.3, hjust = 0)))

w <- 12; h <- 10.5
base <- file.path(OUT, "suppl_fig6")
ggsave(paste0(base, ".pdf"), fig, width = w, height = h, device = "pdf")
ggsave(paste0(base, ".png"), fig, width = w, height = h, dpi = 300, device = "png")
cat("Saved:", paste0(base, ".pdf"), "\n")
cat("Saved:", paste0(base, ".png"), "\n")
