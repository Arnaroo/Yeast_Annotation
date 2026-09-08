#!/usr/bin/env Rscript
# =============================================================================
# plot_suppl_fig7.R
# =============================================================================
# Supplementary Figure S7. Where the benchmark falls on the genes whose
# DRS-derived UTR is shorter than the prior reference, and the mechanism
# behind them.
#
# Panels
#   A  What TIF-seq says about each shortfall gene, three ways: whether
#      the benchmark sits inside both calls, between them, or beyond
#      both. The third class is the only one where DRS demonstrably
#      lost sequence.
#   B  The control. On genes where our call is longer, not shorter, the
#      benchmark switches sides -- so panel A is not TIF-seq preferring
#      the shorter number.
#   C  The mechanism, from the annotation: odds of the beyond-both
#      class on transcript size and on depth. Steep on size at the 5'
#      end, flat at the 3' end -- the processivity signature.
#   D  The mechanism, from the alignments: fraction of reads reaching
#      each annotated end across deciles of transcript length. No
#      annotation quantity enters this panel.
#
# Inputs are the outputs of shortfall.py and read_geometry.py. Nothing
# is recomputed here.
#
# Usage:
#   Rscript plot_suppl_fig7.R --indir Data/Shortfall --outdir Figures_out
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

COL_BELOW <- "#228833"; COL_BTWN <- "#AAAAAA"; COL_ABOVE <- "#AA3377"
COL_DRS <- "#009988"; COL_PRIOR <- "#EE7733"

END_LAB <- c(`5prime` = "5' UTR", `3prime` = "3' UTR", `5'` = "5' UTR", `3'` = "3' UTR")
fend <- function(x) factor(unname(END_LAB[as.character(x)]), levels = c("5' UTR", "3' UTR"))

base_theme <- theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey92", colour = "grey60"),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 10), plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.8, colour = "grey30"),
        plot.tag = element_text(size = 13, face = "bold"))

arb <- read.delim(file.path(IN, "shortfall_arbitration.tsv"))

# -- Panel A. Three-way arbitration of the shortfall ---------------------------
CAT_LAB <- c(below = "TIF-seq inside both calls\nboth annotations too long",
             between = "TIF-seq between the calls\nDRS too short, prior too long",
             above = "TIF-seq beyond both calls\nsequence genuinely lost")
dfA <- arb %>% filter(subset == "drs_shorter") %>%
  select(end, n, below = tifseq_below_both_pct, between = tifseq_between_pct, above = tifseq_above_both_pct) %>%
  pivot_longer(c(below, between, above), names_to = "cat", values_to = "pct") %>%
  mutate(cat = factor(CAT_LAB[cat], levels = CAT_LAB), end = fend(end))
nA <- arb %>% filter(subset == "drs_shorter") %>% select(end, n) %>% mutate(end = fend(end))

pA <- ggplot(dfA, aes(x = end, y = pct, fill = cat)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  geom_text(data = nA, aes(x = end, y = 68, label = paste0("n = ", format(n, big.mark = ","))),
            inherit.aes = FALSE, size = 2.9, colour = "grey30") +
  scale_fill_manual(values = setNames(c(COL_BELOW, COL_BTWN, COL_ABOVE), CAT_LAB), name = NULL,
                    guide = guide_legend(nrow = 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, 72)) +
  labs(title = "Where the benchmark falls on shortfall genes",
       subtitle = "Genes where the DRS call is shorter than the prior reference",
       x = NULL, y = "Shortfall genes (%)", tag = "A") +
  base_theme + theme(legend.position = "bottom", legend.text = element_text(size = 7.6),
                     legend.key.size = unit(0.42, "cm"))

# -- Panel B. The control, genes where the DRS call is longer -----------------
dfB <- arb %>% select(end, subset, DRS = drs_closer_pct, Prior = prior_closer_pct) %>%
  pivot_longer(c(DRS, Prior), names_to = "closer", values_to = "pct") %>%
  mutate(end = fend(end),
         subset = factor(subset, levels = c("drs_shorter", "drs_longer_control"),
                         labels = c("DRS call shorter\n(the shortfall)", "DRS call longer\n(control)")),
         closer = factor(closer, levels = c("DRS", "Prior"),
                         labels = c("DRS call closer to TIF-seq", "Prior reference closer to TIF-seq")))

pB <- ggplot(dfB, aes(x = subset, y = pct, fill = closer)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.8), vjust = -0.4, size = 2.8) +
  geom_hline(yintercept = 50, colour = "grey55", linewidth = 0.4, linetype = "dotted") +
  facet_wrap(~end) +
  scale_fill_manual(values = c(COL_DRS, COL_PRIOR), name = NULL, guide = guide_legend(nrow = 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, 92)) +
  labs(title = "The benchmark is not simply preferring the shorter call",
       subtitle = "When the DRS call is longer instead, TIF-seq switches sides",
       x = NULL, y = "Genes (%)", tag = "B") +
  base_theme + theme(legend.position = "bottom", legend.text = element_text(size = 8),
                     legend.key.size = unit(0.42, "cm"), axis.text.x = element_text(size = 8.2))

# -- Panel C. Odds of the beyond-both class, on size and on depth -------------
TERM_LAB <- c(log10_cds_len = "Coding length\nper 10 fold", log10_depth = "Read depth\nper 10 fold")
dfC <- read.delim(file.path(IN, "shortfall_models.tsv")) %>%
  filter(model == "P(both too short | shortfall)") %>%
  mutate(end = fend(end), term = factor(TERM_LAB[term], levels = rev(TERM_LAB)))

pC <- ggplot(dfC, aes(x = odds_ratio_per_decade, y = term, colour = end)) +
  geom_vline(xintercept = 1, colour = "grey45", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, linewidth = 0.6,
                position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  scale_x_continuous(trans = "log10", breaks = c(0.25, 0.5, 1, 2, 4, 8),
                     labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  scale_colour_manual(values = c("5' UTR" = COL_ABOVE, "3' UTR" = "#0077BB"), name = NULL) +
  labs(title = "The mechanism, from the annotation",
       subtitle = "Odds that TIF-seq lies beyond both calls. Steep on size at the 5' end only.",
       x = "Odds ratio (log scale)", y = NULL, tag = "C") +
  base_theme + theme(legend.position = "bottom", legend.key.size = unit(0.42, "cm"),
                     axis.text.y = element_text(size = 8.4))

# -- Panel D. Reads reaching each end, across transcript length ---------------
dfD <- read.delim(file.path(IN, "read_geometry_strata.tsv")) %>%
  filter(stratum == "tx_len") %>%
  select(bin, lo, hi, n, `5' UTR` = mean_frac_reach5, `3' UTR` = mean_frac_reach3) %>%
  pivot_longer(c(`5' UTR`, `3' UTR`), names_to = "end", values_to = "frac") %>%
  mutate(end = factor(end, levels = c("5' UTR", "3' UTR")), mid = sqrt(lo * hi))

pD <- ggplot(dfD, aes(mid, 100 * frac, colour = end)) +
  geom_line(linewidth = 0.75) + geom_point(size = 2.1) +
  scale_x_continuous(trans = "log10", breaks = c(500, 1000, 2000, 4000, 8000)) +
  scale_y_continuous(limits = c(0, 45), expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = c("5' UTR" = COL_ABOVE, "3' UTR" = "#0077BB"), name = NULL) +
  labs(title = "The mechanism, from the alignments",
       subtitle = "Reads reaching within 10 nt of the annotated end, by decile of transcript length",
       x = "Transcript length (nt, log scale)", y = "Reads reaching the end (%)", tag = "D") +
  base_theme + theme(legend.position = "bottom", legend.key.size = unit(0.42, "cm"))

# -- Assemble ------------------------------------------------------------------
g5 <- arb$tifseq_above_both_pct[arb$end == "5prime" & arb$subset == "drs_shorter"]
g3 <- arb$tifseq_above_both_pct[arb$end == "3prime" & arb$subset == "drs_shorter"]
d5 <- arb$drs_closer_pct[arb$end == "5prime" & arb$subset == "drs_shorter"]
d3 <- arb$drs_closer_pct[arb$end == "3prime" & arb$subset == "drs_shorter"]
n5 <- format(arb$n[arb$end == "5prime" & arb$subset == "drs_shorter"], big.mark = ",")
n3 <- format(arb$n[arb$end == "3prime" & arb$subset == "drs_shorter"], big.mark = ",")

fig <- (pA | pB) / (pC | pD) +
  plot_annotation(caption = paste0(
    "Genes where the DRS-derived UTR is shorter than the prior reference, tested against the ",
    "independent TIF-seq benchmark (Pelechano et al., 2013), which took no part in constructing ",
    "the annotation (n = ", n5, " at the 5' end, n = ", n3, " at the 3' end). ",
    "(A) On most shortfall genes the benchmark is shorter than both calls, so the shorter DRS call ",
    "is the better of the two; the DRS call is closer to TIF-seq for ", sprintf("%.1f", d5),
    "% of genes at the 5' end and ", sprintf("%.1f", d3), "% at the 3' end. Sequence is demonstrably ",
    "lost only in the third class, ", sprintf("%.1f", g5), "% of shortfall genes at the 5' end and ",
    sprintf("%.1f", g3), "% at the 3' end. ",
    "(B) The direction is not an artefact of the comparison: on genes where the DRS call is longer, ",
    "the prior reference is the closer of the two. ",
    "(C) Odds of the third class rise steeply with coding length at the 5' end and not at all at the ",
    "3' end, which is where a 3' to 5' processivity limit predicts them to rise and not to rise. ",
    "(D) The same contrast measured directly from the alignments, with no annotation quantity ",
    "involved: the fraction of reads reaching the annotated 5' end collapses as transcripts get ",
    "longer, while the fraction reaching the 3' end does not."),
    theme = theme(plot.caption = element_text(size = 7.5, colour = "grey40", lineheight = 1.3, hjust = 0)))

w <- 12; h <- 10.5
base <- file.path(OUT, "suppl_fig7")
ggsave(paste0(base, ".pdf"), fig, width = w, height = h, device = "pdf")
ggsave(paste0(base, ".png"), fig, width = w, height = h, dpi = 300, device = "png")
cat("Saved:", paste0(base, ".pdf"), "\n")
cat("Saved:", paste0(base, ".png"), "\n")
