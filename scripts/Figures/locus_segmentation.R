# ============================================================
# locus_segmentation.R
# ============================================================
# Segments the real coverage locus_coverage_real.py built for each of
# the nine Figure 1B / Supplementary Fig. S3 windows, using exactly the
# released boundary-calling algorithm (segmentation_corrected.R).
#
# Segments on `window_position` (1..n in the gene's own transcript
# direction), not on genomic order: for a minus-strand gene the two run
# in opposite directions, and segmenting on genomic order silently
# swaps which side of the call is five and which is three.
#
# Input:   locus_coverage.tsv (locus_coverage_real.py)
#          ../Sensitivity/segmentation_corrected.R (sourced, not copied)
# Outputs: locus_breakpoints.tsv
#          locus_expressed.tsv
#
# Usage:
#   Rscript locus_segmentation.R --indir Data/LocusTracks --outdir Data/LocusTracks
# ============================================================

suppressMessages({
  library(dplyr)
  library(readr)
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

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(HERE, "..", "Sensitivity", "segmentation_corrected.R"))

# Copied verbatim from ReferenceConstruction/02_transcripts_segmentation.R.
perform_segmentation <- function(data) {
  Sn <- Segmentor3IsBack::Segmentor(data$transcript_np, model = 3)
  floor(Segmentor3IsBack::getBreaks(Sn))
}

cov <- read_tsv(file.path(IN, "locus_coverage.tsv"), show_col_types = FALSE)

bp_rows <- list()
exp_rows <- list()

for (pn in unique(cov$panel)) {
  w <- cov %>%
    filter(panel == pn) %>%
    arrange(window_position) %>%
    mutate(position = row_number(), transcript_np = as.integer(depth_sense)) %>%
    as.data.frame()
  target_gene <- unique(w$gene)[1]

  brk <- perform_segmentation(w)
  brk <- unique(c(brk, 1, nrow(w)))
  brk <- brk[brk >= 1 & brk <= nrow(w)]

  best <- select_segments(w$transcript_np, brk, REL_MERGE_WINDOW, REL_EFFECT_FRAC)
  best <- unique(c(best, 1, nrow(w)))

  ex <- detect_expressed_region_corrected(w$transcript_np, brk, best, target_gene)

  gpos <- (cov %>% filter(panel == pn) %>% arrange(window_position))$position

  bp_rows[[pn]] <- bind_rows(
    tibble(panel = pn, stage = "candidate", genomic_position = gpos[sort(brk)]),
    tibble(panel = pn, stage = "selected", genomic_position = gpos[sort(best)])
  )

  if (!is.null(ex)) {
    r <- rle(ex)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    keep <- which(r$values)
    if (length(keep)) {
      exp_rows[[pn]] <- tibble(
        panel = pn,
        expressed_start = pmin(gpos[starts[keep]], gpos[ends[keep]]),
        expressed_end   = pmax(gpos[starts[keep]], gpos[ends[keep]])
      )
    }
  }

  cat(sprintf("%-7s window %5d nt  candidate %4d  selected %3d  expressed segments %d\n",
              pn, nrow(w), length(brk), length(best),
              if (is.null(exp_rows[[pn]])) 0L else nrow(exp_rows[[pn]])))
}

bp  <- bind_rows(bp_rows)
exp <- bind_rows(exp_rows)
write_tsv(bp,  file.path(OUT, "locus_breakpoints.tsv"))
write_tsv(exp, file.path(OUT, "locus_expressed.tsv"))
cat("\nwrote locus_breakpoints.tsv (", nrow(bp), "rows )\n")
cat("wrote locus_expressed.tsv (", nrow(exp), "rows )\n")
