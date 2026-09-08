# ============================================================
# segmentation_corrected.R
# ============================================================
# The boundary-calling algorithm used to build the released annotation,
# factored out into standalone functions so it can be re-run with
# different threshold values without touching ReferenceConstruction/.
#
# This is the same algorithm as
# ReferenceConstruction/02_transcripts_segmentation.R, confirmed to
# reproduce the released final_utr.tsv exactly (both ends, every
# comparable gene, zero deviation). Any script that needs to compute a
# boundary call under a non-released set of thresholds -- as
# sensitivity_sweep.R does -- sources this file rather than
# re-implementing the logic, so there is exactly one place it can drift.
#
# Four thresholds are exposed as arguments to the two functions below:
#   min_reads     gates which genes are attempted at all (checked by the
#                 caller, against each gene's peak coverage)
#   merge_window  breakpoints closer together than this are merged
#   effect_frac   a merged breakpoint survives only if the mean-coverage
#                 step across it exceeds this fraction of the gene's peak
#   expr_frac     a segment counts as expressed if its mean coverage
#                 exceeds this fraction of the peak, capped at CAP reads
#
# Two further constants are part of the algorithm rather than the four
# reviewed thresholds, and are held fixed regardless of which of the
# four is being swept:
#   CAP            the expr_frac / rel10 relative thresholds are capped
#                  at this many reads, so a very highly expressed gene
#                  does not require a proportionally huge segment mean
#                  to be called expressed
#   rel10          the boundary-propagation threshold, fixed at
#                  floor(0.1 * peak) capped at CAP; governs how far an
#                  expressed call is allowed to extend outward one
#                  breakpoint at a time once the core window qualifies
#
# A handful of genes needed a hand-set expression threshold or are
# always called unexpressed at their annotated position, both carried
# over unchanged from the released construction scripts:
#   SPECIAL_THRESH     gene -> fixed relative_value, overriding expr_frac
#   FORCE_UNEXPRESSED  genes whose profile is expression-like but is
#                      known, from manual inspection, not to reflect
#                      the annotated transcript (e.g. a neighbouring
#                      gene's read-through)
#
# Provides: select_segments(), detect_expressed_region_corrected(),
# compute_utrs(), and the constants FLANK, CAP, SPECIAL_THRESH,
# FORCE_UNEXPRESSED.
# ============================================================

FLANK <- 1000
CAP <- 50
REL_MERGE_WINDOW <- 40   # released value of the merge_window threshold
REL_EFFECT_FRAC <- 0.10  # released value of the effect_frac threshold
SPECIAL_THRESH <- c(YAL003W = 250, YDR500C = 250)
FORCE_UNEXPRESSED <- c("YCL058C", "YBL094C")

# Merges Segmentor3IsBack breakpoints closer than merge_window, keeping
# whichever side of a merged pair shows the larger local coverage step,
# then discards any surviving breakpoint whose step is below
# effect_frac of the gene's peak coverage.
select_segments <- function(np, all_breaks, merge_window, effect_frac) {
  n <- length(np)
  all_breaks <- all_breaks[all_breaks >= 1L & all_breaks <= n]
  bp <- sort(unique(c(all_breaks, 1L, n)))
  if (length(bp) <= 2L) return(integer(0))
  bp <- bp[-c(1L, length(bp))]
  filtered <- integer(0)
  repeat {
    filtered <- integer(0); i <- 1L; changed <- FALSE
    while (i <= length(bp)) {
      cur <- bp[i]
      if (i < length(bp) && bp[i + 1L] - cur < merge_window) {
        nxt <- bp[i + 1L]
        if (nxt != n) {
          diff_cur <- abs(mean(np[max(1L, cur - 25L):cur]) - mean(np[(cur + 1L):min(n, cur + 25L)]))
          diff_nxt <- abs(mean(np[max(1L, nxt - 25L):nxt]) - mean(np[(nxt + 1L):min(n, nxt + 25L)]))
          filtered <- c(filtered, if (diff_cur >= diff_nxt) cur else nxt)
          i <- i + 2L; changed <- TRUE
        } else break
      } else { filtered <- c(filtered, cur); i <- i + 1L }
    }
    if (!changed) break
    bp <- filtered
  }
  if (length(filtered) == 0L) return(integer(0))
  relative_threshold <- floor(effect_frac * max(np))
  extended <- c(1L, filtered, n)
  keep <- vapply(seq_along(filtered), function(j) {
    lm <- mean(np[(extended[j] + 1L):filtered[j]])
    rm <- mean(np[(filtered[j] + 1L):extended[j + 2L]])
    isTRUE(abs(lm - rm) > relative_threshold)
  }, logical(1))
  filtered[keep]
}

# Marks the expressed region of a gene's coverage profile: the segments
# between the surviving breakpoints (select_segments()) whose mean
# coverage exceeds expr_frac of peak, extended outward one breakpoint at
# a time wherever coverage stays above the fixed rel10 threshold. Returns
# NULL for genes whose ORF-flanking window is too short to carry a full
# FLANK on both sides (see borne_inf/borne_sup below).
detect_expressed_region_corrected <- function(np, all_breaks, interesting, gene, expr_frac = 0.05) {
  n <- length(np)
  expressed <- logical(n)
  borne_inf <- FLANK; borne_sup <- n - FLANK
  if (borne_sup <= borne_inf) return(NULL)
  all_breaks <- all_breaks[all_breaks >= 1L & all_breaks <= n]
  break_points <- sort(unique(c(all_breaks, 1L, n)))
  ib <- sort(unique(c(interesting, 1L, n)))
  middle <- ib[ib >= borne_inf & ib <= borne_sup]
  if (length(middle) == 0L) return(NULL)
  middle <- sort(c(middle, borne_inf, borne_sup))

  if (gene %in% names(SPECIAL_THRESH)) {
    relative_value <- SPECIAL_THRESH[[gene]]
  } else {
    relative_value <- min(floor(expr_frac * max(np)), CAP)
  }
  for (i in 2:length(middle)) {
    seg_mean <- mean(np[middle[i - 1L]:middle[i]], na.rm = TRUE)
    if (seg_mean > relative_value) expressed[middle[i - 1L]:middle[i]] <- TRUE
  }

  left_flag <- isTRUE(expressed[borne_inf])
  right_flag <- isTRUE(expressed[borne_sup])

  if (left_flag) {
    first_bp <- sort(c(break_points[break_points <= borne_inf], 1L, borne_inf))
    expressed[first_bp[length(first_bp) - 1L]:first_bp[length(first_bp)]] <- TRUE
  }
  if (right_flag) {
    final_bp <- sort(c(break_points[break_points >= borne_sup], n, borne_sup))
    expressed[final_bp[1L]:final_bp[2L]] <- TRUE
  }

  rel10 <- min(floor(0.1 * max(np)), CAP)
  if (left_flag) {
    old_mean <- mean(np[middle[1L]:middle[2L]], na.rm = TRUE)
    first_bp <- sort(c(break_points[break_points <= borne_inf], 1L))
    if (length(first_bp) >= 2L) {
      for (i in (length(first_bp) - 1L):1L) {
        seg_mean <- mean(np[first_bp[i]:first_bp[i + 1L]], na.rm = TRUE)
        if (seg_mean > rel10 && old_mean >= seg_mean) {
          expressed[first_bp[i]:first_bp[i + 1L]] <- TRUE; old_mean <- seg_mean
        } else break
      }
    }
  }
  if (right_flag) {
    old_mean <- mean(np[middle[length(middle) - 1L]:middle[length(middle)]], na.rm = TRUE)
    final_bp <- sort(c(break_points[break_points >= borne_sup], n, borne_sup))
    if (length(final_bp) >= 2L) {
      for (i in 1L:(length(final_bp) - 1L)) {
        seg_mean <- mean(np[final_bp[i]:final_bp[i + 1L]], na.rm = TRUE)
        if (seg_mean > rel10 && (old_mean - rel10) <= seg_mean) {
          expressed[final_bp[i]:final_bp[i + 1L]] <- TRUE; old_mean <- seg_mean
        } else break
      }
    }
  }

  if (gene %in% FORCE_UNEXPRESSED) expressed[] <- FALSE
  if (isTRUE(expressed[1])) expressed[1:borne_inf] <- FALSE
  if (isTRUE(expressed[n])) expressed[borne_sup:n] <- FALSE
  expressed
}

# Converts an expressed-region indicator to 5'/3' UTR lengths relative
# to the FLANK-nt window flanking the annotated ORF.
compute_utrs <- function(expressed, n) {
  pos <- which(expressed)
  if (length(pos) == 0L) return(c(five = 0L, three = 0L))
  c(five  = max(0L, FLANK - min(pos)),
    three = max(0L, max(pos) - (n - FLANK)))
}
