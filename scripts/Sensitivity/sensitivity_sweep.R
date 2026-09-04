# ============================================================
# sensitivity_sweep.R
# ============================================================
# Threshold sensitivity of the boundary-calling step (segmentation_corrected.R),
# measured on the six real construction conditions rather than a proxy
# library, answering the reviewer question of why the four segmentation
# thresholds take the values they do.
#
# Four thresholds govern boundary calling downstream of the raw
# Segmentor3IsBack change-point search: the coverage floor gating which
# genes are attempted (min_reads), the breakpoint-merging window
# (merge_window), the breakpoint effect-size filter (effect_frac), and
# the expressed-segment floor (expr_frac). Each is swept across five
# values spanning well beyond the released setting on either side, for
# 5^4 = 625 combinations per condition (the released combination is an
# exact grid point). Segmentor3IsBack itself takes no threshold
# argument, so the change-point search is run once per gene per
# condition -- over every gene that could qualify at any grid point --
# and the 625 combinations then only re-evaluate the cheap downstream
# functions in segmentation_corrected.R, which is all four thresholds
# actually change.
#
# The released annotation's DRS-derived boundary (DRS_UTR_corrected.tsv)
# is the per-gene maximum across all six conditions, so this script's
# output is the same per-condition quantity that gets combined that way
# downstream (see sensitivity_analyse.py), not a stand-in for it.
#
# Input:  {condition}_combined.tsv for the six conditions NS_PS, NS_RDT,
#         S10_PS, S10_RDT, NS and S10, in the same format produced by
#         ReferenceConstruction/01_aggregate_data.R (columns: gene,
#         position, transcript_np)
# Output: sensitivity_sweep_summary.tsv    one row per condition x combo
#         sensitivity_sweep_per_gene.tsv.gz  one row per condition x combo x gene
#
# Run:  Rscript sensitivity_sweep.R
# ============================================================

suppressPackageStartupMessages({
  library(Segmentor3IsBack)
  library(parallel)
  library(readr)
  library(dplyr)
})

# -- CONFIGURATION -------------------------------------------------------
HERE       <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "."
COMBINED_DIR <- file.path(HERE, "..", "..", "Data", "CombinedProfiles")  # {condition}_combined.tsv
WORK_DIR     <- file.path(HERE, "work")   # breakpoint caches, resumable
OUT_DIR      <- file.path(HERE, "..", "..", "Data", "Sensitivity")
THREADS      <- 16
# -------------------------------------------------------------------------

dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(HERE, "segmentation_corrected.R"))

MIN_READS_SET    <- c(5, 10, 20, 40, 80)
MERGE_WINDOW_SET <- c(10, 20, 40, 80, 160)
EFFECT_FRAC_SET  <- c(0.025, 0.05, 0.10, 0.20, 0.35)
EXPR_FRAC_SET    <- c(0.01, 0.025, 0.05, 0.10, 0.20)
# The released values (min_reads=20, merge_window=40, effect_frac=0.10,
# expr_frac=0.05) sit at index 3 of every set above, so the released
# combination is an exact grid point.

SIX <- c(NS_PS = "NS_PS_combined.tsv", NS_RDT = "NS_RDT_combined.tsv",
          S10_PS = "S10_PS_combined.tsv", S10_RDT = "S10_RDT_combined.tsv",
          NS = "NS_combined.tsv", S10 = "S10_combined.tsv")

say <- function(...) { cat(sprintf(...), "\n", sep = ""); flush.console() }

say("sensitivity_sweep.R")
say("started %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
say("")

read_gp <- function(path) {
  read_tsv(path, col_types = cols(gene = col_character(), position = col_integer(),
                                   transcript_np = col_integer()), progress = FALSE)
}

run_combo <- function(min_reads, merge_window, effect_frac, expr_frac,
                       cand, gene_max, breaks_cache, cov_list) {
  genes <- cand[gene_max[cand] > min_reads]
  res <- mclapply(genes, function(g) {
    ab <- breaks_cache[[g]]
    if (is.null(ab)) return(NULL)
    np <- cov_list[[g]]
    ib <- select_segments(np, ab, merge_window, effect_frac)
    ex <- detect_expressed_region_corrected(np, ab, ib, g, expr_frac)
    if (is.null(ex)) return(NULL)
    u <- compute_utrs(ex, length(np))
    c(five = unname(u["five"]), three = unname(u["three"]))
  }, mc.cores = THREADS)
  ok <- !vapply(res, is.null, logical(1))
  if (!any(ok)) return(NULL)
  data.frame(gene = genes[ok],
             five = vapply(res[ok], function(x) x[["five"]], numeric(1)),
             three = vapply(res[ok], function(x) x[["three"]], numeric(1)),
             stringsAsFactors = FALSE)
}

grid <- expand.grid(min_reads = MIN_READS_SET, merge_window = MERGE_WINDOW_SET,
                    effect_frac = EFFECT_FRAC_SET, expr_frac = EXPR_FRAC_SET,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
say("grid: %d combinations per condition", nrow(grid))
say("")

all_summary <- list()
all_per_gene <- list()

for (cond in names(SIX)) {
  say("======================================================")
  say("CONDITION %s", cond)
  path <- file.path(COMBINED_DIR, SIX[[cond]])
  if (!file.exists(path)) { say("  SKIP, missing %s", path); next }

  db <- read_gp(path)
  cov_list <- split(db$transcript_np, db$gene)
  rm(db); invisible(gc())
  gene_max <- vapply(cov_list, max, integer(1))
  cand <- names(cov_list)[gene_max > min(MIN_READS_SET)]
  say("  genes total %d, candidates (peak > %d) %d", length(cov_list), min(MIN_READS_SET), length(cand))

  # Segmentor3IsBack runs once per gene per condition, over the union of
  # genes that could ever qualify, and is cached so a killed run resumes
  # without repeating a finished condition.
  cache_path <- file.path(WORK_DIR, sprintf("%s_breaks_cache.rds", cond))
  if (file.exists(cache_path)) {
    breaks_cache <- readRDS(cache_path)
    say("  reusing breakpoint cache: %s", cache_path)
  } else {
    t0 <- Sys.time()
    breaks_cache <- mclapply(cand, function(g) {
      np <- cov_list[[g]]
      tryCatch({
        Sn <- Segmentor3IsBack::Segmentor(np, model = 3)
        floor(Segmentor3IsBack::getBreaks(Sn))
      }, error = function(e) NULL)
    }, mc.cores = THREADS, mc.preschedule = FALSE)
    names(breaks_cache) <- cand
    saveRDS(breaks_cache, cache_path)
    say("  change-point search done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  nfail <- sum(vapply(breaks_cache, is.null, logical(1)))
  say("  genes segmented %d, failed %d", length(breaks_cache) - nfail, nfail)

  t0 <- Sys.time()
  per_gene <- vector("list", nrow(grid))
  summ     <- vector("list", nrow(grid))
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    d <- run_combo(g$min_reads, g$merge_window, g$effect_frac, g$expr_frac,
                   cand, gene_max, breaks_cache, cov_list)
    if (is.null(d)) next
    d$combo <- k
    per_gene[[k]] <- d
    summ[[k]] <- data.frame(
      condition = cond, combo = k,
      min_reads = g$min_reads, merge_window = g$merge_window,
      effect_frac = g$effect_frac, expr_frac = g$expr_frac,
      n_genes = nrow(d),
      median_five = stats::median(d$five), median_three = stats::median(d$three),
      stringsAsFactors = FALSE)
  }
  say("  grid evaluated in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

  summary_df <- do.call(rbind, summ)
  gene_df <- do.call(rbind, per_gene)
  gene_df$condition <- cond
  all_summary[[cond]] <- summary_df
  all_per_gene[[cond]] <- gene_df
  say("  combos with output: %d / %d", nrow(summary_df), nrow(grid))
  say("")
}

summary_all <- bind_rows(all_summary)
per_gene_all <- bind_rows(all_per_gene) %>% select(condition, combo, gene, five, three)

write.table(summary_all, file.path(OUT_DIR, "sensitivity_sweep_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
gz <- gzfile(file.path(OUT_DIR, "sensitivity_sweep_per_gene.tsv.gz"), "w")
write.table(per_gene_all, gz, sep = "\t", row.names = FALSE, quote = FALSE)
close(gz)

say("======================================================")
say("wrote %s (%d rows)", file.path(OUT_DIR, "sensitivity_sweep_summary.tsv"), nrow(summary_all))
say("wrote %s (%d rows)", file.path(OUT_DIR, "sensitivity_sweep_per_gene.tsv.gz"), nrow(per_gene_all))
say("finished %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
