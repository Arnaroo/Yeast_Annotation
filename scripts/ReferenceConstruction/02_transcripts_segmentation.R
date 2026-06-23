# ============================================================
# 02_transcripts_segmentation.R
#
# Segment the cumulative DRS read-depth profile of each gene
# using change-point detection, then identify "expressed"
# (high-coverage) regions that define UTR boundaries.
#
# This script is run once per condition using the combined
# profile produced by 01_aggregate_data.R.
#
# Pipeline per gene:
#   1. perform_segmentation()  — Segmentor3IsBack (negative
#      binomial model) to detect all coverage transitions
#   2. select_segments()       — merge nearby breakpoints,
#      keep only those with a substantial coverage change
#   3. detect_expressed_regions() — label positions as
#      expressed if mean coverage exceeds a relative threshold
#
# Input:  {condition}_combined.tsv  (from 01_aggregate_data.R)
# Output: {condition}_expressed.tsv
#
# R packages required: dplyr, readr, Segmentor3IsBack
# Install Segmentor3IsBack from CRAN:
#   install.packages("Segmentor3IsBack")
# ============================================================

library(dplyr)
library(readr)
library(Segmentor3IsBack)

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these before running
# ════════════════════════════════════════════════════════════════════════

# Directory containing {condition}_combined.tsv (output of script 01)
DATA_DIR    <- "/path/to/data/gene_expressed"

# Directory for output files and optional diagnostic plots
RESULTS_DIR <- "/path/to/data/results"

# Condition to process (one at a time; run once per condition)
COND        <- "NS_PS"

# Minimum read count threshold: genes whose maximum coverage across
# all positions is below this value are excluded from segmentation.
MIN_READS   <- 20

# Set to TRUE to save diagnostic PNG plots for a random set of genes.
# Plots are saved to RESULTS_DIR/diagnostic_plots/.
SAVE_PLOTS  <- FALSE
N_PLOT_GENES <- 20    # number of randomly sampled genes to plot

# ════════════════════════════════════════════════════════════════════════

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
if (SAVE_PLOTS) {
  dir.create(file.path(RESULTS_DIR, "diagnostic_plots"), showWarnings = FALSE)
}

# ── Load and filter data ───────────────────────────────────────────────
message("Loading: ", file.path(DATA_DIR, paste0(COND, "_combined.tsv")))
db <- read_tsv(file.path(DATA_DIR, paste0(COND, "_combined.tsv")))

filtered_db <- db %>%
  group_by(gene) %>%
  filter(max(transcript_np) > MIN_READS) %>%
  ungroup() %>%
  mutate(breaks            = NA,
         interisting_breaks = NA,
         expressed          = NA)

message(sprintf("Genes after min-read filter (>%d): %d",
                MIN_READS, length(unique(filtered_db$gene))))

# Optional: sample genes for diagnostic plots
list_gene <- if (SAVE_PLOTS) {
  sample(unique(filtered_db$gene), min(N_PLOT_GENES, length(unique(filtered_db$gene))))
} else {
  character(0)
}

plot_dir <- file.path(RESULTS_DIR, "diagnostic_plots")

# ════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ════════════════════════════════════════════════════════════════════════

# ── Step 1: Segmentation ─────────────────────────────────────────────
# Runs Segmentor3IsBack on the per-position coverage vector.
# Returns all breakpoint positions (indices into data$position).
perform_segmentation <- function(data) {
  Sn     <- Segmentor3IsBack::Segmentor(data$transcript_np, model = 3)
  breaks <- floor(Segmentor3IsBack::getBreaks(Sn))

  if (SAVE_PLOTS && unique(data$gene) %in% list_gene) {
    png(file.path(plot_dir,
                  paste0("1_segmentation_", unique(data$gene), ".png")))
    plot(data$transcript_np,
         main = paste("Segmentation:", unique(data$gene)),
         xlab = "Position", ylab = "Read count")
    abline(v = breaks, col = "purple")
    dev.off()
  }
  breaks
}

# ── Step 2: Breakpoint selection ─────────────────────────────────────
# Merges breakpoints within 40 positions of each other (keeping the one
# with the larger coverage jump), then discards breakpoints where the
# coverage difference between adjacent segments is below 10 % of the
# gene's peak coverage.
select_segments <- function(data) {
  break_points <- data$position[which(data$breaks)]
  break_points <- break_points[-c(1, length(break_points))]

  # Merge nearby breakpoints iteratively
  repeat {
    filtered <- c()
    i <- 1
    changed <- FALSE

    while (i <= length(break_points)) {
      cur  <- break_points[i]
      if (i < length(break_points) && break_points[i + 1] - cur < 40) {
        nxt <- break_points[i + 1]
        if (nxt != length(data$gene)) {
          diff_cur <- abs(mean(data$transcript_np[max(1, cur - 25):cur]) -
                          mean(data$transcript_np[(cur + 1):min(nrow(data), cur + 25)]))
          diff_nxt <- abs(mean(data$transcript_np[max(1, nxt - 25):nxt]) -
                          mean(data$transcript_np[(nxt + 1):min(nrow(data), nxt + 25)]))
          filtered <- c(filtered, if (diff_cur >= diff_nxt) cur else nxt)
          i <- i + 2; changed <- TRUE
        } else break
      } else {
        filtered <- c(filtered, cur); i <- i + 1
      }
    }
    if (!changed) break
    break_points <- filtered
  }

  # Keep only breakpoints where the coverage change exceeds
  # 10 % of the gene's maximum coverage
  relative_threshold <- floor(0.1 * max(data$transcript_np))
  extended <- c(1, filtered, length(data$transcript_np))

  interesting <- sapply(seq_along(filtered), function(j) {
    lm <- mean(data$transcript_np[(extended[j] + 1):filtered[j]])
    rm <- mean(data$transcript_np[(filtered[j] + 1):extended[j + 2]])
    abs(lm - rm) > relative_threshold
  })
  interesting_breaks <- filtered[interesting]

  if (SAVE_PLOTS && unique(data$gene) %in% list_gene) {
    png(file.path(plot_dir,
                  paste0("2_selected_", unique(data$gene), ".png")))
    plot(data$transcript_np,
         main = paste("Selected breakpoints:", unique(data$gene)),
         xlab = "Position", ylab = "Read count")
    abline(v = interesting_breaks, col = "red")
    dev.off()
  }
  interesting_breaks
}

# ── Step 3: Expressed region detection ───────────────────────────────
# Marks positions as expressed if the segment mean exceeds 5 % of
# the gene's peak coverage.  Boundary segments (outside the central
# window 1000–(n-1000)) are propagated inward if they exceed 10 %
# of peak coverage.
detect_expressed_regions <- function(data) {
  data$expressed <- FALSE
  borne_inf <- 1000
  borne_sup <- nrow(data) - 1000

  break_points        <- data$position[which(data$breaks)]
  interesting_breaks  <- data$position[which(data$interisting_breaks)]

  middle_breaks <- interesting_breaks[
    interesting_breaks >= borne_inf & interesting_breaks <= borne_sup
  ]
  if (length(middle_breaks) == 0) return(NULL)

  middle_breaks <- sort(c(middle_breaks, borne_inf, borne_sup))
  relative_value <- min(floor(0.05 * max(data$transcript_np)), 100)

  for (i in 2:length(middle_breaks)) {
    seg_mean <- mean(data$transcript_np[middle_breaks[i-1]:middle_breaks[i]], na.rm = TRUE)
    if (seg_mean > relative_value)
      data$expressed[middle_breaks[i-1]:middle_breaks[i]] <- TRUE
  }

  left_flag  <- isTRUE(data$expressed[which(data$position == borne_inf)])
  right_flag <- isTRUE(data$expressed[which(data$position == borne_sup)])

  # Propagate expressed regions outward (if coverage remains above 10 % of peak)
  rel10 <- min(floor(0.1 * max(data$transcript_np)), 100)

  if (left_flag) {
    first_bp <- sort(c(break_points[break_points <= borne_inf], 1, borne_inf))
    old_mean <- mean(data$transcript_np[middle_breaks[1]:middle_breaks[2]], na.rm = TRUE)
    for (i in (length(first_bp) - 1):1) {
      seg_mean <- mean(data$transcript_np[first_bp[i]:first_bp[i+1]], na.rm = TRUE)
      if (seg_mean > rel10 && old_mean >= seg_mean) {
        data$expressed[first_bp[i]:first_bp[i+1]] <- TRUE
        old_mean <- seg_mean
      } else break
    }
  }

  if (right_flag) {
    final_bp <- sort(c(break_points[break_points >= borne_sup], nrow(data), borne_sup))
    old_mean <- mean(data$transcript_np[middle_breaks[length(middle_breaks)-1]:
                                        middle_breaks[length(middle_breaks)]], na.rm = TRUE)
    for (i in 1:(length(final_bp) - 1)) {
      seg_mean <- mean(data$transcript_np[final_bp[i]:final_bp[i+1]], na.rm = TRUE)
      if (seg_mean > rel10 && old_mean <= seg_mean) {
        data$expressed[final_bp[i]:final_bp[i+1]] <- TRUE
        old_mean <- seg_mean
      } else break
    }
  }

  if (SAVE_PLOTS && unique(data$gene) %in% list_gene) {
    all_bp <- sort(c(1, borne_inf, borne_sup, break_points, nrow(data)))
    png(file.path(plot_dir,
                  paste0("3_expressed_", unique(data$gene), ".png")))
    plot(data$transcript_np,
         main = paste("Expressed regions:", unique(data$gene)),
         xlab = "Position", ylab = "Read count")
    abline(v = c(borne_inf, borne_sup), col = "red")
    for (i in 2:length(all_bp)) {
      s <- all_bp[i-1]; e <- all_bp[i]
      if (all(data$expressed[s:e]))
        rect(s, min(data$transcript_np), e, max(data$transcript_np),
             col = rgb(0, 1, 0, 0.5), border = NA)
    }
    dev.off()
  }
  data
}

# ════════════════════════════════════════════════════════════════════════
# MAIN LOOP — process each gene
# ════════════════════════════════════════════════════════════════════════
unique_genes <- unique(filtered_db$gene)
message(sprintf("Processing %d genes for condition: %s", length(unique_genes), COND))

all_db <- data.frame()
for (gene_name in unique_genes) {
  message("  ", gene_name)
  gene_data <- filtered_db[filtered_db$gene == gene_name, ]
  gene_data$expressed <- NA

  breaks     <- perform_segmentation(gene_data)
  breaks     <- unique(c(breaks, 1, nrow(gene_data)))
  gene_data$breaks[gene_data$position %in% breaks] <- TRUE

  best_breaks <- select_segments(gene_data)
  best_breaks <- unique(c(best_breaks, 1, nrow(gene_data)))
  gene_data$interisting_breaks[gene_data$position %in% best_breaks] <- TRUE

  result <- detect_expressed_regions(gene_data)
  if (!is.null(result)) all_db <- rbind(all_db, result)
}

# ── Save output ────────────────────────────────────────────────────────
out_path <- file.path(RESULTS_DIR, paste0(COND, "_expressed.tsv"))
write_tsv(all_db, out_path)
message(sprintf("Written: %s  (%d rows, %d genes)",
                out_path, nrow(all_db), length(unique(all_db$gene))))
