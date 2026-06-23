# ============================================================
# 03_utr_table.R
#
# From the expressed-region tables produced by script 02,
# compute per-gene UTR lengths, then merge with the Nagalakshmi
# et al. (2008) short-read reference to produce the final
# merged annotation.
#
# Pipeline:
#   1. For each condition, read {cond}_expressed.tsv and
#      derive 5'/3' UTR lengths from the position of the first
#      and last expressed nucleotide.
#   2. Within the NS condition group (NS, NS_PS, NS_RDT) and
#      the S10 condition group (S10, S10_PS, S10_RDT), take
#      the per-boundary maximum across replicates.
#   3. Combine NS and S10 groups into a single DRS UTR table
#      (DRS_UTR_corrected.tsv) by taking the per-boundary
#      maximum across the two groups.
#   4. Merge with Nagalakshmi_UTR.csv: for every gene and
#      every boundary, retain the larger of the two values.
#      This ensures we never shorten a Nagalakshmi boundary.
#
# Inputs:
#   {cond}_expressed.tsv  for each condition (from script 02)
#   Nagalakshmi_UTR.csv   (Gene, five_prime_utr, three_prime_utr)
#
# Outputs (all written to RESULTS_DIR):
#   {cond}_utr.tsv            per-condition UTR table
#   max_NS_utr.tsv            max across NS condition group
#   max_S10_utr.tsv           max across S10 condition group
#   DRS_UTR_corrected.tsv     max across ALL conditions  ← pre-merge DRS annotation
#   final_utr.tsv             merged DRS + Nagalakshmi   ← final annotation
#
# R packages required: dplyr, readr
# ============================================================

library(dplyr)
library(readr)

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these before running
# ════════════════════════════════════════════════════════════════════════

# Directory containing {cond}_expressed.tsv files (output of script 02)
RESULTS_DIR <- "/path/to/data/results"

# Path to the Nagalakshmi et al. reference UTR table
# Required columns: Gene, five_prime_utr, three_prime_utr (NAs allowed)
NAGALAKSHMI_UTR <- "/path/to/Nagalakshmi_UTR.csv"

# Conditions to include (must match the cond names used in scripts 01–02)
NS_CONDITIONS  <- c("NS", "NS_PS", "NS_RDT")
S10_CONDITIONS <- c("S10", "S10_PS", "S10_RDT")
ALL_CONDITIONS <- c(NS_CONDITIONS, S10_CONDITIONS)

# ════════════════════════════════════════════════════════════════════════


# ── Helper: convert expressed positions to UTR lengths ─────────────────
#
# For a gene with n positions (1..n) mapped to a transcript of length n:
#   - Positions 1..1000      correspond to the 5' window
#   - Positions 1001..(n-1000) correspond to the ORF and any UTR
#   - Positions (n-999)..n   correspond to the 3' window
#
# 5'UTR length = max(0, 1000 - first_expressed_position)
# 3'UTR length = max(0, last_expressed_position - (n - 1000))
compute_utrs <- function(db) {
  gene_list <- unique(db$gene)
  utr_rows  <- vector("list", length(gene_list))

  for (j in seq_along(gene_list)) {
    g      <- gene_list[j]
    sub_db <- db[db$gene == g, ]
    n      <- nrow(sub_db)
    borne_inf <- 1000
    borne_sup <- n - 1000

    expressed_pos <- sub_db$position[sub_db$expressed == TRUE]

    if (length(expressed_pos) > 0) {
      five_prime  <- max(0, borne_inf - min(expressed_pos))
      three_prime <- max(0, max(expressed_pos) - borne_sup)
    } else {
      five_prime  <- 0
      three_prime <- 0
    }
    utr_rows[[j]] <- data.frame(gene = g,
                                 five_prime_utr  = five_prime,
                                 three_prime_utr = three_prime,
                                 stringsAsFactors = FALSE)
  }
  bind_rows(utr_rows)
}


# ════════════════════════════════════════════════════════════════════════
# STEP 1 — Compute per-condition UTR tables
# ════════════════════════════════════════════════════════════════════════
message("Step 1: computing per-condition UTR tables...")

cond_utrs <- list()
for (cond in ALL_CONDITIONS) {
  fp <- file.path(RESULTS_DIR, paste0(cond, "_expressed.tsv"))
  if (!file.exists(fp)) {
    warning(sprintf("File not found, skipping condition '%s': %s", cond, fp))
    next
  }
  db       <- read_tsv(fp, col_types = cols()) %>%
    select(gene, position, expressed)
  utr_db   <- compute_utrs(db)
  cond_utrs[[cond]] <- utr_db

  out <- file.path(RESULTS_DIR, paste0(cond, "_utr.tsv"))
  write.table(utr_db, file = out, row.names = FALSE, sep = "\t", quote = FALSE)
  message(sprintf("  Written: %s  (%d genes)", basename(out), nrow(utr_db)))
}


# ── Helper: merge list of UTR tables, keeping per-boundary maximum ─────
merge_max <- function(utr_list, suffix_5, suffix_3) {
  # Rename columns to avoid collisions in join
  named <- mapply(function(df, cond) {
    df %>% rename(!!paste0("five_prime_utr_",  cond) := five_prime_utr,
                  !!paste0("three_prime_utr_", cond) := three_prime_utr)
  }, utr_list, names(utr_list), SIMPLIFY = FALSE)

  merged <- Reduce(function(a, b) full_join(a, b, by = "gene"), named)

  five_cols  <- grep("^five_prime_utr_",  names(merged), value = TRUE)
  three_cols <- grep("^three_prime_utr_", names(merged), value = TRUE)

  merged %>%
    mutate(max_five_prime_utr  = do.call(pmax, c(across(all_of(five_cols)),  list(na.rm = TRUE))),
           max_three_prime_utr = do.call(pmax, c(across(all_of(three_cols)), list(na.rm = TRUE)))) %>%
    select(gene, max_five_prime_utr, max_three_prime_utr)
}


# ════════════════════════════════════════════════════════════════════════
# STEP 2 — Per-group maximum (NS group, S10 group)
# ════════════════════════════════════════════════════════════════════════
message("Step 2: per-group maximum UTR tables...")

ns_conds  <- intersect(NS_CONDITIONS,  names(cond_utrs))
s10_conds <- intersect(S10_CONDITIONS, names(cond_utrs))

max_NS  <- merge_max(cond_utrs[ns_conds],  "five_prime_utr", "three_prime_utr")
max_S10 <- merge_max(cond_utrs[s10_conds], "five_prime_utr", "three_prime_utr")

write.table(max_NS,  file.path(RESULTS_DIR, "max_NS_utr.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)
write.table(max_S10, file.path(RESULTS_DIR, "max_S10_utr.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)
message(sprintf("  max_NS_utr.tsv  : %d genes", nrow(max_NS)))
message(sprintf("  max_S10_utr.tsv : %d genes", nrow(max_S10)))


# ════════════════════════════════════════════════════════════════════════
# STEP 3 — DRS_UTR_corrected.tsv  (max across ALL conditions)
#
# This is the pre-merge DRS annotation used as input to the
# Nagalakshmi merge (Step 4). Column names follow the convention
# expected by plot_figure2.py / plot_drs_vs_nagalakshmi_comparison.R.
# ════════════════════════════════════════════════════════════════════════
message("Step 3: building DRS_UTR_corrected.tsv (max across all conditions)...")

drs_all <- full_join(max_NS, max_S10, by = "gene", suffix = c("_NS", "_S10")) %>%
  mutate(max_five_prime_utr  = pmax(max_five_prime_utr_NS,  max_five_prime_utr_S10,  na.rm = TRUE),
         max_three_prime_utr = pmax(max_three_prime_utr_NS, max_three_prime_utr_S10, na.rm = TRUE)) %>%
  select(gene, max_five_prime_utr, max_three_prime_utr)

write.table(drs_all, file.path(RESULTS_DIR, "DRS_UTR_corrected.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)
message(sprintf("  DRS_UTR_corrected.tsv : %d genes", nrow(drs_all)))


# ════════════════════════════════════════════════════════════════════════
# STEP 4 — final_utr.tsv  (merge DRS with Nagalakshmi)
#
# For every gene and every boundary (5' and 3'), take the
# per-boundary maximum of the DRS annotation and the Nagalakshmi
# reference.  This guarantees we never shorten a Nagalakshmi boundary:
# if the DRS value is larger it wins; if Nagalakshmi is larger it wins.
# Genes present only in Nagalakshmi (not detected in DRS) retain their
# Nagalakshmi boundaries; genes present only in DRS keep their DRS values.
# ════════════════════════════════════════════════════════════════════════
message("Step 4: merging DRS annotation with Nagalakshmi reference...")

nag <- read.csv(NAGALAKSHMI_UTR, stringsAsFactors = FALSE)
colnames(nag)[1:3] <- c("Gene", "nag_five_prime_utr", "nag_three_prime_utr")
nag$nag_five_prime_utr  <- suppressWarnings(as.numeric(nag$nag_five_prime_utr))
nag$nag_three_prime_utr <- suppressWarnings(as.numeric(nag$nag_three_prime_utr))

message(sprintf("  Nagalakshmi: %d genes total, %d with valid 5' UTR, %d with valid 3' UTR",
                nrow(nag),
                sum(!is.na(nag$nag_five_prime_utr)),
                sum(!is.na(nag$nag_three_prime_utr))))

final <- full_join(drs_all, nag, by = c("gene" = "Gene")) %>%
  mutate(
    final_five_prime_utr  = pmax(max_five_prime_utr,  nag_five_prime_utr,  na.rm = TRUE),
    final_three_prime_utr = pmax(max_three_prime_utr, nag_three_prime_utr, na.rm = TRUE)
  ) %>%
  select(Gene = gene, final_five_prime_utr, final_three_prime_utr)

# Sanity check: final boundaries must be >= DRS and >= Nagalakshmi
# (if both are non-NA); any violation indicates a coding error.
joined_check <- full_join(drs_all, nag, by = c("gene" = "Gene")) %>%
  full_join(final, by = c("gene" = "Gene"))

drs_violations <- with(joined_check,
  sum(!is.na(max_five_prime_utr) &
      !is.na(final_five_prime_utr) &
      final_five_prime_utr < max_five_prime_utr, na.rm = TRUE) +
  sum(!is.na(max_three_prime_utr) &
      !is.na(final_three_prime_utr) &
      final_three_prime_utr < max_three_prime_utr, na.rm = TRUE))

nag_violations <- with(joined_check,
  sum(!is.na(nag_five_prime_utr) &
      !is.na(final_five_prime_utr) &
      final_five_prime_utr < nag_five_prime_utr, na.rm = TRUE) +
  sum(!is.na(nag_three_prime_utr) &
      !is.na(final_three_prime_utr) &
      final_three_prime_utr < nag_three_prime_utr, na.rm = TRUE))

if (drs_violations > 0 || nag_violations > 0) {
  stop(sprintf(
    "MERGE VIOLATION: %d final boundaries shorter than DRS; %d shorter than Nagalakshmi",
    drs_violations, nag_violations))
} else {
  message("  Sanity check passed: zero merge violations.")
}

out_final <- file.path(RESULTS_DIR, "final_utr.tsv")
write.table(final, out_final, row.names = FALSE, sep = "\t", quote = FALSE)

message(sprintf("  final_utr.tsv written: %d genes", nrow(final)))
message(sprintf("    Genes with non-zero 5' UTR: %d", sum(final$final_five_prime_utr > 0, na.rm = TRUE)))
message(sprintf("    Genes with non-zero 3' UTR: %d", sum(final$final_three_prime_utr > 0, na.rm = TRUE)))
message("")
message("Output files written to: ", RESULTS_DIR)
message("  DRS_UTR_corrected.tsv  — pre-merge DRS annotation")
message("  final_utr.tsv          — final merged annotation (DRS ∪ Nagalakshmi)")
