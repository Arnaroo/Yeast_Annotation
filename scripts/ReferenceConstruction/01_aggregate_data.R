# ============================================================
# 01_aggregate_data.R
#
# Aggregate per-replicate DRS read-count profiles into one
# combined profile per condition.
#
# Input:  per-replicate TSV files named {condition}_{replicate}.txt
#         Three columns (no header): gene  position  transcript_np
#         All replicates for a condition must be in DATA_DIR.
#
# Output: {condition}_combined.tsv in DATA_DIR
#         Two columns with header: gene  position  transcript_np
#         (transcript_np is the sum across all replicates)
#
# R packages required: dplyr, readr
# ============================================================

library(dplyr)
library(readr)

# ════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these before running
# ════════════════════════════════════════════════════════════════════════

# Directory containing the per-replicate input files
DATA_DIR <- "/path/to/data/gene_expressed"

# Conditions to aggregate.
# For each condition, all files matching {condition}*.txt in DATA_DIR
# are summed. Adjust the list to match your experimental design.
CONDITIONS <- c("NS_PS", "NS_RDT", "S10_PS", "S10_RDT")

# ════════════════════════════════════════════════════════════════════════


aggregate_replicates <- function(cond, data_dir) {
  # List all replicate files for this condition
  file_list <- list.files(
    path       = data_dir,
    pattern    = paste0("^", cond, ".*\\.txt$"),
    full.names = TRUE
  )

  if (length(file_list) == 0) {
    warning(sprintf("No files found for condition '%s' in %s", cond, data_dir))
    return(invisible(NULL))
  }
  message(sprintf("[%s] Found %d replicate file(s):", cond, length(file_list)))
  message(paste(" ", basename(file_list), collapse = "\n"))

  # Read each replicate file and combine
  data_list <- lapply(file_list, function(file) {
    df <- read.table(file, header = FALSE)
    colnames(df) <- c("gene", "position", "transcript_np")
    df
  })

  combined <- bind_rows(data_list)

  # Sum transcript_np across all replicates for each gene × position
  final <- combined %>%
    group_by(gene, position) %>%
    summarise(transcript_np = sum(transcript_np), .groups = "drop")

  output_file <- file.path(data_dir, paste0(cond, "_combined.tsv"))
  write_tsv(final, output_file)
  message(sprintf("[%s] Written: %s  (%d rows)", cond, output_file, nrow(final)))
}


# Run for all conditions
for (cond in CONDITIONS) {
  aggregate_replicates(cond, DATA_DIR)
}

message("Done. Combined profiles written to: ", DATA_DIR)
