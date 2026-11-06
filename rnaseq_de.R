# Recalculate differential expression from the deposited gene-count matrices.
# Each contrast is fitted separately with DESeq2, as in HOMER's pairwise workflow.
# Usage: Rscript --vanilla rnaseq_de.R [bud|seedling]

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))

start_analysis("rnaseq_de", c("readr", "DESeq2"))

args <- commandArgs(TRUE)
tissues <- if (length(args)) args else c("bud", "seedling")

if (!all(tissues %in% c("bud", "seedling"))) {
  stop("Choose bud or seedling, or omit the argument to run both")
}

pair_indices <- list(c(1, 2), c(1, 3), c(1, 5), c(2, 4), c(2, 6),
  c(3, 4), c(3, 5), c(4, 6), c(5, 6))

rrna <- c("AT2G01010", "AT2G01020", "AT3G41768", "AT3G41979", "ATCG00920",
  "ATCG00950", "ATCG00960", "ATCG00970", "ATCG01160", "ATCG01170",
  "ATCG01180", "ATCG01210")

agreement <- list()

for (tissue in tissues) {

  prefix <- if (tissue == "bud") "flowerbuds" else "seedlings"
  x <- readr::read_csv(input_path(paste0(prefix, "_counts.csv")), show_col_types = FALSE)

  if (anyNA(x$gene_id) || anyDuplicated(x$gene_id)) {
    stop("Missing or duplicated gene identifiers in ", prefix, " counts")
  }

  counts <- as.matrix(x[, -1])
  rownames(counts) <- x$gene_id

  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
    stop("DESeq2 requires non-negative integer gene counts")
  }

  storage.mode(counts) <- "integer"
  sample_parts <- strsplit(colnames(counts), "_", fixed = TRUE)
  genotype <- vapply(sample_parts, function(x) x[2], "")
  treatment <- vapply(sample_parts, function(x) x[3], "")

  if (tissue == "bud") {

    genotype <- c(g3 = "g3", mtn11mtn21 = "mtn", wt = "wt")[genotype]
    groups <- paste(genotype, treatment, sep = "_")
    group_order <- c("g3_ms", "g3_spd", "mtn_ms", "mtn_spd", "wt_ms", "wt_spd")
    threshold <- 16
    filter_columns <- seq_len(ncol(counts))

  } else {

    genotype <- c(mtn11mtn25 = "5", mtn11mtn21 = "m", wt = "w")[genotype]
    treatment <- c(ms = "m", spd = "s")[treatment]
    groups <- paste0(genotype, treatment)
    group_order <- c("5m", "5s", "mm", "ms", "wm", "ws")
    threshold <- 14
    filter_columns <- 1:16

  }

  if (anyNA(genotype) || !all(groups %in% group_order)) {
    stop("Unexpected sample labels in ", prefix, " counts")
  }

  # Apply the same whole-experiment low-count and rRNA exclusions as rnaseq.sh.
  keep <- rowSums(counts[, filter_columns, drop = FALSE]) > threshold
  counts <- counts[keep & !rownames(counts) %in% rrna, , drop = FALSE]
  result_table <- data.frame(Gene = rownames(counts))

  for (indices in pair_indices) {

    pair <- group_order[indices]
    comparison <- paste(pair, collapse = " vs. ")
    selected <- which(groups %in% pair)

    message(tissue, ": ", comparison, " (", nrow(counts), " genes)")

    samples <- data.frame(
      row.names = colnames(counts)[selected],
      Condition = factor(groups[selected], levels = pair)
    )

    dds <- DESeq2::DESeqDataSetFromMatrix(counts[, selected, drop = FALSE],
      samples, design = ~Condition)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    result <- DESeq2::results(dds, contrast = c("Condition", pair[2], pair[1]),
      alpha = .1)

    result_table[[paste(comparison, "Log2 Fold Change")]] <- result$log2FoldChange
    result_table[[paste(comparison, "p-value")]] <- result$pvalue
    result_table[[paste(comparison, "adj. p-value")]] <- result$padj

  }

  # Fit before excluding novel StringTie identifiers: they contribute to fitting
  # and multiplicity adjustment, but are not retained in the gene figure tables.
  result_table <- result_table[!grepl("^STR", result_table$Gene), ]
  readr::write_tsv(result_table, paste0(tissue, ".tsv"))

  # TPMs are only needed to pass these tables to the plotting script.
  tpm_file <- file.path(data_dir, paste0(prefix, "_tpm.csv"))

  if (file.exists(tpm_file)) {

    tpm <- readr::read_csv(tpm_file, show_col_types = FALSE)
    columns <- character()

    if (tissue == "bud") {

      for (genotype in c("g3", "mtn11mtn21", "wt")) {
        for (treatment in c("ms", "spd")) {
          columns <- c(columns, paste0("bud_", genotype, "_", treatment, "_", 1:4))
        }
      }

    } else {

      for (treatment in c("ms", "spd")) {
        columns <- c(columns, paste0("sdl_mtn11mtn21_", treatment, "_", 1:3))
      }

      for (genotype in c("mtn11mtn25", "wt")) {
        for (replicate in 1:3) {
          columns <- c(columns,
            paste0("sdl_", genotype, "_", c("ms", "spd"), "_", replicate))
        }
      }

    }

    if (!all(columns %in% names(tpm))) {
      stop("Unexpected TPM sample names for ", tissue)
    }

    missing_tpm <- sum(!result_table$Gene %in% tpm$gene_id)

    if (missing_tpm) {
      message(missing_tpm, " genes lack deposited TPMs; their TPM values will remain NA")
    }

    for (column in columns) {
      result_table[[paste(column, "TPM")]] <- tpm[[column]][match(result_table$Gene, tpm$gene_id)]
    }

  } else {

    message("TPM matrix absent: DE results will be written without correlation-plot inputs")

  }

  readr::write_tsv(result_table, paste0(tissue, ".tsv"))

  # Optional comparison against the supplied statistics, keeping NA tests as NA.
  reference_file <- file.path(data_dir, "rnaseq", paste0(tissue, ".tsv"))

  if (file.exists(reference_file)) {

    reference <- readr::read_tsv(reference_file, show_col_types = FALSE)
    index <- match(result_table$Gene, reference$Gene)

    for (indices in pair_indices) {

      comparison <- paste(group_order[indices], collapse = " vs. ")
      fc_name <- paste(comparison, "Log2 Fold Change")
      fdr_name <- paste(comparison, "adj. p-value")

      fc <- result_table[[fc_name]]
      expected_fc <- reference[[fc_name]][index]
      fdr <- result_table[[fdr_name]]
      expected_fdr <- reference[[fdr_name]][index]
      finite <- is.finite(fc) & is.finite(expected_fc)

      significant <- is.finite(fdr) & fdr < .05 & is.finite(fc) & abs(fc) >= 1
      expected_significant <- is.finite(expected_fdr) & expected_fdr < .05 &
        is.finite(expected_fc) & abs(expected_fc) >= 1

      agreement[[length(agreement) + 1L]] <- data.frame(
        Tissue = tissue,
        Comparison = comparison,
        MatchedGenes = sum(!is.na(index)),
        FiniteFoldChanges = sum(finite),
        FoldChangeCorrelation = cor(fc[finite], expected_fc[finite]),
        MaximumFoldChangeDifference = max(abs(fc[finite] - expected_fc[finite])),
        DifferentSignificanceCalls = sum(significant != expected_significant, na.rm = TRUE)
      )

    }

  }

}

if (length(agreement)) {
  readr::write_tsv(do.call(rbind, agreement), "comparison_with_supplied_results.tsv")
}

writeLines(capture.output(sessionInfo()), "DESeq2_sessionInfo.txt")
