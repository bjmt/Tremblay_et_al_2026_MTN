# Call DMRs from the Bismark cytosine reports using the original DMRcaller workflow.
# Usage: Rscript --vanilla dmr_calling.R [comparison]

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("dmr_calling", c("DMRcaller", "betareg", "readr"))

library(DMRcaller)
library(betareg)

comparisons <- list(
  wtms_vs_mtnms = c("wm", "mm"),
  wtms_vs_g3ms = c("wm", "gm"),
  mtnms_vs_g3ms = c("mm", "gm"),
  wtspd_vs_mtnspd = c("ws", "ms"),
  wtspd_vs_g3spd = c("ws", "gs"),
  mtnspd_vs_g3spd = c("ms", "gs"),
  wtms_vs_wtspd = c("wm", "ws"),
  mtnms_vs_mtnspd = c("mm", "ms"),
  g3ms_vs_g3spd = c("gm", "gs")
)

args <- commandArgs(TRUE)

if (length(args)) {

  if (!all(args %in% names(comparisons))) {
    stop("Unknown comparison. Choose: ", paste(names(comparisons), collapse = ", "))
  }

  comparisons <- comparisons[args]

}

folder <- file.path(output_dir, "preprocessing", "bs_bud", "methcalls")
cores <- as.integer(Sys.getenv("THREADS", "6"))

for (comparison in names(comparisons)) {

  groups <- comparisons[[comparison]]
  samples <- c(paste0(groups[1], 1:3), paste0(groups[2], 1:3))
  reports <- character()

  for (sample in samples) {

    files <- list.files(file.path(folder, sample),
      "CX_report[.]txt([.]gz)?$", full.names = TRUE)

    if (length(files) != 1L) {
      stop("Expected one Bismark CX report for ", sample, "; run bisulfite.sh extract first")
    }

    reports <- c(reports, files)

  }

  message("Reading methylation calls for ", comparison)
  replicates <- lapply(reports, DMRcaller::readBismark)
  methylation <- replicates[[1]]

  for (i in seq_along(replicates)[-1]) {
    methylation <- DMRcaller::joinReplicates(methylation, replicates[[i]])
  }

  rm(replicates)
  tables <- list()
  failed <- FALSE

  for (context in c("CG", "CHG", "CHH")) {

    message("Calling ", comparison, " ", context)

    # The recovered analysis uses DMRcaller defaults with bins and six replicates.
    regions <- tryCatch(
      DMRcaller::computeDMRsReplicates(methylation,
        context = context, condition = c(rep("a", 3), rep("b", 3)),
        method = "bins", cores = cores),
      error = function(error) {
        warning(comparison, " ", context, ": ", conditionMessage(error))
        NULL
      }
    )

    if (is.null(regions)) {
      failed <- TRUE
      next
    }

    if (length(regions)) {
      tables[[context]] <- as.data.frame(regions)
    }

  }

  if (length(tables)) {

    result <- do.call(rbind, tables)

  } else {

    if (failed) {
      warning("No table written for ", comparison, ": DMR calling failed")
      rm(methylation)
      next
    }

    result <- data.frame(seqnames = character(), start = integer(), end = integer(),
      context = character(), regionType = character())

  }

  readr::write_tsv(result, paste0(comparison, ".tsv"))
  rm(methylation)

}
