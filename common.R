# Data and outputs are kept beside the code folder.
# Set MTN_DATA and MTN_OUTPUT to use different locations.

entry <- grep("^--file=", commandArgs(FALSE), value = TRUE)
entry <- sub("^--file=", "", entry)

if (length(entry) != 1L) {
  stop("Run scripts with Rscript --vanilla.")
}

code_dir <- dirname(normalizePath(entry))

data_dir <- Sys.getenv("MTN_DATA", file.path(code_dir, "..", "mtn-data"))
output_dir <- Sys.getenv("MTN_OUTPUT", file.path(code_dir, "..", "mtn-output"))

data_dir <- normalizePath(data_dir, mustWork = FALSE)
output_dir <- normalizePath(output_dir, mustWork = FALSE)

need_packages <- function(packages) {

  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  missing <- packages[!installed]

  if (length(missing)) {
    stop("Install these R packages first: ", paste(missing, collapse = ", "))
  }

}

input_path <- function(...) {

  path <- file.path(data_dir, ...)

  if (!file.exists(path)) {
    stop("Missing input: ", path, ". See README.md.")
  }

  path

}

derived_path <- input_path

rna_count_file <- function(tissue, required = TRUE) {

  assay <- switch(tissue, bud = "rna_bud", seedling = "rna_seedling")

  paths <- c(
    file.path(data_dir, "rnaseq", assay, "gene_count_matrix.csv"),
    file.path(output_dir, "preprocessing", assay, "gene_count_matrix.csv")
  )

  available <- paths[file.exists(paths)]

  if (length(available)) {
    return(available[1])
  }

  if (required) {
    stop("Missing original ", tissue, " count matrix. Supply ", paths[1],
      " or run rnaseq.sh ", assay, " counts first. See README.md.")
  }

  NULL

}

read_rna_counts <- function(tissue) {

  path <- rna_count_file(tissue)
  x <- readr::read_csv(path, show_col_types = FALSE)

  if (names(x)[1] != "gene_id" || anyNA(x$gene_id) || anyDuplicated(x$gene_id)) {
    stop("Missing or duplicated gene identifiers in ", path)
  }

  if (tissue == "bud") {
    groups <- c("gm", "gs", "mm", "ms", "wm", "ws")
    replicates <- 1:4
  } else {
    groups <- c("5m", "5s", "mm", "ms", "wm", "ws")
    replicates <- 1:3
  }

  samples <- unlist(lapply(groups, function(group) paste0(group, replicates)))

  if (!setequal(names(x)[-1], samples)) {
    stop("Unexpected sample columns in ", path,
      ". Use the original short sample names listed in README.md.")
  }

  if (!all(vapply(x[, samples], is.numeric, logical(1)))) {
    stop("Non-numeric gene counts in ", path)
  }

  counts <- as.matrix(x[, samples])
  rownames(counts) <- x$gene_id

  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
    stop("Gene counts must be finite, non-negative integers in ", path)
  }

  storage.mode(counts) <- "integer"
  counts

}

start_analysis <- function(name, packages = character()) {

  need_packages(packages)

  out <- file.path(output_dir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  setwd(out)

  options(device = function(...) grDevices::pdf(file = NULL, ...))
  set.seed(1)

  message("Writing ", out)

}
