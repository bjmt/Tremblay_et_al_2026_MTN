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

start_analysis <- function(name, packages = character()) {

  need_packages(packages)

  out <- file.path(output_dir, name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  setwd(out)

  options(device = function(...) grDevices::pdf(file = NULL, ...))
  set.seed(1)

  message("Writing ", out)

}
