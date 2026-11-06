# Annotate gene-level differential expression and add StringTie TPM values.
# Usage: Rscript annotate_rnaseq.R assay work_directory annotation.gff3

args <- commandArgs(TRUE)

if (length(args) != 3L) {
  stop("Usage: annotate_rnaseq.R assay work_directory annotation.gff3")
}

assay <- args[1]

if (!assay %in% c("rna_bud", "rna_seedling")) {
  stop("Unknown assay")
}

work <- normalizePath(args[2], mustWork = TRUE)
annotation <- normalizePath(args[3], mustWork = TRUE)

for (package in c("readr", "rtracklayer")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing package: ", package)
  }
}

setwd(work)

# Keep reference genes and the genotype/treatment contrasts used in the analyses.

x <- readr::read_tsv("gene_count_diff.tsv", show_col_types = FALSE)
x <- x[!grepl("^STR", x[[1]]), ]

if (assay == "rna_bud") {
  groups <- c("g3_ms", "g3_spd", "mtn_ms", "mtn_spd", "wt_ms", "wt_spd")
} else {
  groups <- c("5m", "5s", "mm", "ms", "wm", "ws")
}

pairs <- list(c(1, 2), c(1, 3), c(1, 5), c(2, 4), c(2, 6),
  c(3, 4), c(3, 5), c(4, 6), c(5, 6))

columns <- character()

for (pair in pairs) {
  comparison <- paste(groups[pair], collapse = " vs. ")
  columns <- c(columns, paste(comparison,
    c("Log2 Fold Change", "p-value", "adj. p-value")))
}

if (!all(columns %in% names(x))) {
  stop("Unexpected HOMER output column names")
}

x <- as.data.frame(x[, c(names(x)[1], columns)])
names(x)[1] <- "Gene"

# Add TPMs from each sample, matching by gene identifier.

files <- sort(list.files("stringtie", "[.]gene_abund[.]txt$", full.names = TRUE))

for (file in files) {

  sample <- sub("[.]gene_abund[.]txt$", "", basename(file))
  abundance <- readr::read_tsv(file, show_col_types = FALSE)

  x[[paste(sample, "TPM")]] <- abundance$TPM[match(x$Gene, abundance$`Gene ID`)]

}

# Keep missing statistical tests as NA.

annotation <- rtracklayer::import(annotation)

fields <- c(
  symbol = "Gene symbol",
  full_name = "Gene name",
  locus_type = "Gene type",
  computational_description = "Gene description"
)

for (field in names(fields)) {

  values <- as.character(S4Vectors::mcols(annotation)[[field]])
  names(values) <- annotation$Name

  x[[fields[[field]]]] <- values[x$Gene]

}

readr::write_tsv(x, "gene_count_diff.annotated.tsv")
