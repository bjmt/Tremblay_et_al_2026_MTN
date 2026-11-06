# Read the published supplementary workbook and Zenodo measurements.

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))

need_packages(c("readxl", "readr"))

workbook <- input_path("supplementary_datasets.xlsx")

# Metabolite measurements: Dataset S1.

met <- readxl::read_xlsx(workbook, sheet = "Dataset S1", skip = 1)
met <- met[!is.na(met$Genotype), ]

names(met)[1:6] <- c("genotype", "treatment", "tissue", "SUL", "PHO", "NIT")

met$genotype <- tolower(met$genotype)
met$genotype[grepl("mtn1-1mtn2-1", met$genotype)] <- "mtn1121"
met$genotype[met$genotype == "mtn1-1mtn2-5"] <- "mtn1125"
met$treatment <- tolower(met$treatment)
met$tissue <- tolower(met$tissue)

met_folder <- file.path(data_dir, "metabolites")
dir.create(met_folder, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(met, file.path(met_folder, "metabolite_measurements.csv"))

# RNA-seq contrasts: Datasets S2 and S3.
# HOMER names a contrast "A vs. B" with log2 fold change B/A.

pairs <- list(
  bud = c("g3_ms vs. g3_spd", "g3_ms vs. mtn_ms", "g3_ms vs. wt_ms",
    "g3_spd vs. mtn_spd", "g3_spd vs. wt_spd", "mtn_ms vs. mtn_spd",
    "mtn_ms vs. wt_ms", "mtn_spd vs. wt_spd", "wt_ms vs. wt_spd"),
  seedling = c("5m vs. 5s", "5m vs. mm", "5m vs. wm",
    "5s vs. ms", "5s vs. ws", "mm vs. ms", "mm vs. wm",
    "ms vs. ws", "wm vs. ws")
)

rna_folder <- file.path(data_dir, "rnaseq")
dir.create(rna_folder, recursive = TRUE, showWarnings = FALSE)

for (tissue in c("bud", "seedling")) {

  sheet <- if (tissue == "bud") "Dataset S3" else "Dataset S2"
  x <- readxl::read_xlsx(workbook, sheet = sheet, skip = 1)

  if (names(x)[1] != "Gene" || !all(grepl("Log2 Fold Change", names(x)[seq(2, 26, 3)]))) {
    stop("Unexpected columns in ", sheet)
  }

  columns <- "Gene"

  for (comparison in pairs[[tissue]]) {
    columns <- c(columns, paste(comparison,
      c("Log2 Fold Change", "p-value", "adj. p-value")))
  }

  names(x)[1:28] <- columns

  if (anyNA(x$Gene) || anyDuplicated(x$Gene)) {
    stop("Missing or duplicated gene identifiers in ", sheet)
  }

  readr::write_tsv(x, file.path(rna_folder, paste0(tissue, ".tsv")))

}

# Stack the GO term tables in Dataset S4. This dataset contains term P-values,
# rather than the gene fractions needed to reproduce the published dot plot.

go_sheet <- readxl::read_xlsx(workbook, sheet = "Dataset S4",
  col_names = FALSE, col_types = "text", .name_repair = "minimal")

go <- list()

for (column in seq_len(ncol(go_sheet) - 3L)) {

  if (!identical(go_sheet[[column]][3], "p_value")) {
    next
  }

  table <- data.frame(
    Comparison = go_sheet[[column]][2],
    Pvalue = suppressWarnings(as.numeric(go_sheet[[column]][-(1:3)])),
    TermID = go_sheet[[column + 1L]][-(1:3)],
    Source = go_sheet[[column + 2L]][-(1:3)],
    Term = go_sheet[[column + 3L]][-(1:3)]
  )

  go[[length(go) + 1L]] <- table[!is.na(table$Pvalue), ]

}

readr::write_tsv(do.call(rbind, go), file.path(rna_folder, "go_terms.tsv"))

# Optional count-based correlations from the two GEO count matrices.

bud_counts <- file.path(data_dir, "flowerbuds_counts.csv")
seedling_counts <- file.path(data_dir, "seedlings_counts.csv")

if (file.exists(bud_counts) && file.exists(seedling_counts)) {

  bud <- readr::read_csv(bud_counts, show_col_types = FALSE)
  seedling <- readr::read_csv(seedling_counts, show_col_types = FALSE)

  bud_columns <- character()
  seedling_columns <- character()

  for (genotype in c("wt", "mtn11mtn21", "g3")) {
    for (treatment in c("ms", "spd")) {
      bud_columns <- c(bud_columns, paste0("bud_", genotype, "_", treatment, "_", 1:4))
    }
  }

  for (genotype in c("mtn11mtn21", "mtn11mtn25", "wt")) {
    for (treatment in c("ms", "spd")) {
      seedling_columns <- c(seedling_columns,
        paste0("sdl_", genotype, "_", treatment, "_", 1:3))
    }
  }

  genes <- bud$gene_id[bud$gene_id %in% seedling$gene_id]

  counts <- cbind(
    data.frame(gene_id = genes),
    as.data.frame(bud[match(genes, bud$gene_id), bud_columns]),
    as.data.frame(seedling[match(genes, seedling$gene_id), seedling_columns])
  )

  readr::write_csv(counts, file.path(rna_folder, "counts.csv"))

}

# Phenotype tables from the Zenodo downloads.

zenodo <- file.path(data_dir, "zenodo")
plant_file <- file.path(zenodo, "Plant Developmental monitoring.xlsx")
root_file <- file.path(zenodo, "Root lengths (treatments).csv")

if (file.exists(plant_file) && file.exists(root_file)) {

  folder <- file.path(data_dir, "phenotypes")
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)

  roots <- read.csv(root_file, header = FALSE, stringsAsFactors = FALSE)

  if (trimws(roots[14, 8]) != "WT" || trimws(roots[30, 8]) != "aabb") {
    stop("Unexpected layout in Root lengths (treatments).csv")
  }

  wt <- roots[16:25, 8:12]
  mtn <- roots[32:41, 8:12]

  roots <- rbind(data.frame(Genotype = "WT", wt), data.frame(Genotype = "aabb", mtn))
  names(roots) <- c("Genotype", "MS", "NA", "Spd", "Spm", "Tspm")

  for (column in 2:6) {
    roots[[column]] <- as.numeric(roots[[column]])
  }

  readr::write_tsv(roots, file.path(folder, "root_lengths.tsv"))

  seeds <- readxl::read_xlsx(plant_file, sheet = "NA+polyamines",
    col_types = "text", .name_repair = "minimal")
  seeds <- seeds[, c(1:3, 5:7)]
  seeds <- seeds[!is.na(seeds[[1]]), ]

  names(seeds) <- c("WT", "WT_weight", "WT_count", "aabb", "aabb_weight", "aabb_count")

  for (column in c(2, 3, 5, 6)) {
    seeds[[column]] <- as.numeric(seeds[[column]])
  }

  readr::write_tsv(seeds, file.path(folder, "seed_yields.tsv"))

  development <- list()

  for (sheet in c("WT MS", "WT Spd", "dm MS", "dm Spd", "G3 MS", "G3 Spd")) {

    x <- readxl::read_xlsx(plant_file, sheet = sheet,
      col_types = "text", .name_repair = "minimal")
    names(x) <- trimws(names(x))

    # The seed-count and bud-count columns occur in a different order in dm MS.
    weight_column <- grep("^weight of seeds", names(x))[1]
    bud_column <- which(names(x) %in% c("Buds cluster", "cluster of buds"))[1]

    columns <- c(1, match(c("Planted", "Bolted", "Days to Bolt", "Leaves at Bolting",
      "Size of Plant", "Number of Siliques", "Size of Silique"), names(x)),
      weight_column, match("Number of Seeds", names(x)), bud_column)

    if (anyNA(columns)) {
      stop("Missing phenotype columns in ", sheet)
    }

    x <- x[, columns]

    names(x) <- c("Number", "Planted", "Bolted", "Days to Bolt", "Leaves at Bolting",
      "Size of Plant", "Number of Siliques", "Size of Silique",
      "Seed weight", "Number of Seeds", "Buds cluster")

    x$Number <- suppressWarnings(as.numeric(x$Number))
    x <- x[!is.na(x$Number), ]

    for (column in 4:11) {
      x[[column]] <- suppressWarnings(as.numeric(x[[column]]))
    }

    genotype <- strsplit(sheet, " ", fixed = TRUE)[[1]][1]
    treatment <- strsplit(sheet, " ", fixed = TRUE)[[1]][2]

    if (genotype == "dm") {
      genotype <- "mtn1-1mtn2-1"
    }

    x <- cbind(data.frame(Genotype = genotype, Treatment = treatment), as.data.frame(x))
    development[[sheet]] <- x

  }

  development <- do.call(rbind, development)
  readr::write_tsv(development, file.path(folder, "development.tsv"))

}

message("Prepared tables in ", data_dir)
