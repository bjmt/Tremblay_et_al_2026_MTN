# Analyse the DMR intervals supplied in Supplementary Dataset S5.
# Intersections use identical chromosome/start/end coordinates, as in the
# supplementary UpSet plots. They do not merge partially overlapping intervals.

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))

start_analysis("dmrs", c("readr", "UpSetR"))

files <- list.files(input_path("dmrs"), "[.]tsv$", full.names = TRUE)
tables <- lapply(files, readr::read_tsv, show_col_types = FALSE)
names(tables) <- sub("[.]tsv$", "", basename(files))

counts <- list()

for (comparison in names(tables)) {
  for (context in c("CG", "CHG", "CHH")) {
    for (direction in c("gain", "loss")) {

      x <- tables[[comparison]]
      x <- x[x$context == context & x$regionType == direction, ]

      counts[[length(counts) + 1L]] <- data.frame(
        Comparison = comparison,
        Context = context,
        Direction = direction,
        Regions = nrow(unique(x[, c("seqnames", "start", "end")]))
      )

    }
  }
}

readr::write_tsv(do.call(rbind, counts), "dmr_counts.tsv")

comparisons <- list(
  genotype_ms = c(
    "wtms_vs_mtnms" = "mtn vs WT",
    "wtms_vs_g3ms" = "G3 vs WT",
    "mtnms_vs_g3ms" = "G3 vs mtn"
  ),
  genotype_spd = c(
    "wtspd_vs_mtnspd" = "mtn vs WT",
    "wtspd_vs_g3spd" = "G3 vs WT",
    "mtnspd_vs_g3spd" = "G3 vs mtn"
  ),
  treatment = c(
    "mtnms_vs_mtnspd" = "mtn +Spd vs -Spd",
    "g3ms_vs_g3spd" = "G3 +Spd vs -Spd"
  )
)

intersections <- list()

for (comparison_group in names(comparisons)) {
  for (context in c("CG", "CHG", "CHH")) {

    sets <- list()

    for (comparison in names(comparisons[[comparison_group]])) {
      for (direction in c("gain", "loss")) {

        x <- tables[[comparison]]
        x <- x[x$context == context & x$regionType == direction, ]
        label <- paste(comparisons[[comparison_group]][comparison],
          if (direction == "gain") "hyper" else "hypo")

        sets[[label]] <- unique(with(x, paste(seqnames, start, end, sep = ":")))

      }
    }

    for (pair in combn(names(sets), 2, simplify = FALSE)) {

      intersections[[length(intersections) + 1L]] <- data.frame(
        ComparisonGroup = comparison_group,
        Context = context,
        Set1 = pair[1],
        Set2 = pair[2],
        SharedIntervals = length(intersect(sets[[pair[1]]], sets[[pair[2]]]))
      )

    }

    sets <- sets[lengths(sets) > 0]

    if (!length(sets)) {
      next
    }

    if (length(sets) == 1L) {

      pdf(paste0(comparison_group, "_", context, ".pdf"), width = 7, height = 4)
      barplot(lengths(sets), ylab = "DMRs")
      dev.off()

    } else {

      plot <- UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets),
        order.by = "freq", keep.order = TRUE)

      pdf(paste0(comparison_group, "_", context, ".pdf"), width = 8, height = 4)
      print(plot)
      dev.off()

    }

  }
}

readr::write_tsv(do.call(rbind, intersections), "dmr_pair_overlaps.tsv")

# Optional gene lists connecting WT/mtn RNA-seq changes to CG and CHG DMRs.

annotation_file <- file.path(data_dir, "reference/Araport11_GFF3_genes_transposons.gff3")

if (file.exists(annotation_file)) {

  need_packages(c("rtracklayer", "GenomicRanges"))

  annotation <- rtracklayer::import(annotation_file)
  genes <- annotation[annotation$type == "gene"]
  GenomeInfoDb::seqlevels(genes) <- sub("^Chr", "", GenomeInfoDb::seqlevels(genes))

  rna <- readr::read_tsv(input_path("rnaseq", "bud.tsv"), show_col_types = FALSE)
  fc <- rna[["mtn_ms vs. wt_ms Log2 Fold Change"]]
  fdr <- rna[["mtn_ms vs. wt_ms adj. p-value"]]

  # The stored coefficient is WT/mtn, so negative values mean higher in mtn.
  selected <- list(
    mtn_up = rna$Gene[which(is.finite(fdr) & fdr < .05 & fc <= -1)],
    mtn_down = rna$Gene[which(is.finite(fdr) & fdr < .05 & fc >= 1)]
  )

  for (context in c("CG", "CHG")) {

    x <- tables[["wtms_vs_mtnms"]]
    regions <- GenomicRanges::makeGRangesFromDataFrame(x[x$context == context, ])
    overlap <- IRanges::overlapsAny(genes, regions, ignore.strand = TRUE)
    overlapping_genes <- as.character(genes$ID[overlap])

    for (direction in names(selected)) {
      readr::write_lines(intersect(selected[[direction]], overlapping_genes),
        paste0(direction, "_", context, "_dmr_genes.txt"))
    }

  }

} else {

  message("Skipping DMR/gene overlaps: supply the Araport11 GFF3 under reference/.")

}
