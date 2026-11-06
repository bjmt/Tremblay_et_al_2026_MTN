# Plot full GO enrichment exports, or recalculate enrichment with g:Profiler.
# Usage: Rscript --vanilla rnaseq_go.R [--recalculate] [contrast_folder]

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("rnaseq_go", c("ggplot2", "readr"))

library(ggplot2)

args <- commandArgs(TRUE)
recalculate <- "--recalculate" %in% args
args <- args[args != "--recalculate"]
contrast_folder <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  file.path(data_dir, "rnaseq")

comparisons <- data.frame(
  File = c("WTms_WTspd", "mtnms_mtnspd", "G3ms_G3spd", "mtnms_WTms",
    "G3ms_WTms", "G3ms_mtnms", "mtnspd_WTspd", "G3spd_WTspd", "G3spd_mtnspd"),
  Contrast = c("wt_ms vs. wt_spd", "mtn_ms vs. mtn_spd", "g3_ms vs. g3_spd",
    "mtn_ms vs. wt_ms", "g3_ms vs. wt_ms", "g3_ms vs. mtn_ms",
    "mtn_spd vs. wt_spd", "g3_spd vs. wt_spd", "g3_spd vs. mtn_spd"),
  Label = c("WT +Spd vs -Spd", "mtn +Spd vs -Spd", "G3 +Spd vs -Spd",
    "mtn vs WT -Spd", "G3 vs WT -Spd", "G3 vs mtn -Spd",
    "mtn vs WT +Spd", "G3 vs WT +Spd", "G3 vs mtn +Spd"),
  Reverse = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
)

if (recalculate) {

  need_packages("gprofiler2")
  bud <- readr::read_tsv(file.path(contrast_folder, "bud.tsv"), show_col_types = FALSE)
  folder <- file.path(getwd(), "enrichment")
  dir.create(folder, showWarnings = FALSE)

} else {

  folder <- input_path("rnaseq", "enrichment")

}

tables <- list()

for (i in seq_len(nrow(comparisons))) {
  for (direction in c("up", "down")) {

    file <- file.path(folder, paste0("Bud_", comparisons$File[i], "_", direction, ".txt"))

    if (recalculate) {

      fc <- bud[[paste(comparisons$Contrast[i], "Log2 Fold Change")]]
      fdr <- bud[[paste(comparisons$Contrast[i], "adj. p-value")]]
      selected <- if (direction == "up") fc >= 1 else fc <= -1
      genes <- bud$Gene[which(is.finite(fdr) & fdr < .05 & selected)]

      if (!length(genes)) {
        next
      }

      result <- gprofiler2::gost(genes, organism = "athaliana")

      if (is.null(result$result)) {
        next
      }

      x <- result$result[, c("term_id", "term_name", "p_value", "source",
        "query_size", "intersection_size")]
      readr::write_tsv(x, file)

    } else {

      if (!file.exists(file)) {
        next
      }

      x <- readr::read_tsv(file, show_col_types = FALSE)

    }

    x$Comparison <- comparisons$Label[i]
    up <- xor(direction == "up", comparisons$Reverse[i])
    x$Direction <- if (up) "Up-regulated" else "Down-regulated"
    x$Fraction <- x$intersection_size / x$query_size
    tables[[length(tables) + 1L]] <- x[, c("Comparison", "Direction", "term_id",
      "term_name", "p_value", "Fraction")]

  }
}

if (!length(tables)) {
  stop("No GO enrichment results found")
}

go <- do.call(rbind, tables)
readr::write_tsv(go, "go_terms.tsv")

terms <- c("pollination", "extracellular region", "regulation of fertilization",
  "response to red or far red light", "cell wall organization",
  "pectinesterase activity", "response to abiotic stimulus",
  "response to light stimulus", "cellular response to hypoxia", "photosystem I")

go <- go[go$term_name %in% terms, ]
go$term_name <- factor(go$term_name, levels = rev(sort(terms)))
go$Comparison <- factor(go$Comparison, levels = comparisons$Label)
go$Direction <- factor(go$Direction, levels = c("Up-regulated", "Down-regulated"))

plot <- ggplot(go, aes(Comparison, term_name, size = -log10(p_value), fill = Fraction)) +
  geom_point(shape = 21) +
  facet_wrap(~Direction) +
  scale_fill_gradient2(low = "pink", mid = "red", high = "black",
    midpoint = mean(range(go$Fraction)), name = "Fraction of genes") +
  scale_size(name = expression(-log[10](P))) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("go_terms.pdf", plot, width = 11, height = 5)
