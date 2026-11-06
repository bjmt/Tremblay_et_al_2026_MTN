# PCA panels for internode lengths and metabolite measurements.
# Usage: Rscript --vanilla pca.R [--recalculate]

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("pca", c("ggplot2", "readr"))

library(ggplot2)

recalculate <- "--recalculate" %in% commandArgs(TRUE)
colours <- c("wt" = "grey60", "mtn1125" = "orange", "mtn1121" = "skyblue",
  "g2" = "seagreen3", "g3" = "yellow")
labels <- c("wt" = "WT", "mtn1125" = "mtn1-1mtn2-5",
  "mtn1121" = "mtn1-1mtn2-1", "g2" = "G2", "g3" = "G3")

plot_scores <- function(x, xlab = "PC1", ylab = "PC2") {

  x$genotype <- factor(x$genotype, levels = names(colours))
  x$treatment <- factor(tolower(x$treatment), levels = c("ms", "spd"))

  ggplot(x, aes(PC1, PC2, fill = genotype, shape = treatment)) +
    geom_point(size = 2.5) +
    scale_fill_manual(values = colours, labels = labels, name = NULL) +
    scale_shape_manual(values = c(21, 24), labels = c("-Spd", "+Spd"), name = NULL) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = xlab, y = ylab) +
    theme_bw()

}

internodes <- readr::read_csv(input_path("phenotypes", "internode_pca.csv"),
  show_col_types = FALSE)
values <- as.matrix(internodes[, as.character(1:10)])

if (anyNA(values)) {
  stop("Missing internode lengths in the PCA matrix")
}

model <- prcomp(values, center = TRUE, scale. = TRUE)
variance <- 100 * model$sdev^2 / sum(model$sdev^2)
parts <- strsplit(internodes$Sample, "_", fixed = TRUE)
scores <- data.frame(
  genotype = vapply(parts, function(x) x[1], ""),
  treatment = vapply(parts, function(x) x[2], ""),
  model$x[, 1:2]
)
scores$genotype <- c("WT" = "wt", "dm" = "mtn1121", "G3" = "g3")[scores$genotype]

plot <- plot_scores(scores, sprintf("PC1 (%.0f%%)", variance[1]),
  sprintf("PC2 (%.0f%%)", variance[2]))
ggsave("internodes.pdf", plot, width = 5, height = 4)
readr::write_tsv(scores, "internode_scores.tsv")

for (tissue in c("seedling", "bud")) {

  if (recalculate) {

    need_packages("pcaMethods")
    x <- readr::read_csv(input_path("metabolites", "metabolite_measurements.csv"),
      show_col_types = FALSE)
    x <- x[x$tissue == tissue, ]
    values <- as.matrix(x[, -(1:3)])

    # Unit-variance scaling and NIPALS are used in the original analysis script.
    values <- pcaMethods::prep(values, scale = "uv")
    model <- pcaMethods::pca(values, method = "nipals", nPcs = 2)
    scores <- data.frame(x[, c("genotype", "treatment")], pcaMethods::scores(model))

  } else {

    scores <- readr::read_table(input_path("metabolites",
      paste0(tissue, "_pca_scores.csv")), show_col_types = FALSE)
    parts <- strsplit(scores$group, "x", fixed = TRUE)
    scores$genotype <- tolower(vapply(parts, function(x) x[2], ""))
    scores$treatment <- vapply(parts, function(x) x[3], "")
    scores$genotype[scores$genotype == "aabb"] <- "mtn1121"

  }

  plot <- plot_scores(scores)
  ggsave(paste0(tissue, "_metabolites.pdf"), plot, width = 5, height = 4)
  readr::write_tsv(scores, paste0(tissue, "_metabolite_scores.tsv"))

}
