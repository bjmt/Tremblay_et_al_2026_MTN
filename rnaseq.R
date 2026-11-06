script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
args <- commandArgs(TRUE)
contrast_folder <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  file.path(data_dir, "rnaseq")
analysis_name <- if (length(args)) "rnaseq_recalculated" else "rnaseq"

start_analysis(analysis_name, c("dplyr", "ggplot2", "eulerr", "ComplexHeatmap", "matrixStats", "rtracklayer", "UpSetR", "readr", "circlize", "viridis"))
dir.create("gene_lists", showWarnings = FALSE)
dir.create("venns", showWarnings = FALSE)
library(dplyr)
library(ggplot2)
library(eulerr)
library(ComplexHeatmap)
library(matrixStats)

library(rtracklayer)

library(UpSetR)

plot_euler <- function(x) {
  colnames(x) <- sub("_(up|down)$", "\n\\1", colnames(x))
  plot(euler(x), quantities = list(cex = .8), labels = list(cex = .7))
}

bud <- readr::read_tsv(file.path(contrast_folder, "bud.tsv"), show_col_types = FALSE)
sdl <- readr::read_tsv(file.path(contrast_folder, "seedling.tsv"), show_col_types = FALSE)

if (ncol(bud) < 52 || ncol(sdl) < 46) {
  stop("Correlation plots require the GEO TPM matrices; rerun rnaseq_de.R with those files present")
}

# Missing tests are not significant. Use 1 only in the local plotting copies;
# leave NA values intact in the statistical result files.
for (column in grep("adj. p-value", names(bud), fixed = TRUE)) {
  bud[[column]][!is.finite(bud[[column]])] <- 1
}

for (column in grep("adj. p-value", names(sdl), fixed = TRUE)) {
  sdl[[column]][!is.finite(sdl[[column]])] <- 1
}

bud_msspd_mat_up <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT", "mtn", "G3")))
bud_msspd_mat_up[, 1][bud[["wt_ms vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["wt_ms vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_msspd_mat_up[, 2][bud[["mtn_ms vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. mtn_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_msspd_mat_up[, 3][bud[["g3_ms vs. g3_spd adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. g3_spd Log2 Fold Change"]] >= 1] <- TRUE

bud_msspd_mat_down <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT", "mtn", "G3")))
bud_msspd_mat_down[, 1][bud[["wt_ms vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["wt_ms vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_msspd_mat_down[, 2][bud[["mtn_ms vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. mtn_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_msspd_mat_down[, 3][bud[["g3_ms vs. g3_spd adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. g3_spd Log2 Fold Change"]] <= -1] <- TRUE

bud_ms_mat_up <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT_mtn", "WT_g3", "mtn_g3")))
bud_ms_mat_up[, 1][bud[["mtn_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. wt_ms Log2 Fold Change"]] <= -1] <- TRUE
bud_ms_mat_up[, 2][bud[["g3_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. wt_ms Log2 Fold Change"]] <= -1] <- TRUE
bud_ms_mat_up[, 3][bud[["g3_ms vs. mtn_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. mtn_ms Log2 Fold Change"]] <= -1] <- TRUE

bud_ms_mat_down <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT_mtn", "WT_g3", "mtn_g3")))
bud_ms_mat_down[, 1][bud[["mtn_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. wt_ms Log2 Fold Change"]] >= 1] <- TRUE
bud_ms_mat_down[, 2][bud[["g3_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. wt_ms Log2 Fold Change"]] >= 1] <- TRUE
bud_ms_mat_down[, 3][bud[["g3_ms vs. mtn_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. mtn_ms Log2 Fold Change"]] >= 1] <- TRUE

bud_spd_mat_up <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT_mtn", "WT_g3", "mtn_g3")))
bud_spd_mat_up[, 1][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat_up[, 2][bud[["g3_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat_up[, 3][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] <= -1] <- TRUE

bud_spd_mat_down <- matrix(FALSE, nrow = nrow(bud), ncol = 3,
  dimnames = list(bud$Gene, c("WT_mtn", "WT_g3", "mtn_g3")))
bud_spd_mat_down[, 1][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_spd_mat_down[, 2][bud[["g3_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_spd_mat_down[, 3][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] >= 1] <- TRUE

bud_spd_mat <- matrix(FALSE, nrow = nrow(bud), ncol = 6,
  dimnames = list(bud$Gene, c("WT_mtn_up", "WT_g3_up", "mtn_g3_up",
      "WT_mtn_down", "WT_g3_down", "mtn_g3_down")))
bud_spd_mat[, 1][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat[, 2][bud[["g3_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat[, 3][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat[, 4][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_spd_mat[, 5][bud[["g3_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_spd_mat[, 6][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] >= 1] <- TRUE

bud_ms_mat2 <- matrix(FALSE, nrow = nrow(bud), ncol = 4,
  dimnames = list(bud$Gene, c("WT_mtn_up", "mtn_g3_up",
      "WT_mtn_down", "mtn_g3_down")))
bud_ms_mat2[, 1][bud[["mtn_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. wt_ms Log2 Fold Change"]] <= -1] <- TRUE
bud_ms_mat2[, 2][bud[["g3_ms vs. mtn_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. mtn_ms Log2 Fold Change"]] <= -1] <- TRUE
bud_ms_mat2[, 3][bud[["mtn_ms vs. wt_ms adj. p-value"]] < 0.05 &
  bud[["mtn_ms vs. wt_ms Log2 Fold Change"]] >= 1] <- TRUE
bud_ms_mat2[, 4][bud[["g3_ms vs. mtn_ms adj. p-value"]] < 0.05 &
  bud[["g3_ms vs. mtn_ms Log2 Fold Change"]] >= 1] <- TRUE

bud_spd_mat2 <- matrix(FALSE, nrow = nrow(bud), ncol = 4,
  dimnames = list(bud$Gene, c("WT_mtn_up", "mtn_g3_up",
      "WT_mtn_down", "mtn_g3_down")))
bud_spd_mat2[, 1][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat2[, 2][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] <= -1] <- TRUE
bud_spd_mat2[, 3][bud[["mtn_spd vs. wt_spd adj. p-value"]] < 0.05 &
  bud[["mtn_spd vs. wt_spd Log2 Fold Change"]] >= 1] <- TRUE
bud_spd_mat2[, 4][bud[["g3_spd vs. mtn_spd adj. p-value"]] < 0.05 &
  bud[["g3_spd vs. mtn_spd Log2 Fold Change"]] >= 1] <- TRUE

pdf("Venn_Bud_MS_wt_mtn_g3.pdf", width = 3.7, height = 3.3)
plot_euler(bud_ms_mat2)
dev.off()

pdf("venns/MS_WT_mtn_up__mtn_g3_up.pdf", width = 3.7, height = 3.3)
plot_euler(bud_ms_mat2[, 1:2])
dev.off()
pdf("venns/MS_WT_mtn_down__mtn_g3_up.pdf", width = 3.7, height = 3.3)
plot_euler(bud_ms_mat2[, 3:2])
dev.off()
pdf("venns/MS_WT_mtn_up__mtn_g3_down.pdf", width = 3.7, height = 3.3)
plot_euler(bud_ms_mat2[, c(1, 4)])
dev.off()
pdf("venns/MS_WT_mtn_down__mtn_g3_down.pdf", width = 3.7, height = 3.3)
plot_euler(bud_ms_mat2[, c(3, 4)])
dev.off()

pdf("Venn_Bud_Spd_wt_mtn_g3.pdf", width = 3.7, height = 3.3)
plot_euler(bud_spd_mat2)
dev.off()

pdf("venns/Spd_WT_mtn_up__mtn_g3_up.pdf", width = 3.7, height = 3.3)
plot_euler(bud_spd_mat2[, 1:2])
dev.off()
pdf("venns/Spd_WT_mtn_up__mtn_g3_down.pdf", width = 3.7, height = 3.3)
plot_euler(bud_spd_mat2[, c(1, 4)])
dev.off()
pdf("venns/Spd_WT_mtn_down__mtn_g3_down.pdf", width = 3.7, height = 3.3)
plot_euler(bud_spd_mat2[, c(3, 4)])
dev.off()
pdf("venns/Spd_WT_mtn_down__mtn_g3_up.pdf", width = 3.7, height = 3.3)
plot_euler(bud_spd_mat2[, 3:2])
dev.off()

readr::write_lines(
  rownames(bud_spd_mat2)[bud_spd_mat2[, 1] & bud_spd_mat2[, 4] &
      !bud_spd_mat2[, 2] & !bud_spd_mat2[, 3]],
  "gene_lists/Bud_Spd_wtmtnUp_mtng3Down.txt")

readr::write_lines(
  rownames(bud_spd_mat2)[!bud_spd_mat2[, 1] & !bud_spd_mat2[, 4] &
      bud_spd_mat2[, 2] & bud_spd_mat2[, 3]],
  "gene_lists/Bud_Spd_wtmtnDown_mtng3Up.txt")

readr::write_lines(
  rownames(bud_ms_mat2)[bud_ms_mat2[, 1] & bud_ms_mat2[, 4] &
      !bud_ms_mat2[, 2] & !bud_ms_mat2[, 3]],
  "gene_lists/Bud_MS_wtmtnUp_mtng3Down.txt")

readr::write_lines(
  rownames(bud_ms_mat2)[!bud_ms_mat2[, 1] & !bud_ms_mat2[, 4] &
      bud_ms_mat2[, 2] & bud_ms_mat2[, 3]],
  "gene_lists/Bud_MS_wtmtnDown_mtng3Up.txt")

bud_m <- as.matrix(bud[, 29:52])
rownames(bud_m) <- bud$Gene
sdl_m <- as.matrix(sdl[, 29:46])
rownames(sdl_m) <- sdl$Gene

colnames(bud_m) <- gsub(" TPM", "", colnames(bud_m), TRUE)
colnames(sdl_m) <- gsub(" TPM", "", colnames(sdl_m), TRUE)

bud_selected <- head(order(rowVars(bud_m), decreasing = TRUE, na.last = NA), 1000)
sdl_selected <- head(order(rowVars(sdl_m), decreasing = TRUE, na.last = NA), 2000)

bud_pcc2 <- cor(scale(bud_m[bud_selected, c(17:24, 9:16, 1:8)]))
sdl_pcc2 <- cor(scale(sdl_m[sdl_selected,
  c(13, 15, 17, 14, 16, 18, 7, 9, 11, 8, 10, 12, 1, 3, 5, 2, 4, 6)]))

counts_file <- file.path(data_dir, "rnaseq/counts.csv")
have_counts <- file.exists(counts_file)
if (have_counts) {

  gcm <- readr::read_csv(counts_file)
  bud_pcc <- cor(scale(as.matrix(gcm[order(rowVars(as.matrix(gcm[, 2:25])),
          decreasing = TRUE)[1:1000], 2:25])))
  sdl_pcc <- cor(scale(as.matrix(gcm[order(rowVars(as.matrix(gcm[, c(38:43, 32:37, 26:31)])),
          decreasing = TRUE)[1:1000], c(38:43, 32:37, 26:31)])))

} else {
  message("Skipping count-based correlations: supply the GEO counts and rerun prepare_data.R.")
}

comparisons_bud <- list(
  "G3ms_G3spd" = c(2, 4),
  "G3ms_mtnms" = c(5, 7),
  "G3ms_WTms" = c(8, 10),
  "G3spd_mtnspd" = c(11, 13),
  "G3spd_WTspd" = c(14, 16),
  "mtnms_mtnspd" = c(17, 19),
  "mtnms_WTms" = c(20, 22),
  "mtnspd_WTspd" = c(23, 25),
  "WTms_WTspd" = c(26, 28)
)

comparisons_sdl <- list(
  "m1125ms_m1125spd" = c(2, 4),
  "m1125ms_m1121ms" = c(5, 7),
  "m1125ms_WTms" = c(8, 10),
  "m1125spd_m1121spd" = c(11, 13),
  "m1125spd_WTspd" = c(14, 16),
  "m1121ms_m1121spd" = c(17, 19),
  "m1121ms_WTms" = c(20, 22),
  "m1121spd_WTspd" = c(23, 25),
  "WTms_WTspd" = c(26, 28)
)

comparisons_bud_up_g <- comparisons_bud
comparisons_sdl_up_g <- comparisons_sdl
comparisons_bud_down_g <- comparisons_bud
comparisons_sdl_down_g <- comparisons_sdl

for (i in seq_along(comparisons_bud)) {
  comparisons_bud_up_g[[i]] <- bud$Gene[bud[[comparisons_bud[[i]][1]]] >= 1 &
    bud[[comparisons_bud[[i]][2]]] < 0.05]
  readr::write_lines(comparisons_bud_up_g[[i]],
    paste0("gene_lists/Bud_", names(comparisons_bud)[i], "_up.txt")
  )
  comparisons_bud_down_g[[i]] <- bud$Gene[bud[[comparisons_bud[[i]][1]]] <= -1 &
    bud[[comparisons_bud[[i]][2]]] < 0.05]
  readr::write_lines(comparisons_bud_down_g[[i]],
    paste0("gene_lists/Bud_", names(comparisons_bud)[i], "_down.txt")
  )
}

for (i in seq_along(comparisons_sdl)) {
  comparisons_sdl_up_g[[i]] <- sdl$Gene[sdl[[comparisons_sdl[[i]][1]]] >= 1 &
    sdl[[comparisons_sdl[[i]][2]]] < 0.05]
  readr::write_lines(comparisons_sdl_up_g[[i]],
    paste0("gene_lists/Seedling_", names(comparisons_sdl)[i], "_up.txt")
  )
  comparisons_sdl_down_g[[i]] <- sdl$Gene[sdl[[comparisons_sdl[[i]][1]]] <= -1 &
    sdl[[comparisons_sdl[[i]][2]]] < 0.05]
  readr::write_lines(comparisons_sdl_down_g[[i]],
    paste0("gene_lists/Seedling_", names(comparisons_sdl)[i], "_down.txt")
  )
}

# HOMER coefficients are second group / first group.
# Match the displayed mutant-versus-reference direction.
ups_bud_ms <- list(
  "mtn1-1mtn2-1 vs WT up-reg" = readr::read_lines("gene_lists/Bud_mtnms_WTms_down.txt"),
  "mtn1-1mtn2-1 vs WT down-reg" = readr::read_lines("gene_lists/Bud_mtnms_WTms_up.txt"),
  "G3 vs WT up-reg" = readr::read_lines("gene_lists/Bud_G3ms_WTms_down.txt"),
  "G3 vs WT down-reg" = readr::read_lines("gene_lists/Bud_G3ms_WTms_up.txt"),
  "G3 vs mtn1-1mtn2-1 up-reg" = readr::read_lines("gene_lists/Bud_G3ms_mtnms_down.txt"),
  "G3 vs mtn1-1mtn2-1 down-reg" = readr::read_lines("gene_lists/Bud_G3ms_mtnms_up.txt")
)
upset_plot <- upset(fromList(ups_bud_ms), length(ups_bud_ms), order.by="freq")#, empty.intersections="on")
pdf("upset_bud_ms.pdf", width = 6.5, height = 3.5)
print(upset_plot)
dev.off()

ups_bud_spd <- list(
  "WT -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Bud_WTms_WTspd_up.txt"),
  "WT -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Bud_WTms_WTspd_down.txt"),
  "mtn1-1mtn2-1 -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Bud_mtnms_mtnspd_up.txt"),
  "mtn1-1mtn2-1 -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Bud_mtnms_mtnspd_down.txt"),
  "G3 -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Bud_G3ms_G3spd_up.txt"),
  "G3 -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Bud_G3ms_G3spd_down.txt")
)
upset_plot <- upset(fromList(ups_bud_spd), length(ups_bud_spd), order.by="freq")
pdf("upset_bud_spd.pdf", width = 5.0, height = 3.5)
print(upset_plot)
dev.off()

# HOMER coefficients are second group / first group.
# Match the displayed mutant-versus-reference direction.
ups_sdl_ms <- list(
  "mtn1-1mtn2-1 vs WT up-reg" = readr::read_lines("gene_lists/Seedling_m1121ms_WTms_down.txt"),
  "mtn1-1mtn2-1 vs WT down-reg" = readr::read_lines("gene_lists/Seedling_m1121ms_WTms_up.txt"),
  "mtn1-1mtn2-5 vs WT up-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_WTms_down.txt"),
  "mtn1-1mtn2-5 vs WT down-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_WTms_up.txt"),
  "mtn1-1mtn2-5 vs mtn1-1mtn2-1 up-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_m1121ms_down.txt"),
  "mtn1-1mtn2-5 vs mtn1-1mtn2-1 down-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_m1121ms_up.txt")
)
upset_plot <- upset(fromList(ups_sdl_ms), length(ups_sdl_ms), order.by="freq")
pdf("upset_seedling_ms.pdf", width = 7.0, height = 3.5)
print(upset_plot)
dev.off()

ups_sdl_spd <- list(
  "WT -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Seedling_WTms_WTspd_up.txt"),
  "WT -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Seedling_WTms_WTspd_down.txt"),
  "mtn1-1mtn2-1 -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Seedling_m1121ms_m1121spd_up.txt"),
  "mtn1-1mtn2-1 -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Seedling_m1121ms_m1121spd_down.txt"),
  "mtn1-1mtn2-5 -Spd vs +Spd up-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_m1125spd_up.txt"),
  "mtn1-1mtn2-5 -Spd vs +Spd down-reg" = readr::read_lines("gene_lists/Seedling_m1125ms_m1125spd_down.txt")
)
upset_plot <- upset(fromList(ups_sdl_spd), length(ups_sdl_spd), order.by="freq")
pdf("upset_seedling_spd.pdf", width = 6.0, height = 3.5)
print(upset_plot)
dev.off()

if (have_counts) {

  pdf("Bud_RNAseq_PCC.pdf", width = 5, height = 4)
  ComplexHeatmap::draw(Heatmap(bud_pcc, name = "Pearson R",
    col = circlize::colorRamp2(c(0.6, 0.7, 0.8, 0.9, 1), viridis::inferno(5)),
    cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = FALSE, show_row_names = FALSE,
    row_split = c("WT\n-Spd", "WT\n+Spd", "mm\n-Spd", "mm\n+Spd", "G3\n-Spd", "G3\n+Spd")[
      rep(1:6, each = 4)],
    column_split = c("WT\n-Spd", "WT\n+Spd", "mm\n-Spd", "mm\n+Spd", "G3\n-Spd", "G3\n+Spd")[
      rep(1:6, each = 4)]
  ))
  dev.off()

  pdf("Seedling_RNAseq_PCC.pdf", width = 5, height = 4)
  ComplexHeatmap::draw(Heatmap(sdl_pcc, name = "Pearson R",
    col = circlize::colorRamp2(c(0.6, 0.7, 0.8, 0.9, 1), viridis::inferno(5)),
    cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = FALSE, show_row_names = FALSE,
    row_split = c("WT\n-Spd", "WT\n+Spd", "m5\n-Spd", "m5\n+Spd", "mm\n-Spd", "mm\n+Spd")[
      rep(1:6, each = 3)],
    column_split = c("WT\n-Spd", "WT\n+Spd", "m5\n-Spd", "m5\n+Spd", "mm\n-Spd", "mm\n+Spd")[
      rep(1:6, each = 3)]
  ))
  dev.off()

}

pdf("Bud_RNAseq_PCC_tpm.pdf", width = 5, height = 4)
ComplexHeatmap::draw(Heatmap(bud_pcc2, name = "Pearson R",
  col = circlize::colorRamp2(c(0.7, 0.8, 0.9, 1), viridis::inferno(4)),
  cluster_columns = FALSE, cluster_rows = FALSE,
  show_column_names = FALSE, show_row_names = FALSE,
  row_split = c("WT\n-Spd", "WT\n+Spd", "mm\n-Spd", "mm\n+Spd", "G3\n-Spd", "G3\n+Spd")[
    rep(1:6, each = 4)],
  column_split = c("WT\n-Spd", "WT\n+Spd", "mm\n-Spd", "mm\n+Spd", "G3\n-Spd", "G3\n+Spd")[
    rep(1:6, each = 4)]
))
dev.off()

pdf("Seedling_RNAseq_PCC_tpm.pdf", width = 5, height = 4)
ComplexHeatmap::draw(Heatmap(sdl_pcc2, name = "Pearson R",
  col = circlize::colorRamp2(c(0.7, 0.8, 0.9, 1), viridis::inferno(4)),
  cluster_columns = FALSE, cluster_rows = FALSE,
  show_column_names = FALSE, show_row_names = FALSE,
  row_split = c("WT\n-Spd", "WT\n+Spd", "m5\n-Spd", "m5\n+Spd", "mm\n-Spd", "mm\n+Spd")[
    rep(1:6, each = 3)],
  column_split = c("WT\n-Spd", "WT\n+Spd", "m5\n-Spd", "m5\n+Spd", "mm\n-Spd", "mm\n+Spd")[
    rep(1:6, each = 3)]
))
dev.off()

annotation <- file.path(data_dir, "reference/Araport11_GFF3_genes_transposons.gff3")
if (file.exists(annotation)) {

  a11 <- import(file.path(data_dir, "reference/Araport11_GFF3_genes_transposons.gff3"))

  g2t <- structure(a11$locus_type[!is.na(a11$locus_type)], names = a11$ID[!is.na(a11$locus_type)])

  te_wt_mtn_sdl <- data.frame(
    Gene = sdl$Gene,
    IsTE = g2t[sdl$Gene] == "transposable_element_gene",
    FC = -sdl[["mm vs. wm Log2 Fold Change"]],
    Pval = sdl[["mm vs. wm adj. p-value"]]
  )
  te_wt_m5_sdl <- data.frame(
    Gene = sdl$Gene,
    IsTE = g2t[sdl$Gene] == "transposable_element_gene",
    FC = -sdl[["5m vs. wm Log2 Fold Change"]],
    Pval = sdl[["5m vs. wm adj. p-value"]]
  )
  te_wt_mtn_bud <- data.frame(
    Gene = bud$Gene,
    IsTE = g2t[bud$Gene] == "transposable_element_gene",
    FC = -bud[["mtn_ms vs. wt_ms Log2 Fold Change"]],
    Pval = bud[["mtn_ms vs. wt_ms adj. p-value"]]
  )
  te_wt_g3_bud <- data.frame(
    Gene = bud$Gene,
    IsTE = g2t[bud$Gene] == "transposable_element_gene",
    FC = -bud[["g3_ms vs. wt_ms Log2 Fold Change"]],
    Pval = bud[["g3_ms vs. wt_ms adj. p-value"]]
  )

  te_wt_mtn_sdl$FC[te_wt_mtn_sdl$FC > 10] <- 10
  te_wt_mtn_sdl$FC[te_wt_mtn_sdl$FC < -10] <- -10
  te_wt_mtn_sdl$Pval[te_wt_mtn_sdl$Pval < 10^-50] <- 10^-50
  with(te_wt_mtn_sdl, {
    print(sum(IsTE & Pval < .05 & FC >= 1) / sum(Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1) / sum(Pval < .05 & FC <= -1))
    print(sum(IsTE & Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1))
  })
  with(te_wt_mtn_sdl, {
    pdf("TE_genes_Seedlings_WTms_mtnms.pdf", width = 4.2, height = 3.6)
    plot(FC, -log10(Pval), pch = 19, cex = 0.4,
      bty = "n", ylab = expression(-log[10](Pvalue)),
      main = expression(paste("WT vs. ", italic("mtn1-1mtn2-1"), " (Seedlings)")),
      xlab = expression(log[2](FC)), col = "grey"
    )
    abline(v = c(-1, 1), lty = 2, col = "grey")
    abline(h = -log10(0.05), lty = 2, col = "grey")
    points(FC[Pval < .05 & abs(FC) >= 1],
      -log10(Pval[Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "black")
    points(FC[IsTE & Pval < .05 & abs(FC) >= 1],
      -log10(Pval[IsTE & Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "red")
    legend(-10, 45, legend = c("TE gene", "Sig. DE", "N.S."), col = c("red", "black", "grey"),
      pch = 19, cex = 0.8, bty = "n")
    dev.off()
  })

  te_wt_mtn_bud$FC[te_wt_mtn_bud$FC > 10] <- 10
  te_wt_mtn_bud$FC[te_wt_mtn_bud$FC < -10] <- -10
  te_wt_mtn_bud$Pval[te_wt_mtn_bud$Pval < 10^-50] <- 10^-50
  with(te_wt_mtn_bud, {
    print(sum(IsTE & Pval < .05 & FC >= 1) / sum(Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1) / sum(Pval < .05 & FC <= -1))
    print(sum(IsTE & Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1))
  })
  with(te_wt_mtn_bud, {
    pdf("TE_genes_Buds_WTms_mtnms.pdf", width = 4.2, height = 3.6)
    plot(FC, -log10(Pval), pch = 19, cex = 0.4,
      bty = "n", ylab = expression(-log[10](Pvalue)),
      main = expression(paste("WT vs. ", italic("mtn1-1mtn2-1"), " (Buds)")),
      xlab = expression(log[2](FC)), col = "grey"
    )
    abline(v = c(-1, 1), lty = 2, col = "grey")
    abline(h = -log10(0.05), lty = 2, col = "grey")
    points(FC[Pval < .05 & abs(FC) >= 1],
      -log10(Pval[Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "black")
    points(FC[IsTE & Pval < .05 & abs(FC) >= 1],
      -log10(Pval[IsTE & Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "red")
    legend(-10, 45, legend = c("TE gene", "Sig. DE", "N.S."), col = c("red", "black", "grey"),
      pch = 19, cex = 0.8, bty = "n")
    dev.off()
  })

  te_wt_m5_sdl$FC[te_wt_m5_sdl$FC > 10] <- 10
  te_wt_m5_sdl$FC[te_wt_m5_sdl$FC < -10] <- -10
  te_wt_m5_sdl$Pval[te_wt_m5_sdl$Pval < 10^-50] <- 10^-50
  with(te_wt_m5_sdl, {
    print(sum(IsTE & Pval < .05 & FC >= 1) / sum(Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1) / sum(Pval < .05 & FC <= -1))
    print(sum(IsTE & Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1))
  })
  with(te_wt_m5_sdl, {
    pdf("TE_genes_Seedlings_WTms_m5ms.pdf", width = 4.2, height = 3.6)
    plot(FC, -log10(Pval), pch = 19, cex = 0.4,
      bty = "n", ylab = expression(-log[10](Pvalue)),
      main = expression(paste("WT vs. ", italic("mtn1-1mtn2-5"), " (Seedlings)")),
      xlab = expression(log[2](FC)), col = "grey"
    )
    abline(v = c(-1, 1), lty = 2, col = "grey")
    abline(h = -log10(0.05), lty = 2, col = "grey")
    points(FC[Pval < .05 & abs(FC) >= 1],
      -log10(Pval[Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "black")
    points(FC[IsTE & Pval < .05 & abs(FC) >= 1],
      -log10(Pval[IsTE & Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "red")
    legend(-10, 45, legend = c("TE gene", "Sig. DE", "N.S."), col = c("red", "black", "grey"),
      pch = 19, cex = 0.8, bty = "n")
    dev.off()
  })

  te_wt_g3_bud$FC[te_wt_g3_bud$FC > 10] <- 10
  te_wt_g3_bud$FC[te_wt_g3_bud$FC < -10] <- -10
  te_wt_g3_bud$Pval[te_wt_g3_bud$Pval < 10^-50] <- 10^-50
  with(te_wt_g3_bud, {
    print(sum(IsTE & Pval < .05 & FC >= 1) / sum(Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1) / sum(Pval < .05 & FC <= -1))
    print(sum(IsTE & Pval < .05 & FC >= 1))
    print(sum(IsTE & Pval < .05 & FC <= -1))
  })
  with(te_wt_g3_bud, {
    pdf("TE_genes_Buds_WTms_G3ms.pdf", width = 4.2, height = 3.6)
    plot(FC, -log10(Pval), pch = 19, cex = 0.4,
      bty = "n", ylab = expression(-log[10](Pvalue)),
      main = expression("WT vs. G3 (Buds)"),
      xlab = expression(log[2](FC)), col = "grey"
    )
    abline(v = c(-1, 1), lty = 2, col = "grey")
    abline(h = -log10(0.05), lty = 2, col = "grey")
    points(FC[Pval < .05 & abs(FC) >= 1],
      -log10(Pval[Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "black")
    points(FC[IsTE & Pval < .05 & abs(FC) >= 1],
      -log10(Pval[IsTE & Pval < .05 & abs(FC) >= 1]), pch = 19, cex = 0.4, col = "red")
    legend(-10, 45, legend = c("TE gene", "Sig. DE", "N.S."), col = c("red", "black", "grey"),
      pch = 19, cex = 0.8, bty = "n")
    dev.off()
  })

} else {
  message("Skipping TE volcano plots: add the Araport11 GFF3 under reference/.")
}
