# Plot HOMER gene and TE methylation profiles.
# Optional argument: the folder containing gene.txt and transposable_element.txt.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
args <- commandArgs(TRUE)
folder <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  file.path(output_dir, "preprocessing/bs_bud/profiles")
start_analysis("methylation_profiles")
palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
groups <- c("wm", "ws", "mm", "ms", "gm", "gs")
labels <- c("WT -Spd", "WT +Spd", "mtn1-1mtn2-1 -Spd",
            "mtn1-1mtn2-1 +Spd", "G3 -Spd", "G3 +Spd")
for (feature in c("gene", "transposable_element")) {
  prof <- read.delim(file.path(folder, paste0(feature, ".txt")),
                     check.names = FALSE, quote = "", comment.char = "")
  for (context in if (feature == "gene") "CG" else c("CG", "CHG", "CHH")) {
    # Plot rows 11:231: 1 kb upstream, the scaled body, 1 kb downstream.
    if (nrow(prof) < 231L) stop("Profile has fewer than 231 rows: ", feature)
    xi <- 11:231
    means <- sapply(groups, function(group) {
      expected <- paste0("/", group, 1:3, "_", context, "-tagDir Coverage")
      columns <- vapply(expected, function(suffix) {
        hit <- which(endsWith(names(prof), suffix))
        if (length(hit) != 1L) stop("Expected exactly one profile column ending in ", suffix)
        hit
      }, integer(1))
      rowMeans(prof[, columns, drop = FALSE])[xi]
    })
    prefix <- if (feature == "gene") "gene" else "te"
    pdf(paste0(prefix, "_", tolower(context), ".pdf"), width = 7, height = 3.5)
    par(mar = c(4, 4, 1, 10), xpd = NA)
    matplot(prof[[1]][xi], means, type = "l", lwd = 2, lty = 1,
            col = palette, bty = "n", xaxt = "n", xlab = NA,
            ylab = paste0("Fraction m", context), ylim = range(c(0, means)))
    axis(1, at = c(-1000, 0, 4000, 5000), labels = c("-1 kb", "Start", "End", "+1 kb"))
    usr <- par("usr")
    legend(usr[2] + .04 * diff(usr[1:2]), usr[4], xjust = 0, yjust = 1,
           bty = "n", lwd = 2, col = palette, legend = labels, cex = .65)
    dev.off()
    write.table(data.frame(position = prof[[1]][xi], means),
      paste0(prefix, "_", tolower(context), "_means.tsv"), sep = "\t",
      quote = FALSE, row.names = FALSE)
  }
}
