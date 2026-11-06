script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("metabolites", c("ggplot2", "dplyr", "readr", "tidyr", "rstatix", "multcompView", "cowplot"))
dir.create("pvalues", showWarnings = FALSE)
library(ggplot2)
library(dplyr)

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

hplc <- readr::read_csv(input_path("metabolites", "metabolite_measurements.csv"))
hplc$genotype[hplc$genotype == "mtn1121 (B)"] <- "mtn1121"
hplc$genotype[hplc$genotype == "mtn1121 (A)"] <- "mtn1121"

hplc_bud <- hplc[hplc$tissue == "bud", ]
str(cor.test(hplc_bud$SAH, hplc_bud$MTA))
str(cor.test(hplc_bud$Met, hplc_bud$MTA))
sort(cor(as.data.frame(hplc_bud[, -(1:3)]), use="pairwise.complete.obs")[, "MTA"], decreasing=TRUE)

mtapcc <- data.frame(Metabolite = colnames(hplc_bud[, -(1:3)]), Corr = NA, Pvalue = NA)
for (i in 1:nrow(mtapcc)) {
  i_ <- cor.test(hplc_bud[[mtapcc$Metabolite[i]]], hplc_bud$MTA)
  mtapcc$Corr[i] <- i_$estimate
  mtapcc$Pvalue[i] <- i_$p.value
}
mtapcc <- mtapcc[mtapcc$Metabolite != "MTA", ]
mtapcc <- mtapcc[order(mtapcc$Corr, decreasing = TRUE), ]
readr::write_tsv(mtapcc, "MTA_bud_correlations.tsv")

hplc_sld <- hplc[hplc$tissue == "seedling", ]
sort(cor(as.data.frame(hplc_sld[, -(1:3)]), use="pairwise.complete.obs")[, "MTA"], decreasing=TRUE)

hplc_bud <- hplc[hplc$tissue == "bud", -3]
hplc_sdl <- hplc[hplc$tissue == "seedling", -3]

hplc_bud2 <- hplc_bud %>%
  group_by(genotype, treatment) %>%
  tidyr::gather("metabolite", "amount", 3:31)
hplc_bud2$GenoTreat <- paste0(hplc_bud2$genotype, ":", hplc_bud2$treatment)

hplc_sdl2 <- hplc_sdl %>%
  group_by(genotype, treatment) %>%
  tidyr::gather("metabolite", "amount", 3:31)
hplc_sdl2$GenoTreat <- paste0(hplc_sdl2$genotype, ":", hplc_sdl2$treatment)

mets <- colnames(hplc)[-(1:3)]

for (i in seq_along(mets)) {
  res <- rstatix::dunn_test(hplc_bud2[hplc_bud2$metabolite == mets[i], 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
  readr::write_tsv(res,
    paste0("pvalues/HPLC_buds_", mets[i], ".tsv"))
  res <- rstatix::dunn_test(hplc_sdl2[hplc_sdl2$metabolite == mets[i], 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
  readr::write_tsv(res,
    paste0("pvalues/HPLC_seedlings_", mets[i], ".tsv"))
}

prep_data <- function(x) {
  group_by(x[, c("genotype", "treatment", "amount")], genotype, treatment) %>%
    summarise(
      lower = quantile(amount, .25, na.rm = TRUE),
      mid = quantile(amount, .5, na.rm = TRUE),
      upper = quantile(amount, .75, na.rm = TRUE)
    )
}

tt <- rstatix::dunn_test(hplc_bud2[hplc_bud2$metabolite == "Cys", 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
multcompView::multcompLetters(structure(tt$p.adj, names=paste0(tt$group1,"-",tt$group2)))

tt <- rstatix::dunn_test(hplc_sdl2[hplc_sdl2$metabolite == "Spd", 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
multcompView::multcompLetters(structure(tt$p.adj, names=paste0(tt$group1,"-",tt$group2)))

tt <- rstatix::dunn_test(hplc_sdl2[hplc_sdl2$metabolite == "Put", 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
multcompView::multcompLetters(structure(tt$p.adj, names=paste0(tt$group1,"-",tt$group2)))

tt <- rstatix::dunn_test(hplc_bud2[hplc_bud2$metabolite == "Put", 5:4],
    amount ~ GenoTreat, p.adjust.method = "BH")
multcompView::multcompLetters(structure(tt$p.adj, names=paste0(tt$group1,"-",tt$group2)))

plot_met <- function(x, ylim2 = NA, title = NULL) {
  x$genotype <- factor(x$genotype, levels = rev(c("wt", "mtn1125", "mtn1121", "g2", "g3")))
  x$treatment <- factor(x$treatment, levels = rev(c("ms", "spd")))
  x_summ <- prep_data(x)
  ggplot() +
    geom_point(aes(genotype, amount, colour = treatment, group = treatment), data = x,
        size = 0.85, position = position_jitterdodge(jitter.width = 0.08)) +
      geom_errorbar(aes(genotype, ymin = lower, ymax = upper, group = treatment),
        data = x_summ, width = 0.4, position = position_dodge(0.75)) +
      geom_errorbar(aes(genotype, ymin = mid, ymax = mid, group = treatment),
        data = x_summ, width = 0.5, position = position_dodge(0.75)) +
      scale_colour_manual(values = rev(c("black", "brown")),
        labels = rev(c("-Spd", "+Spd")), name = element_blank()) +
      scale_x_discrete(labels = rev(c("WT", bquote(italic("mtn1-1mtn2-5")),
          bquote(italic("mtn1-1mtn2-1")), "G2", "G3"))) +
      ylim(c(0, ylim2)) +
      xlab(NULL) +
      ylab(expression(paste("pmol ", mg^{-1}))) +
      coord_flip() +
      theme_bw() +
      ggtitle(title) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.75),
        rect = element_rect(size = 0.9),
        axis.ticks = element_line(size = 0.75),
        text = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 6),
        axis.text.x = element_text(colour = "black", size = 6),
        panel.grid = element_blank()
      )
}

cowplot::save_plot(
  "metabolites.pdf",
  base_width = 5, base_height = 3,
  cowplot::plot_grid(ncol = 2, byrow = TRUE,
    plot_met(hplc_bud2[hplc_bud2$metabolite == "MTA", ], 175, "MTA (Buds)"),
    plot_met(hplc_sdl2[hplc_sdl2$metabolite == "MTA", ], 28, "MTA (Seedlings)"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "Spd", ], 700, "Spd (Buds)"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "NIC", ], 1000, "NA (Buds)")
  )
)

cowplot::save_plot(
  "seedling_metabolites.pdf",
  base_width = 5, base_height = 1.5,
  cowplot::plot_grid(ncol = 2, byrow = TRUE,
    plot_met(hplc_sdl2[hplc_sdl2$metabolite == "Spd", ], 230, "Spd (Seedlings)"),
    plot_met(hplc_sdl2[hplc_sdl2$metabolite == "NIC", ], 475, "NA (Seedlings)")
  )
)

cowplot::save_plot(
  "putrescine.pdf",
  base_width = 5, base_height = 1.5,
  cowplot::plot_grid(ncol = 2, byrow = TRUE,
    plot_met(hplc_sdl2[hplc_sdl2$metabolite == "Put", ], 90, "Put (Seedlings)"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "Put", ], 1000, "Put (Buds)")
  )
)

cowplot::save_plot(
  "sulfur_metabolites.pdf",
  base_width = 5, base_height = 4.5,
  cowplot::plot_grid(ncol = 2, byrow = TRUE,
    plot_met(hplc_bud2[hplc_bud2$metabolite == "SUL", ], 8, "Sulfate"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "SAM", ], 110, "SAM"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "SAH", ], 12, "SAH"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "Met", ], 50, "Met"),
    plot_met(hplc_bud2[hplc_bud2$metabolite == "Cys", ], 50, "Cys")
  )
)

aa <- sort(c("Cys", "Asp", "Glu", "Asn", "Ser", "Gln", "Gly", "His",
    "Arg", "Thr", "Ala", "Pro", "Tyr", "Val", "Met", "Ile", "Leu", "Lys",
    "Phe"))

hplc_aa <- hplc[, c("genotype", "treatment", "tissue", aa)]

# Each genotype is compared with WT of the same tissue and Spd treatment.

means <- hplc_aa %>%
  tidyr::pivot_longer(all_of(aa), names_to = "AA", values_to = "Amount") %>%
  group_by(genotype, treatment, tissue, AA) %>%
  summarise(Mean = mean(Amount, na.rm = TRUE), .groups = "drop")

control <- means[means$genotype == "wt", c("treatment", "tissue", "AA", "Mean")]
names(control)[4] <- "WT_mean"

fold_changes <- means[means$genotype %in% c("mtn1121", "g3"), ] %>%
  left_join(control, by = c("treatment", "tissue", "AA")) %>%
  mutate(FC = Mean / WT_mean, log2FC = log2(FC))

if (any(!is.finite(fold_changes$log2FC))) {
  stop("Non-finite amino-acid fold changes: check the WT denominators")
}

readr::write_tsv(fold_changes, "amino_acid_fold_changes.tsv")

fold_changes$Group <- with(fold_changes, paste(genotype, treatment, tissue, sep = ":"))
fold_changes$AA <- factor(fold_changes$AA, levels = aa)

plot <- ggplot(fold_changes, aes(AA, Group, colour = log2FC, size = abs(log2FC))) +
  geom_point() +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "log2 fold change") +
  scale_size_area(max_size = 6, name = "Absolute log2 fold change") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("amino_acid_fold_changes.pdf", plot, width = 10, height = 4)
