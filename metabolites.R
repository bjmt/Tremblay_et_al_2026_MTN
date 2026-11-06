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
aa_base <- hplc_aa[hplc_aa$genotype == "wt" & hplc_aa$treatment == "ms", ]
aa_base <- aa_base[, -2] %>%
  tidyr::gather("aa", "amount", 3:ncol(.)) %>%
  group_by(genotype, tissue, aa) %>%
  summarise(Mean = mean(amount, na.rm = TRUE))

hplc_fc <- hplc_aa[hplc_aa$genotype != "wt" & hplc_aa$treatment != "ms", ]
hplc_fc <- hplc_fc %>%
  tidyr::gather("aa", "amount", 4:ncol(.)) %>%
  group_by(genotype, treatment, tissue, aa) %>%
  summarise(Mean = mean(amount, na.rm = TRUE))
hplc_fc$FC <- NA_real_
for (i in 1:nrow(hplc_fc)) {
  i_met <- hplc_fc$aa[i]
  i_tiss <- hplc_fc$tissue[i]
  hplc_fc$FC[i] <- hplc_fc$Mean[i] /
    aa_base$Mean[aa_base$tissue == i_tiss & aa_base$aa == i_met]
}
hplc_fc <- hplc_fc[, -5] %>%
  tidyr::spread("aa", "FC")
readr::write_tsv(hplc_fc, "HPLC_AA_foldchange.tsv")

