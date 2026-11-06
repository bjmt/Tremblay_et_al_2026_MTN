script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("tor", c("ggplot2", "dplyr", "readr", "tidyr", "rstatix", "multcompView", "cowplot"))
dir.create("pvalues", showWarnings = FALSE)
library(ggplot2)
library(dplyr)

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tor <- readr::read_csv2(input_path("zenodo", "TOR measurements.csv"))
tor$S6K <- as.numeric(tor$S6K)
tor$S6Kp <- as.numeric(tor$S6Kp)
tor$Ratio <- as.numeric(tor$Ratio)
tor$genotype <- factor(tor$genotype, levels = rev(c("WT", "mtn", "G3")))
tor$treatment <- factor(tor$treatment, levels = rev(c("MS", "Spd")))
tor$GenoTreat <- paste0(as.character(tor$genotype), "_", as.character(tor$treatment))

ttor <- rstatix::dunn_test(tor, Ratio ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(ttor, "tor_ratio_pvalues.tsv")
multcompView::multcompLetters(structure(ttor$p.adj, names=paste0(ttor$group1,"-",ttor$group2)))

ttor2 <- rstatix::dunn_test(tor, S6K ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(ttor2, "tor_s6k_pvalues.tsv")
multcompView::multcompLetters(structure(ttor2$p.adj, names=paste0(ttor2$group1,"-",ttor2$group2)))
ttor3 <- rstatix::dunn_test(tor, S6Kp ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(ttor3, "tor_s6kp_pvalues.tsv")
multcompView::multcompLetters(structure(ttor3$p.adj, names=paste0(ttor3$group1,"-",ttor3$group2)))

tor_s <- tor %>%
  group_by(genotype, treatment) %>%
  summarise(lower = quantile(Ratio, .25, na.rm=T),
    mid = quantile(Ratio, .5, na.rm=T),
    upper = quantile(Ratio, .75, na.rm = TRUE)
  )

tor_s6k <- tor %>%
  group_by(genotype, treatment) %>%
  summarise(lower = quantile(S6K, .25, na.rm=T),
    mid = quantile(S6K, .5, na.rm=T),
    upper = quantile(S6K, .75, na.rm = TRUE)
  )
tor_s6kp <- tor %>%
  group_by(genotype, treatment) %>%
  summarise(lower = quantile(S6Kp, .25, na.rm=T),
    mid = quantile(S6Kp, .5, na.rm=T),
    upper = quantile(S6Kp, .75, na.rm = TRUE)
  )

p1 <- ggplot() +
    geom_point(aes(genotype, Ratio, colour = treatment, group = treatment), data = tor,
        size = 0.85, position = position_jitterdodge(jitter.width = 0.08)) +
      geom_errorbar(aes(genotype, ymin = lower, ymax = upper, group = treatment),
        data = tor_s, width = 0.4, position = position_dodge(0.75)) +
      geom_errorbar(aes(genotype, ymin = mid, ymax = mid, group = treatment),
        data = tor_s, width = 0.5, position = position_dodge(0.75)) +
      scale_colour_manual(values = rev(c("black", "brown")),
        labels = rev(c("-Spd", "+Spd")), name = element_blank()) +
      scale_x_discrete(labels = rev(c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3"))) +
      ylim(c(0, NA)) +
      xlab(NULL) +
      ylab("S6K-p / S6K ratio") +
      coord_flip() +
      theme_bw() +
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

p1.1 <- ggplot() +
    geom_point(aes(genotype, S6K, colour = treatment, group = treatment), data = tor,
        size = 0.85, position = position_jitterdodge(jitter.width = 0.08)) +
      geom_errorbar(aes(genotype, ymin = lower, ymax = upper, group = treatment),
        data = tor_s6k, width = 0.4, position = position_dodge(0.75)) +
      geom_errorbar(aes(genotype, ymin = mid, ymax = mid, group = treatment),
        data = tor_s6k, width = 0.5, position = position_dodge(0.75)) +
      scale_colour_manual(values = rev(c("black", "brown")),
        labels = rev(c("-Spd", "+Spd")), name = element_blank()) +
      scale_x_discrete(labels = rev(c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3"))) +
      ylim(c(0, NA)) +
      xlab(NULL) +
      ylab("Relative S6K") +
      coord_flip() +
      theme_bw() +
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

p1.2 <- ggplot() +
    geom_point(aes(genotype, S6Kp, colour = treatment, group = treatment), data = tor,
        size = 0.85, position = position_jitterdodge(jitter.width = 0.08)) +
      geom_errorbar(aes(genotype, ymin = lower, ymax = upper, group = treatment),
        data = tor_s6kp, width = 0.4, position = position_dodge(0.75)) +
      geom_errorbar(aes(genotype, ymin = mid, ymax = mid, group = treatment),
        data = tor_s6kp, width = 0.5, position = position_dodge(0.75)) +
      scale_colour_manual(values = rev(c("black", "brown")),
        labels = rev(c("-Spd", "+Spd")), name = element_blank()) +
      scale_x_discrete(labels = rev(c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3"))) +
      ylim(c(0, NA)) +
      xlab(NULL) +
      ylab("Relative S6K-p") +
      coord_flip() +
      theme_bw() +
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

endo <- readr::read_tsv(input_path("zenodo", "Endoreduplication index.txt"))

endo <- tidyr::gather(endo, "GenoTreat", "Index", 1:6)
endo$genotype <- vapply(strsplit(endo$GenoTreat, "_", fixed=T), function(x) x[1], "")
endo$treatment <- vapply(strsplit(endo$GenoTreat, "_", fixed=T), function(x) x[2], "")

tendo <- rstatix::dunn_test(endo, Index ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(tendo, "endoredup_index_pvalues.tsv")
multcompView::multcompLetters(structure(tendo$p.adj, names=paste0(tendo$group1,"-",tendo$group2)))

endo$genotype <- factor(endo$genotype, levels = rev(c("WT", "dm", "G3")))
endo$treatment <- factor(endo$treatment, levels = rev(c("MS", "Spd")))

endo_s <- endo %>%
  group_by(genotype, treatment) %>%
  summarise(
    lower = quantile(Index, .25, na.rm=T),
    mid = quantile(Index, .5, na.rm=T),
    upper = quantile(Index, .75, na.rm = TRUE)
  )

p2 <- ggplot() +
    geom_point(aes(genotype, Index, colour = treatment, group = treatment), data = endo,
        size = 0.85, position = position_jitterdodge(jitter.width = 0.08)) +
      geom_errorbar(aes(genotype, ymin = lower, ymax = upper, group = treatment),
        data = endo_s, width = 0.4, position = position_dodge(0.75)) +
      geom_errorbar(aes(genotype, ymin = mid, ymax = mid, group = treatment),
        data = endo_s, width = 0.5, position = position_dodge(0.75)) +
      scale_colour_manual(values = rev(c("black", "brown")),
        labels = rev(c("-Spd", "+Spd")), name = element_blank()) +
      scale_x_discrete(labels = rev(c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3"))) +
      ylim(c(0, NA)) +
      xlab(NULL) +
      ylab("Endoreduplication index") +
      coord_flip() +
      theme_bw() +
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

cowplot::save_plot("tor_and_endoreduplication.pdf",
  base_width = 4.5, base_height = 1.2,
  cowplot::plot_grid(ncol = 2, p1, p2))

cowplot::save_plot("tor_protein_levels.pdf",
  base_width = 4.5, base_height = 1.2,
  cowplot::plot_grid(ncol = 2, p1.1, p1.2))

