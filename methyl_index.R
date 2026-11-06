script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("methyl_index", c("ggplot2", "dplyr", "readr", "tidyr", "rstatix", "multcompView", "cowplot"))

library(ggplot2)
library(cowplot)
library(dplyr)

met <- readr::read_csv(input_path("metabolites", "metabolite_measurements.csv"))
met$tissue <- factor(met$tissue, levels = c("bud", "seedling"))
levels(met$tissue) <- c("Bud", "Seedling")
met$genotype[grepl("mtn1121", met$genotype)] <- "mtn1121"
met$treatment <- factor(met$treatment, levels = c("spd", "ms"))
levels(met$treatment) <- c("+Spd", "-Spd")
met <- met[, c("genotype", "treatment", "tissue", "SAM", "SAH")]
met$Ratio <- met$SAM / met$SAH

met2 <- met %>%
  select(-SAM, -SAH) %>%
  group_by(genotype, treatment, tissue) %>%
  summarise(
    Low = quantile(Ratio, .25, na.rm = TRUE),
    Mid = quantile(Ratio, .5, na.rm = TRUE),
    High = quantile(Ratio, .75, na.rm = TRUE)
  )

met_s <- filter(met, tissue == "Seedling")
met2_s <- filter(met2, tissue == "Seedling")
met_b <- filter(met, tissue == "Bud")
met2_b <- filter(met2, tissue == "Bud")

met_s$GenoTreat <- paste0(met_s$genotype, "_", gsub("-", "_", met_s$treatment, fixed=T))
met_b$GenoTreat <- paste0(met_b$genotype, "_", gsub("-", "_", met_b$treatment, fixed=T))

t1 <- rstatix::dunn_test(met_s, Ratio ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(t1, "methyl_index_seedling_pvalues.tsv")
multcompView::multcompLetters(structure(t1$p.adj, names=paste0(t1$group1,"-",t1$group2)))

t2 <- rstatix::dunn_test(met_b, Ratio ~ GenoTreat, p.adjust.method = "BH")
readr::write_tsv(t2, "methyl_index_bud_pvalues.tsv")
multcompView::multcompLetters(structure(t2$p.adj, names=paste0(t2$group1,"-",t2$group2)))

p1 <- ggplot() +
  geom_point(aes(genotype, Ratio, colour = treatment),
    data = met_s, size = 1, position = position_dodge(0.7)) +
  geom_errorbar(aes(genotype, ymin = Low, ymax = High, colour = treatment),
    data = met2_s, width = 0.7, size = 1, position = position_dodge(0.7)) +
  geom_errorbar(aes(genotype, ymin = Mid, ymax = Mid, colour = treatment),
    data = met2_s, width = 0.5, size = 1, position = position_dodge(0.7)) +
  scale_colour_manual(values = c("brown", "black"), name = element_blank()) +
  scale_x_discrete(limits = c("g3", "g2", "mtn1121", "mtn1125", "wt"),
    labels = c("G3", "G2", bquote(italic("mtn1-1\nmtn2-1")), bquote(italic("mtn1-1\nmtn2-5")), "WT")) +
  xlab(NULL) +
  ylab("SAM/SAH ratio") +
  ylim(c(0, 40)) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 8, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.y.left = element_line(size = 0.5),
    axis.line.x.bottom = element_line(size = 0.5)
  )

p2 <- ggplot() +
  geom_point(aes(genotype, Ratio, colour = treatment),
    data = met_b, size = 1, position = position_dodge(0.7)) +
  geom_errorbar(aes(genotype, ymin = Low, ymax = High, colour = treatment),
    data = met2_b, width = 0.7, size = 1, position = position_dodge(0.7)) +
  geom_errorbar(aes(genotype, ymin = Mid, ymax = Mid, colour = treatment),
    data = met2_b, width = 0.5, size = 1, position = position_dodge(0.7)) +
  scale_colour_manual(values = c("brown", "black"), name = element_blank()) +
  scale_x_discrete(limits = c("g3", "g2", "mtn1121", "mtn1125", "wt"),
    labels = c("G3", "G2", bquote(italic("mtn1-1\nmtn2-1")), bquote(italic("mtn1-1\nmtn2-5")), "WT")) +
  xlab(NULL) +
  ylab("SAM/SAH ratio") +
  ylim(c(0, 40)) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 8),
    axis.title = element_text(size = 8, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.y.left = element_line(size = 0.5),
    axis.line.x.bottom = element_line(size = 0.5)
  )

cowplot::save_plot("methyl_index.pdf",
  cowplot::plot_grid(p1, p2, nrow = 1),
  base_height = 1.8, base_width = 6
)

