script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("phenotypes_supplement", c("ggplot2", "dplyr", "readr", "tidyr", "rstatix", "multcompView", "cowplot"))
dir.create("pvalues", showWarnings = FALSE)
library(ggplot2)
library(dplyr)

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- function(x) gsub("-", "_", x, fixed = TRUE)

mtn <- readr::read_csv(input_path("zenodo", "MTN enzyme activity.txt"), comment = "#")
mtnStats <- mtn %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    lower = quantile(activity, .25),
    mid = quantile(activity, .5),
    upper = quantile(activity, .75),
    mean = mean(activity),
    stdev = sd(activity)
  )

test1 <- rstatix::t_test(mtn, activity ~ genotype)
multcompView::multcompLetters(structure(test1$p.adj, names=paste0(g(test1$group1),"-",g(test1$group2))))

readr::write_tsv(test1,
  "pvalues/MTN_assay.tsv")

mtnStats$genotype <- factor(mtnStats$genotype, levels = c("wt", "m11m21", "m11m25"))

s2p1 <- ggplot() +
  geom_errorbar(aes(genotype, ymin = mean - stdev, ymax = mean + stdev), data = mtnStats,
    width = .55) +
  geom_bar(aes(genotype, mean, fill = genotype), data = mtnStats,
    stat = "identity", width = .8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5)) +
  scale_fill_manual(values = colorBlindBlack8) +
  scale_x_discrete(
    labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), bquote(italic("mtn1-1\nmtn2-5")))) +
  xlab(NULL) +
  ylab(expression(paste("MTN activity (nmol ", mg^{-1}, "prot ", min^{-1}, ")"))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_line(size = 0.75),
    rect = element_rect(size = 0.9),
    axis.ticks = element_line(size = 0.75),
    text = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 6),
    axis.text.x = element_text(colour = "black", size = 6),
    panel.grid = element_blank()
  ) +
  coord_flip()

rls <- readr::read_csv(input_path("zenodo", "Root lengths (MTA).csv"))
rls <- rls[rls$genotype %in% c("aabb", "aacc", "WT"), ]
rls <- rls[rls$treatment %in% c("MgSO4", "100 MTA", "250 MTA", "500 MTA"), ]

rlsStats <- rls %>%
  group_by(genotype, treatment) %>%
  summarise(
    n = n(),
    lower = quantile(Length, .25),
    mid = quantile(Length, .5),
    upper = quantile(Length, .75),
    mean = mean(Length),
    stdev = sd(Length),
    stderr = sd(Length) / sqrt(n())
  )

rls$GenoTreat <- paste0(rls$genotype, ":", rls$treatment)
test2 <- rstatix::dunn_test(rls, Length ~ GenoTreat)
multcompView::multcompLetters(structure(test2$p.adj, names=paste0(g(test2$group1),"-",g(test2$group2))))

readr::write_tsv(test2,
  "pvalues/MTA_rootlengths.tsv"
)

rlsStats$mta <- c("MgSO4" = 0, "100 MTA" = 100, "250 MTA" = 250, "500 MTA" = 500)[rlsStats$treatment]
rlsStats$genotype <- factor(rlsStats$genotype, levels = c("WT", "aabb", "aacc"))

s2p2 <- ggplot() +
  geom_line(aes(mta, mean, group = genotype), data = rlsStats) +
  geom_errorbar(aes(mta, ymin = mean - stderr, ymax = mean + stderr, group = genotype),
    data = rlsStats, width = 15) +
  geom_point(aes(mta, mean, fill = genotype), data = rlsStats, pch = 21, size = 1.) +
  scale_fill_manual(values = colorBlindBlack8,
    labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), bquote(italic("mtn1-1\nmtn2-5"))),
    name = element_blank()) +
  scale_x_continuous(breaks = c(0, 100, 250, 500)) +
  ylim(c(0, 2.75)) +
  ylab("Root length (cm)") +
  xlab(expression(paste("MTA (", "\U003BC", "M)"))) +
  theme_bw() +
  theme(
    rect = element_rect(size = 0.9),
    axis.ticks = element_line(size = 0.75),
    text = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 6),
    axis.text.x = element_text(colour = "black", size = 6),
    panel.grid = element_blank()
  )

cowplot::save_plot("mtn_activity_and_mta_roots.pdf", cowplot::plot_grid(s2p1, s2p2, nrow = 1), base_width = 6, base_height = 2)
