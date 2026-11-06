script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))
start_analysis("phenotypes", c("ggplot2", "dplyr", "readr", "tidyr", "rstatix", "multcompView", "cowplot"))
dir.create("pvalues", showWarnings = FALSE)
library(ggplot2)
library(dplyr)

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

rootlens <- readr::read_tsv(input_path("phenotypes", "root_lengths.tsv")) %>%
  tidyr::gather("Treatment", "Length", 2:6)
rootlens$Genotype <- factor(rootlens$Genotype, levels = c("WT", "aabb"))

rootlensStats <- rootlens %>%
  group_by(Genotype, Treatment) %>%
  summarise(
    n = n(),
    lower = quantile(Length, .25),
    mid = quantile(Length, .5),
    upper = quantile(Length, .75)
  )

test1 <- rstatix::dunn_test(rootlens[rootlens$Genotype == "WT", ], Length ~ Treatment)
multcompView::multcompLetters(structure(test1$p.adj, names=paste0(test1$group1,"-",test1$group2)))

test2 <- rstatix::dunn_test(rootlens[rootlens$Genotype == "aabb", ], Length ~ Treatment)
multcompView::multcompLetters(structure(test2$p.adj, names=paste0(test2$group1,"-",test2$group2)))

readr::write_tsv(test1,
  "pvalues/Rootlengths_treatments_WT.tsv")
readr::write_tsv(test2,
  "pvalues/Rootlengths_treatments_aabb.tsv")

p1 <- ggplot() +
  geom_point(aes(Genotype, Length, colour = Treatment, group = Treatment), data = rootlens,
    size = 0.75, position = position_jitterdodge(jitter.width = .08)) +
  geom_errorbar(aes(Genotype, ymin = lower, ymax = upper, group = Treatment),
    data = rootlensStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = mid, ymax = mid, group = Treatment),
    data = rootlensStats, width = 0.5, position = position_dodge(.75)) +
  geom_vline(xintercept = 1.5, lty = 2, col = "grey") +
  theme_bw() +
  scale_colour_manual(values = colorBlindBlack8,
    name = element_blank()) +
  ylab("Root length (cm)") +
  xlab(NULL) +
  ylim(c(0, 7.25)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")))) +
  theme(
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

seedyields <- readr::read_tsv(input_path("phenotypes", "seed_yields.tsv"))
seedyields1 <- seedyields[, 1:3]
seedyields2 <- seedyields[, 4:6]
colnames(seedyields1) <- c("Treatment", "Weight", "Count")
colnames(seedyields2) <- c("Treatment", "Weight", "Count")
seedyields1$Genotype <- "WT"
seedyields2$Genotype <- "aabb"
seedyields <- rbind(seedyields1, seedyields2)
seedyields$Treatment <- gsub("_\\d", "", seedyields$Treatment)
seedyields$Genotype <- factor(seedyields$Genotype, levels = c("WT", "aabb"))

seedyieldsStats <- seedyields %>%
  group_by(Genotype, Treatment) %>%
  summarise(
    n = n(),
    Weight.lower = quantile(Weight, .25),
    Weight.mid = quantile(Weight, .5),
    Weight.upper = quantile(Weight, .75),
    Count.lower = quantile(Count, .25),
    Count.mid = quantile(Count, .5),
    Count.upper = quantile(Count, .75)
  )

test3 <- rstatix::dunn_test(seedyields[seedyields$Genotype=="aabb", ], Weight ~ Treatment)
multcompView::multcompLetters(structure(test3$p.adj, names=paste0(test3$group1,"-",test3$group2)))
test4 <- rstatix::dunn_test(seedyields[seedyields$Genotype=="WT", ], Weight ~ Treatment)
multcompView::multcompLetters(structure(test4$p.adj, names=paste0(test4$group1,"-",test4$group2)))

readr::write_tsv(test3,
  "pvalues/SeedWeights_treatments_aabb.tsv")
readr::write_tsv(test4,
  "pvalues/SeedWeights_treatments_WT.tsv")

test5 <- rstatix::dunn_test(seedyields[seedyields$Genotype=="aabb", ], Count ~ Treatment)
multcompView::multcompLetters(structure(test5$p.adj, names=paste0(test5$group1,"-",test5$group2)))
test6 <- rstatix::dunn_test(seedyields[seedyields$Genotype=="WT", ], Count ~ Treatment)
multcompView::multcompLetters(structure(test6$p.adj, names=paste0(test6$group1,"-",test6$group2)))

readr::write_tsv(test5,
  "pvalues/SeedCounts_treatments_aabb.tsv")
readr::write_tsv(test6,
  "pvalues/SeedCounts_treatments_WT.tsv")

p2 <- ggplot() +
  geom_point(aes(Genotype, Weight, colour = Treatment, group = Treatment), data = seedyields,
    size = 0.75, position = position_jitterdodge(jitter.width = .08)) +
  geom_errorbar(aes(Genotype, ymin = Weight.lower, ymax = Weight.upper, group = Treatment),
    data = seedyieldsStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = Weight.mid, ymax = Weight.mid, group = Treatment),
    data = seedyieldsStats, width = 0.5, position = position_dodge(.75)) +
  geom_vline(xintercept = 1.5, lty = 2, col = "grey") +
  theme_bw() +
  scale_colour_manual(values = colorBlindBlack8,
    name = element_blank()) +
  ylab("Seed yield (mg)") +
  xlab(NULL) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")))) +
  ylim(c(0, 225)) +
  theme(
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

p3 <- ggplot() +
  geom_point(aes(Genotype, Count, colour = Treatment, group = Treatment), data = seedyields,
    size = 0.75, position = position_jitterdodge(jitter.width = .08)) +
  geom_errorbar(aes(Genotype, ymin = Count.lower, ymax = Count.upper, group = Treatment),
    data = seedyieldsStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = Count.mid, ymax = Count.mid, group = Treatment),
    data = seedyieldsStats, width = 0.5, position = position_dodge(.75)) +
  geom_vline(xintercept = 1.5, lty = 2, col = "grey") +
  theme_bw() +
  scale_colour_manual(values = colorBlindBlack8,
    name = element_blank()) +
  ylab("Total seed count") +
  xlab(NULL) +
  ylim(c(0, 9500)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")))) +
  theme(
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

g3fert <- readr::read_tsv(input_path("phenotypes", "development.tsv"))
g3fert$Genotype <- factor(g3fert$Genotype, levels = c("WT", "mtn1-1mtn2-1", "G3"))

g3fertStats <- g3fert %>%
  group_by(Genotype, Treatment) %>%
  summarise(
    DaysBolt.lower = quantile(`Days to Bolt`, 0.25),
    DaysBolt.mid = quantile(`Days to Bolt`, 0.5),
    DaysBolt.upper = quantile(`Days to Bolt`, 0.75),
    LeavesBolt.lower = quantile(`Leaves at Bolting`, 0.25),
    LeavesBolt.mid = quantile(`Leaves at Bolting`, 0.5),
    LeavesBolt.upper = quantile(`Leaves at Bolting`, 0.75),
    Size.lower = quantile(`Size of Plant`, 0.25),
    Size.mid = quantile(`Size of Plant`, 0.5),
    Size.upper = quantile(`Size of Plant`, 0.75),
    SiliqueNum.lower = quantile(`Number of Siliques`, 0.25),
    SiliqueNum.mid = quantile(`Number of Siliques`, 0.5),
    SiliqueNum.upper = quantile(`Number of Siliques`, 0.75),
    WeightSeeds.lower = quantile(`Seed weight`, 0.25, na.rm = TRUE),
    WeightSeeds.mid = quantile(`Seed weight`, 0.5, na.rm = TRUE),
    WeightSeeds.upper = quantile(`Seed weight`, 0.75, na.rm = TRUE),
    SeedsNum.lower = quantile(`Number of Seeds`, 0.25),
    SeedsNum.mid = quantile(`Number of Seeds`, 0.5),
    SeedsNum.upper = quantile(`Number of Seeds`, 0.75),
    BudNum.lower = quantile(`Buds cluster`, 0.25),
    BudNum.mid = quantile(`Buds cluster`, 0.5),
    BudNum.upper = quantile(`Buds cluster`, 0.75)
  )

g3fert$GenoTreat <- paste0(g3fert$Genotype, ":", g3fert$Treatment)

g <- function(x) gsub("-", "_", x, fixed = TRUE)

test7 <- rstatix::dunn_test(g3fert, `Days to Bolt` ~ GenoTreat)
multcompView::multcompLetters(structure(test7$p.adj, names=paste0(g(test7$group1),"-",g(test7$group2))))
test8 <- rstatix::dunn_test(g3fert, `Leaves at Bolting` ~ GenoTreat)
multcompView::multcompLetters(structure(test8$p.adj, names=paste0(g(test8$group1),"-",g(test8$group2))))
test9 <- rstatix::dunn_test(g3fert, `Size of Plant` ~ GenoTreat)
multcompView::multcompLetters(structure(test9$p.adj, names=paste0(g(test9$group1),"-",g(test9$group2))))
test10 <- rstatix::dunn_test(g3fert, `Number of Siliques` ~ GenoTreat)
multcompView::multcompLetters(structure(test10$p.adj, names=paste0(g(test10$group1),"-",g(test10$group2))))
test11 <- rstatix::dunn_test(g3fert, `Seed weight` ~ GenoTreat)
multcompView::multcompLetters(structure(test11$p.adj, names=paste0(g(test11$group1),"-",g(test11$group2))))
test12 <- rstatix::dunn_test(g3fert, `Number of Seeds` ~ GenoTreat)
multcompView::multcompLetters(structure(test12$p.adj, names=paste0(g(test12$group1),"-",g(test12$group2))))
test13 <- rstatix::dunn_test(g3fert, `Buds cluster` ~ GenoTreat)
multcompView::multcompLetters(structure(test13$p.adj, names=paste0(g(test13$group1),"-",g(test13$group2))))

readr::write_tsv(test7,
  "pvalues/G3dev_DaysBolt.tsv")
readr::write_tsv(test8,
  "pvalues/G3dev_LeavesBolt.tsv")
readr::write_tsv(test9,
  "pvalues/G3dev_SizePlant.tsv")
readr::write_tsv(test10,
  "pvalues/G3dev_NumSiliques.tsv")
readr::write_tsv(test11,
  "pvalues/G3dev_SeedWeight.tsv")
readr::write_tsv(test12,
  "pvalues/G3dev_SeedNum.tsv")
readr::write_tsv(test13,
  "pvalues/G3dev_BudsNum.tsv")

g3fert$Treatment <- factor(g3fert$Treatment, levels = c("MS", "Spd"))

f1 <- ggplot() +
  geom_point(aes(Genotype, `Number of Seeds`, colour = Treatment, group = Treatment),
    data = g3fert, size = 0.75, position = position_jitterdodge(jitter.width = .1)) +
  geom_errorbar(aes(Genotype, ymin = SeedsNum.lower, ymax = SeedsNum.upper, group = Treatment),
    data = g3fertStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = SeedsNum.mid, ymax = SeedsNum.mid, group = Treatment),
    data = g3fertStats, width = 0.5, position = position_dodge(.75)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "brown"), limits = c("Spd", "MS"),
    labels = c("-Spd", "+Spd"),
    name = element_blank()) +
  ylab("Total seed count") +
  xlab(NULL) +
  ylim(c(0, 18000)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3")) +
  theme(
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

f2 <- ggplot() +
  geom_point(aes(Genotype, `Seed weight`, colour = Treatment, group = Treatment),
    data = g3fert, size = 0.75, position = position_jitterdodge(jitter.width = .1)) +
  geom_errorbar(aes(Genotype, ymin = WeightSeeds.lower, ymax = WeightSeeds.upper, group = Treatment),
    data = g3fertStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = WeightSeeds.mid, ymax = WeightSeeds.mid, group = Treatment),
    data = g3fertStats, width = 0.5, position = position_dodge(.75)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "brown"),#limits = c("Spd", "MS"),
    labels = c("-Spd", "+Spd"),
    name = element_blank()) +
  ylab("Seed yield (mg)") +
  xlab(NULL) +
  ylim(c(0, 450)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3")) +
  theme(
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

f3 <- ggplot() +
  geom_point(aes(Genotype, `Days to Bolt`, colour = Treatment, group = Treatment),
    data = g3fert, size = 0.75, position = position_jitterdodge(jitter.width = .1)) +
  geom_errorbar(aes(Genotype, ymin = DaysBolt.lower, ymax = DaysBolt.upper, group = Treatment),
    data = g3fertStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = DaysBolt.mid, ymax = DaysBolt.mid, group = Treatment),
    data = g3fertStats, width = 0.5, position = position_dodge(.75)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "brown"),#limits = c("Spd", "MS"),
    labels = c("-Spd", "+Spd"),
    name = element_blank()) +
  ylab("Days to bolt") +
  xlab(NULL) +
  ylim(c(0, 60)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3")) +
  theme(
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

f4 <- ggplot() +
  geom_point(aes(Genotype, `Leaves at Bolting`, colour = Treatment, group = Treatment),
    data = g3fert, size = 0.75, position = position_jitterdodge(jitter.width = .1)) +
  geom_errorbar(aes(Genotype, ymin = LeavesBolt.lower, ymax = LeavesBolt.upper, group = Treatment),
    data = g3fertStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = LeavesBolt.mid, ymax = LeavesBolt.mid, group = Treatment),
    data = g3fertStats, width = 0.5, position = position_dodge(.75)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "brown"),#limits = c("Spd", "MS"),
    labels = c("-Spd", "+Spd"),
    name = element_blank()) +
  ylab("Leaves at bolting") +
  xlab(NULL) +
  ylim(c(0, 39)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3")) +
  theme(
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

f5 <- ggplot() +
  geom_point(aes(Genotype, `Size of Plant`, colour = Treatment, group = Treatment),
    data = g3fert, size = 0.75, position = position_jitterdodge(jitter.width = .1)) +
  geom_errorbar(aes(Genotype, ymin = Size.lower, ymax = Size.upper, group = Treatment),
    data = g3fertStats, width = 0.4, position = position_dodge(.75)) +
  geom_errorbar(aes(Genotype, ymin = Size.mid, ymax = Size.mid, group = Treatment),
    data = g3fertStats, width = 0.5, position = position_dodge(.75)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "brown"),#limits = c("Spd", "MS"),
    labels = c("-Spd", "+Spd"),
    name = element_blank()) +
  ylab("Plant height (cm)") +
  xlab(NULL) +
  ylim(c(0, 50)) +
  scale_x_discrete(labels = c("WT", bquote(italic("mtn1-1\nmtn2-1")), "G3")) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(size = 0.75),
    rect = element_rect(size = 0.9),
    axis.ticks = element_line(size = 0.75),
    text = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(colour = "black", size = 6),
    axis.text.x = element_text(colour = "black", size = 6),
    panel.grid = element_blank()
  ) + coord_flip()

cowplot::save_plot(
  "phenotypes.pdf",
  cowplot::plot_grid(ncol = 2, byrow = FALSE, align = "vh", axis = "tblr",
    f1, f2, f3, f4, f5,
    p1, p2, p3, NULL, NULL
  ),
  base_height = 6, base_width = 5
)

