# Meristem, embryo and internode measurements from the Zenodo deposit.

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))

start_analysis("development", c("ggplot2", "dplyr", "readr", "rstatix"))

library(ggplot2)
library(dplyr)

measurements <- list(
  meristem = c("Meristem sizes.csv", "Size", "Meristem length (micrometres)"),
  columella = c("Embryo columella cells.csv", "Columella Cells", "Columella cells"),
  divisions = c("Embryo divisions.csv", "Division", "QC cell divisions")
)

for (name in names(measurements)) {

  specification <- measurements[[name]]
  x <- readr::read_csv(input_path("zenodo", specification[1]), show_col_types = FALSE)

  x$Genotype[x$Genotype == "dm"] <- "mtn1-1mtn2-1"
  x$Genotype <- factor(x$Genotype, levels = c("WT", "mtn1-1mtn2-1", "G3"))
  x$Treatment <- factor(x$Treatment, levels = c("MS", "Spd"))
  x$Value <- x[[specification[2]]]
  x$Group <- paste(x$Genotype, x$Treatment, sep = ":")

  if (anyNA(x[, c("Genotype", "Treatment", "Value")])) {
    stop("Missing measurements or unknown sample labels in ", specification[1])
  }

  test <- rstatix::dunn_test(x, Value ~ Group, p.adjust.method = "BH")
  readr::write_tsv(test, paste0(name, "_dunn.tsv"))
  readr::write_tsv(rstatix::kruskal_test(x, Value ~ Group),
    paste0(name, "_kruskal.tsv"))

  summary <- x %>%
    group_by(Genotype, Treatment) %>%
    summarise(
      n = n(),
      lower = quantile(Value, .25),
      median = median(Value),
      upper = quantile(Value, .75),
      .groups = "drop"
    )

  readr::write_tsv(summary, paste0(name, "_summary.tsv"))

  plot <- ggplot(x, aes(Genotype, Value, colour = Treatment, group = Treatment)) +
    geom_point(position = position_jitterdodge(jitter.width = .1), size = 1) +
    geom_errorbar(aes(Genotype, ymin = lower, ymax = upper, group = Treatment),
      data = summary, inherit.aes = FALSE,
      width = .3, position = position_dodge(.75), colour = "black") +
    geom_errorbar(aes(Genotype, ymin = median, ymax = median, group = Treatment),
      data = summary, inherit.aes = FALSE,
      width = .4, position = position_dodge(.75), colour = "black") +
    scale_colour_manual(values = c("MS" = "black", "Spd" = "brown"),
      labels = c("-Spd", "+Spd")) +
    labs(x = NULL, y = specification[3], colour = NULL) +
    theme_bw()

  ggsave(paste0(name, ".pdf"), plot, width = 5, height = 3.5)

}

# The deposited long table contains the five representative plants per group.
# It does not contain all 15 plants per group used for the published PCA.

internodes <- readr::read_csv(input_path("zenodo", "Internode measurements.csv"),
  col_select = c(Genotype, Treatment, Plant, Internode, Length),
  show_col_types = FALSE, name_repair = "minimal")

internodes$Genotype[internodes$Genotype == "dm"] <- "mtn1-1mtn2-1"
internodes$Genotype <- factor(internodes$Genotype,
  levels = c("WT", "mtn1-1mtn2-1", "G3"))
internodes$Treatment <- factor(internodes$Treatment, levels = c("MS", "Spd"))

key <- with(internodes, paste(Genotype, Treatment, Plant, Internode))

if (anyNA(internodes) || anyDuplicated(key)) {
  stop("Missing or duplicated internode measurements")
}

plot <- ggplot(internodes, aes(Internode, Length, group = Plant, colour = Plant)) +
  geom_line() +
  geom_point(size = .8) +
  facet_grid(Genotype ~ Treatment) +
  scale_x_continuous(breaks = 1:10) +
  labs(x = "Internode", y = "Length (cm)", colour = "Plant") +
  theme_bw()

ggsave("internode_profiles.pdf", plot, width = 7, height = 6)
readr::write_tsv(internodes, "internode_measurements.tsv")
