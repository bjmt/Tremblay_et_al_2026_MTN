# Plot the GO term P-values in Supplementary Dataset S4.

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_file)
script_folder <- dirname(normalizePath(script_file))

source(file.path(script_folder, "common.R"))

start_analysis("rnaseq_go", c("ggplot2", "readr"))

library(ggplot2)

go <- readr::read_tsv(input_path("rnaseq", "go_terms.tsv"), show_col_types = FALSE)

terms <- c("pollination", "extracellular region", "regulation of fertilization",
  "response to red or far red light", "cell wall organization",
  "pectinesterase activity", "response to abiotic stimulus",
  "response to light stimulus", "cellular response to hypoxia", "photosystem I")

go <- go[go$Term %in% terms, ]
go$Term <- factor(go$Term, levels = rev(sort(terms)))

ggplot(go, aes(Comparison, Term, size = -log10(Pvalue))) +
  geom_point() +
  labs(x = NULL, y = NULL, size = expression(-log[10](P))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("go_terms.pdf", width = 13, height = 7)
