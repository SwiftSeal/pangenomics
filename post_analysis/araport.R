library(tidyverse)
library(ggtree)
library(ggpubr)

paired_colours = c(
  "CN" = "#a6cee3", 
  "CNL" = "#1f78b4",
  "N" = "#b2df8a",
  "NL" = "#33a02c",
  "RN" = "#fb9a99",
  "RNL" = "#e31a1c",
  "TN" = "#fdbf6f",
  "TNL" = "#ff7f00"
)

results <- read.table(
  "araport11_benchmark/resistify_output/results.tsv",
  header=TRUE,
  sep="\t"
) |>
  mutate(Sequence = paste0(Sequence, "_1"))

tree <- read.tree("araport11_benchmark/resistify_output/nbarc.tree")
tree <- treeio::root(tree, outgroup = "ced4")

plot_tree <- ggtree(tree, layout = "circular") %<+% results

plot_tree <- plot_tree +
  geom_tippoint(aes(colour = Classification), size = 1) +
  scale_colour_manual(values = paired_colours)


plot_classification <- ggplot(results, aes(x = Classification, fill = Classification)) +
  geom_bar() +
  scale_fill_manual(values = paired_colours) +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "Classification", y = "# NLRs")

plot_motifs <- ggplot(results, aes(x = as.factor(NBARC_motifs), fill = Classification)) +
  geom_bar() +
  scale_fill_manual(values = paired_colours) +
  theme_bw(base_size = 8) + theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = "# NB-ARC motifs", y = "# NLRs")

plot_merged <- gridExtra::grid.arrange(
  plot_classification,
  plot_motifs,
  plot_tree,
  ncol = 2,
  nrow = 2,
  layout_matrix = rbind(c(1, 3), c(2, 3))
)

plot_merged <- as_ggplot(plot_merged) +
  cowplot::draw_plot_label(label = c("A", "B", "C"), size = 8, x = c(0, 0, 0.5), y = c(1, 0.5, 1))

ggsave("araport11.png", plot_merged, width = 5.9, height = 3, dpi = 600)


summary_table <- results %>%
  group_by(Classification) %>%
  # calculate percentage which has TRUE in MADA column
  summarise(count = n(), 
            percentage_mada = 100 * sum(MADA == "True") / count,
            percentage_CJID = 100 * sum(CJID == "True") / count,
            percentage_functional = 100 * sum(Functionality == "Likely") / count)

write.table(summary_table, "results/summary_table.tsv", sep="\t", quote=FALSE, row.names=FALSE)

motif_translation = c(
  "extEDVID" = "CC",
  "bA" = "TIR",
  "aA" = "TIR",
  "bC" = "TIR",
  "aC" = "TIR",
  "bDaD1" = "TIR",
  "aD3" = "TIR",
  "VG" = "NB-ARC",
  "P-loop" = "NB-ARC",
  "RNSB-A" = "NB-ARC",
  "Walker-B" = "NB-ARC",
  "RNSB-B" = "NB-ARC",
  "RNSB-C" = "NB-ARC",
  "RNSB-D" = "NB-ARC",
  "GLPL" = "NB-ARC",
  "MHD" = "NB-ARC",
  "LxxLxL" = "LRR"
)

results <- read_tsv("araport11_benchmark/resistify_output/results.tsv")
subgroup <- sample(results$Sequence, 15)
results <- results |>
  filter(Sequence %in% subgroup)
domains <- read_tsv("araport11_benchmark/resistify_output/domains.tsv") |>
  filter(Sequence %in% subgroup)
motifs <- read_tsv("araport11_benchmark/resistify_output/motifs.tsv") |>
  mutate(Domain = motif_translation[Motif]) |>
  filter(Sequence %in% subgroup)

myplot <- ggplot() +
  geom_segment(data = results, aes(y = Sequence, yend = Sequence, x = 0, xend = Length)) +
  geom_segment(data = domains, aes(y = Sequence, yend = Sequence, x = Start, xend = End, colour = Domain)) +
  geom_point(data = motifs, aes(y = Sequence, x = Position, colour = Domain), size = 1) +
  scale_colour_tableau() +
  theme_bw() + theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)

ggsave("testplot.png", myplot, width = 5.9, height = 3, dpi = 600)
