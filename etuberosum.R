library(ggtree)
library(dplyr)

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


# Load data 
tree <- read.tree("results/resistify/etuberosum/nbarc.tree")
data <- read.table("results/resistify/etuberosum/results.tsv", header=TRUE, sep="\t")
bed <- read.table("results/bed/etuberosum.bed", header=FALSE, sep="\t")

# join by Sequence and V1
joined <- data %>%
    left_join(bed, by = c("Sequence" = "V4")) %>%
    filter(MADA == "True") %>%
    filter(V1 == "chr04") %>%
    filter(V2 < 2000000)

# add "_1" to the end of data$Sequence
data$Sequence <- paste(data$Sequence, "_1", sep="")

# create df of MADA column for heatmaps
heatmap <- data.frame(data$Sequence, data$MADA, data$CJID)

p <- ggtree(tree, layout = "circular") %<+% data +
    geom_tiplab(aes(alpha = MADA), align = TRUE, size = 2) +
    geom_tippoint(aes(color = Classification)) +
    scale_color_manual(values = paired_colours)

ggsave("tree.png", p, width = 12, height = 12, units = "in")

ggplot(joined, aes(xmin = V2, xmax = V3, y = V1, fill = Classification)) +
    geom_gene_arrow()
