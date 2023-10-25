library(ggtree)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)

# colours
# plot styling
paired_colours = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00")
# load the tree
tree <- read.tree("results/results/SpeciesTree_rooted.txt")

# load the resistify results
genomes <- list.dirs("results/resistify", recursive = FALSE)
resistify <- map(genomes, ~ read.table(paste0(.x, "/results.tsv"), header = TRUE) %>%
    mutate(Genome = basename(.x))
)
resistify <- bind_rows(resistify)
resistify_summary <- resistify %>%
    group_by(Genome, Classification) %>%
    summarise(count = n())

# create tree
tree_plot <- ggtree(tree) +
    geom_tiplab(align = TRUE, color = "white") +
    theme(plot.margin = unit(c(0,0,0.35,0), "in"))
order <- get_taxa_name(plot)
get_taxa_name(plot)
resistify_summary$Genome <- factor(resistify_summary$Genome, levels = order)

# plot facet of resistify summaries
resistify_plot <- ggplot(resistify_summary %>% mutate(dummy = "dummy"), aes(y = dummy, x = count, fill = Classification)) +
    geom_bar(stat = "identity") +
    facet_grid(Genome ~ .) +
    theme(legend.position = c(0.75, 0.8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()
    ) +
    scale_fill_manual(values = paired_colours)

# plot a facet of just the genome names
names_plot <- ggplot(resistify_summary) +
    geom_text(aes(x = 0, y = 0, label = Genome), size = 2, hjust = 0) +
    xlim(0,1) +
    facet_grid(Genome ~ .) +
    theme(plot.margin = unit(c(0.06,0,0.4,0), "in"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()
    )

merged_plot <- plot_grid(tree_plot,
                         names_plot,
                         resistify_plot,
                         ncol = 3,
                         rel_widths = c(2, 1, 3))

ggsave("results/plot.png", merged_plot, width = 6, height = 6, units = "in")

