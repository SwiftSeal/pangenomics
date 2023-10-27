library(ggtree)
library(gggenes)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)

nice_names <- c(
    "candolleanum" = "S. candolleanum",
    "etuberosum" = "S. etuberosum",
    "lycopersicum" = "S. lycopersicum",
    "tuberosum_phureja_E4_63" = "S. tuberosum\n(phureja E4-63)",
    "tuberosum_tuberosum_RH10_15" = "S. tuberosum\n(tuberosum RH10-15)",
    "chilense" = "S. chilense",
    "galapagense" = "S. galapagense",
    "neorickii" = "S. neorickii",
    "tuberosum_phureja_E86_69" = "S. tuberosum\n(phureja E86-69)",
    "tuberosum_tuberosum_RH" = "S. tuberosum\n(tuberosum RH)",
    "chmielewskii" = "S. chmielewskii",
    "habrochaites" = "S. habrochaites",
    "peruvianum" = "S. peruvianum",
    "tuberosum_stenotomum_A6_26" = "S. tuberosum\n(stenotomum A6-26)",
    "verrucosum" = "S. verrucosum",
    "corneliomulleri" = "S. corneliomulleri",
    "lycopersicoides" = "S. lycopersicoides",
    "pimpinellifolium" = "S. pimpinellifolium",
    "tuberosum_stenotomum_PG6359" = "S. tuberosum\n(stenotomum PG6359)"
)

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

# of the CNLs, how many are True for MADA?
resistify %>%
    filter(Classification == "CNL") %>%
    summarise(count = sum(MADA == "True"), total = n())

# create column of nice names
resistify_summary$GenomeFormatted <- nice_names[match(resistify_summary$Genome, names(nice_names))]

# create tree
tree_plot <- ggtree(tree) +
    geom_tiplab(align = TRUE, color = "white") +
    theme(plot.margin = unit(c(0,0,0.35,0), "in"))
order <- get_taxa_name(tree_plot)
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
    geom_text(aes(x = 0, y = 0, label = GenomeFormatted, fontface = "italic"), size = 2, hjust = 0) +
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
                         rel_widths = c(1, 0.7, 3))

ggsave("results/plot.png", merged_plot, width = 6, height = 6, units = "in")


